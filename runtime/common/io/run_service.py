#!/usr/bin/env python3
# Main entry point for BioGPU Streaming Service

import asyncio
import signal
import sys
import logging
import argparse
from pathlib import Path

# Import our service modules
from biogpu_service import BioGPUConfig, ProductionJobDispatcher
from biogpu_streaming import BioGPUStreamingService, start_websocket_server

# Configure logging
def setup_logging(log_level: str = 'INFO'):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('biogpu_service.log')
        ]
    )
    
    # Set specific log levels for noisy libraries
    logging.getLogger('kafka').setLevel(logging.WARNING)
    logging.getLogger('websockets').setLevel(logging.WARNING)

class ServiceManager:
    """Manages all service components"""
    
    def __init__(self, config: BioGPUConfig):
        self.config = config
        self.dispatcher = None
        self.streaming_service = None
        self.websocket_server = None
        self.running = False
        
    async def start(self):
        """Start all services"""
        logging.info("Starting BioGPU Service Manager...")
        
        try:
            # Initialize streaming service
            self.streaming_service = BioGPUStreamingService(
                self.config.redis_url,
                self.config.kafka_servers
            )
            await self.streaming_service.initialize()
            
            # Start WebSocket server
            self.websocket_server = await start_websocket_server(
                self.config,
                self.streaming_service
            )
            
            # Initialize and start job dispatcher
            self.dispatcher = ProductionJobDispatcher(self.config)
            await self.dispatcher.initialize()
            
            self.running = True
            logging.info("All services started successfully")
            
            # Start consuming jobs
            await self.dispatcher.start_consuming()
            
        except Exception as e:
            logging.error(f"Failed to start services: {str(e)}")
            await self.shutdown()
            raise
    
    async def shutdown(self):
        """Gracefully shutdown all services"""
        logging.info("Shutting down services...")
        self.running = False
        
        # Close WebSocket server
        if self.websocket_server:
            self.websocket_server.close()
            await self.websocket_server.wait_closed()
        
        # Close streaming service
        if self.streaming_service:
            await self.streaming_service.close()
        
        # Close dispatcher (would need to implement this)
        # if self.dispatcher:
        #     await self.dispatcher.close()
        
        logging.info("All services shut down")

async def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='BioGPU Streaming Service')
    parser.add_argument(
        '--config',
        type=str,
        help='Path to configuration file (JSON format)',
        default=None
    )
    parser.add_argument(
        '--log-level',
        type=str,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging level'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Validate configuration without starting services'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_level)
    
    # Load configuration
    config = BioGPUConfig()
    
    if args.config:
        # Load from config file if provided
        import json
        config_path = Path(args.config)
        if config_path.exists():
            with open(config_path) as f:
                config_data = json.load(f)
                for key, value in config_data.items():
                    if hasattr(config, key):
                        setattr(config, key, value)
            logging.info(f"Loaded configuration from {config_path}")
        else:
            logging.error(f"Configuration file not found: {config_path}")
            sys.exit(1)
    
    # Log configuration
    logging.info("Configuration:")
    for field in config.__dataclass_fields__:
        value = getattr(config, field)
        # Mask sensitive values
        if 'secret' in field.lower() or 'password' in field.lower():
            value = '***'
        logging.info(f"  {field}: {value}")
    
    if args.dry_run:
        logging.info("Dry run complete, configuration is valid")
        return
    
    # Create service manager
    service_manager = ServiceManager(config)
    
    # Set up signal handlers
    def signal_handler(sig, frame):
        logging.info(f"Received signal {sig}")
        asyncio.create_task(service_manager.shutdown())
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    try:
        # Start services
        await service_manager.start()
    except KeyboardInterrupt:
        logging.info("Keyboard interrupt received")
    except Exception as e:
        logging.error(f"Service error: {str(e)}")
        sys.exit(1)
    finally:
        await service_manager.shutdown()

if __name__ == "__main__":
    try:
        asyncio.run(main())
    except Exception as e:
        logging.critical(f"Fatal error: {str(e)}")
        sys.exit(1)