# BioGPU Streaming Service
# Core streaming functionality for real-time updates

import asyncio
import json
import logging
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, asdict
import redis.asyncio as aioredis
from kafka import KafkaProducer
from kafka.errors import KafkaError
import websockets
from datetime import datetime
import jwt

logger = logging.getLogger(__name__)

class AuthenticationError(Exception):
    """Raised when authentication fails"""
    pass

class BioGPUStreamingService:
    """Core streaming service for publishing real-time updates"""
    
    def __init__(self, redis_url: str, kafka_servers: List[str]):
        self.redis_url = redis_url
        self.kafka_servers = kafka_servers
        self.redis_client = None
        self.kafka_producer = None
        self._initialized = False
        
    async def initialize(self):
        """Initialize Redis and Kafka connections"""
        try:
            # Initialize Redis with new API
            self.redis_client = await aioredis.from_url(
                self.redis_url,
                decode_responses=True
            )
            
            # Initialize Kafka producer
            self.kafka_producer = KafkaProducer(
                bootstrap_servers=self.kafka_servers,
                value_serializer=lambda v: json.dumps(v).encode('utf-8'),
                acks='all',
                retries=3,
                max_in_flight_requests_per_connection=5
            )
            
            self._initialized = True
            logger.info("BioGPU Streaming Service initialized successfully")
            
        except Exception as e:
            logger.error(f"Failed to initialize streaming service: {str(e)}")
            raise
    
    async def publish_update(self, update):
        """Publish an update to all streaming channels"""
        if not self._initialized:
            raise RuntimeError("Streaming service not initialized")
        
        try:
            # Convert update to dict
            update_dict = asdict(update) if hasattr(update, '__dataclass_fields__') else update.__dict__
            update_dict['timestamp'] = update_dict.get('timestamp') or datetime.utcnow().isoformat()
            
            # Publish to Redis for WebSocket subscribers
            channel = f"job_updates:{update.job_id}"
            await self.redis_client.publish(channel, json.dumps(update_dict))
            
            # Also publish to a general updates channel
            await self.redis_client.publish("all_job_updates", json.dumps(update_dict))
            
            # Store latest status in Redis
            status_key = f"job_status:{update.job_id}"
            await self.redis_client.setex(
                status_key, 
                3600,  # 1 hour TTL
                json.dumps(update_dict)
            )
            
            # Send to Kafka for persistent storage/processing
            self.kafka_producer.send('job_updates', update_dict)
            
            logger.debug(f"Published update for job {update.job_id}: {update.stage}")
            
        except Exception as e:
            logger.error(f"Failed to publish update: {str(e)}")
            raise
    
    async def get_job_status(self, job_id: str) -> Optional[Dict]:
        """Get the latest status of a job"""
        if not self._initialized:
            raise RuntimeError("Streaming service not initialized")
        
        status_key = f"job_status:{job_id}"
        status_json = await self.redis_client.get(status_key)
        
        if status_json:
            return json.loads(status_json)
        return None
    
    async def subscribe_to_job(self, job_id: str):
        """Subscribe to updates for a specific job"""
        if not self._initialized:
            raise RuntimeError("Streaming service not initialized")
        
        # Create a new pubsub instance for this subscription
        pubsub = self.redis_client.pubsub()
        channel = f"job_updates:{job_id}"
        await pubsub.subscribe(channel)
        return pubsub
    
    async def close(self):
        """Clean up connections"""
        if self.redis_client:
            await self.redis_client.close()
        
        if self.kafka_producer:
            self.kafka_producer.close()
        
        self._initialized = False
        logger.info("Streaming service closed")


class SecureWebSocketHandler:
    """Handles secure WebSocket connections with JWT authentication"""
    
    def __init__(self, config, streaming_service: BioGPUStreamingService):
        self.config = config
        self.streaming_service = streaming_service
        self.connected_clients: Dict[str, Set[websockets.WebSocketServerProtocol]] = {}
        
    def verify_jwt_token(self, token: str) -> Dict:
        """Verify JWT token and return payload"""
        try:
            payload = jwt.decode(
                token,
                self.config.jwt_secret,
                algorithms=['HS256']
            )
            
            # Check expiration
            if 'exp' in payload:
                exp = datetime.fromtimestamp(payload['exp'])
                if exp < datetime.utcnow():
                    raise AuthenticationError("Token has expired")
            
            return payload
            
        except jwt.ExpiredSignatureError:
            raise AuthenticationError("Token has expired")
        except jwt.InvalidTokenError as e:
            raise AuthenticationError(f"Invalid token: {str(e)}")
    
    async def handle_connection(self, websocket, path):
        """Handle incoming WebSocket connection"""
        client_id = None
        subscribed_jobs = set()
        
        try:
            # Authentication phase
            auth_message = await asyncio.wait_for(websocket.recv(), timeout=5.0)
            auth_data = json.loads(auth_message)
            
            if auth_data.get('type') != 'auth':
                await websocket.send(json.dumps({
                    'type': 'error',
                    'message': 'First message must be authentication'
                }))
                return
            
            # Verify JWT token
            token = auth_data.get('token')
            if not token:
                raise AuthenticationError("No token provided")
            
            user_info = self.verify_jwt_token(token)
            client_id = user_info.get('user_id')
            
            # Send auth success
            await websocket.send(json.dumps({
                'type': 'auth_success',
                'user_id': client_id
            }))
            
            # Add to connected clients
            if client_id not in self.connected_clients:
                self.connected_clients[client_id] = set()
            self.connected_clients[client_id].add(websocket)
            
            logger.info(f"Client {client_id} connected via WebSocket")
            
            # Handle client messages
            async for message in websocket:
                try:
                    data = json.loads(message)
                    
                    if data.get('type') == 'subscribe':
                        job_id = data.get('job_id')
                        if job_id:
                            subscribed_jobs.add(job_id)
                            
                            # Send current status
                            status = await self.streaming_service.get_job_status(job_id)
                            if status:
                                await websocket.send(json.dumps({
                                    'type': 'status_update',
                                    'data': status
                                }))
                            
                            # Start forwarding updates
                            asyncio.create_task(
                                self._forward_job_updates(websocket, job_id)
                            )
                    
                    elif data.get('type') == 'unsubscribe':
                        job_id = data.get('job_id')
                        if job_id in subscribed_jobs:
                            subscribed_jobs.remove(job_id)
                    
                    elif data.get('type') == 'ping':
                        await websocket.send(json.dumps({'type': 'pong'}))
                        
                except json.JSONDecodeError:
                    await websocket.send(json.dumps({
                        'type': 'error',
                        'message': 'Invalid JSON'
                    }))
                except Exception as e:
                    logger.error(f"Error handling message: {str(e)}")
                    await websocket.send(json.dumps({
                        'type': 'error',
                        'message': str(e)
                    }))
        
        except AuthenticationError as e:
            await websocket.send(json.dumps({
                'type': 'auth_error',
                'message': str(e)
            }))
        
        except asyncio.TimeoutError:
            await websocket.send(json.dumps({
                'type': 'error',
                'message': 'Authentication timeout'
            }))
        
        except websockets.exceptions.ConnectionClosed:
            logger.info(f"Client {client_id} disconnected")
        
        except Exception as e:
            logger.error(f"WebSocket error: {str(e)}")
            await websocket.send(json.dumps({
                'type': 'error',
                'message': 'Internal server error'
            }))
        
        finally:
            # Clean up
            if client_id and client_id in self.connected_clients:
                self.connected_clients[client_id].discard(websocket)
                if not self.connected_clients[client_id]:
                    del self.connected_clients[client_id]
    
    async def _forward_job_updates(self, websocket, job_id: str):
        """Forward job updates to WebSocket client"""
        try:
            pubsub = await self.streaming_service.subscribe_to_job(job_id)
            
            async for message in pubsub.listen():
                if websocket.closed:
                    break
                
                # Skip subscription confirmation messages
                if message['type'] == 'message':
                    update = json.loads(message['data'])
                    await websocket.send(json.dumps({
                        'type': 'job_update',
                        'data': update
                    }))
                
        except Exception as e:
            logger.error(f"Error forwarding updates: {str(e)}")
    
    async def broadcast_to_user(self, user_id: str, message: Dict):
        """Broadcast message to all connections for a user"""
        if user_id in self.connected_clients:
            message_json = json.dumps(message)
            
            # Send to all user's connections
            disconnected = []
            for websocket in self.connected_clients[user_id]:
                try:
                    await websocket.send(message_json)
                except websockets.exceptions.ConnectionClosed:
                    disconnected.append(websocket)
            
            # Clean up disconnected clients
            for ws in disconnected:
                self.connected_clients[user_id].discard(ws)


async def start_websocket_server(config, streaming_service: BioGPUStreamingService):
    """Start the WebSocket server"""
    handler = SecureWebSocketHandler(config, streaming_service)
    
    server = await websockets.serve(
        handler.handle_connection,
        "0.0.0.0",
        8765,
        ping_interval=30,
        ping_timeout=10
    )
    
    logger.info("WebSocket server started on port 8765")
    return server