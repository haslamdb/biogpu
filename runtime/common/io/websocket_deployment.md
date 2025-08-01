# BioGPU WebSocket Streaming Service - Deployment Guide

## Overview

The BioGPU streaming service provides real-time updates for bioinformatics pipeline processing through WebSocket connections. This guide covers both development and production deployment scenarios.

## Architecture Components

### Infrastructure Services
- **Redis** (port 6379): Pub/sub messaging and job status caching
- **Kafka** (port 9092): Reliable job queue and event streaming
- **Zookeeper** (port 2181): Kafka coordination service
- **WebSocket Server** (port 8765): Real-time client connections
- **API Server** (port 8000): REST endpoints (planned)
- **Metrics Server** (port 8001): Prometheus metrics (planned)

### Service Flow
1. Jobs submitted to Kafka queue
2. Job dispatcher processes jobs from queue
3. Progress updates published to Redis channels
4. WebSocket clients receive real-time updates via subscriptions

## Development Setup

### Prerequisites
```bash
# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # On Linux/Mac
# or
venv\Scripts\activate     # On Windows

# Install dependencies
pip install -r requirements.txt
```

### Starting Services

#### Option 1: Using Docker Compose (Recommended)
```bash
# Start all infrastructure services
docker-compose up -d

# The service will be available with:
# - Redis: localhost:6379
# - Kafka: localhost:9092
# - WebSocket: localhost:8765
# - Kafka UI: localhost:8080
```

#### Option 2: Manual Setup
```bash
# Ensure virtual environment is activated
source venv/bin/activate

# Start infrastructure (Redis, Kafka, Zookeeper)
docker-compose up -d redis kafka zookeeper

# Run the service
python3 run_service.py

# Check if WebSocket server is listening
ss -tulpn | grep 8765
```

### Development WebSocket Connection
```javascript
// Connect from local development
const ws = new WebSocket('ws://localhost:8765');

ws.onopen = () => {
    // Authenticate with JWT
    ws.send(JSON.stringify({
        type: 'auth',
        token: 'your-jwt-token'
    }));
    
    // Subscribe to job updates
    ws.send(JSON.stringify({
        type: 'subscribe',
        job_id: 'job-123'
    }));
};

ws.onmessage = (event) => {
    const update = JSON.parse(event.data);
    console.log('Job update:', update);
};
```

## Production Deployment

### Network Configuration
The WebSocket server binds to `0.0.0.0:8765`, making it accessible from any network interface. This is production-ready but should be secured properly.

### SSL/TLS Setup with Reverse Proxy

#### Nginx Configuration
```nginx
server {
    listen 443 ssl http2;
    server_name api.biogpu.com;

    ssl_certificate /etc/ssl/certs/biogpu.crt;
    ssl_certificate_key /etc/ssl/private/biogpu.key;
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers HIGH:!aNULL:!MD5;

    # WebSocket endpoint
    location /ws {
        proxy_pass http://biogpu-service:8765;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # Timeouts for long-running connections
        proxy_read_timeout 3600s;
        proxy_send_timeout 3600s;
    }
    
    # Optional: REST API endpoint
    location /api {
        proxy_pass http://biogpu-service:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

#### Caddy Configuration (Alternative)
```
api.biogpu.com {
    reverse_proxy /ws biogpu-service:8765 {
        header_up Upgrade {http.request.header.Upgrade}
        header_up Connection {http.request.header.Connection}
    }
    
    reverse_proxy /api/* biogpu-service:8000
}
```

### Production WebSocket Connection
```javascript
// Connect from remote client with SSL
const ws = new WebSocket('wss://api.biogpu.com/ws');

// Or direct connection with port (if SSL configured on service)
const ws = new WebSocket('wss://biogpu.example.com:8765');
```

### Environment Variables for Production

Create a `.env` file or set these in your deployment:
```bash
# Redis Configuration
REDIS_URL=redis://redis-prod:6379

# Kafka Configuration  
KAFKA_SERVERS=kafka-prod-1:9092,kafka-prod-2:9092,kafka-prod-3:9092

# Security
JWT_SECRET=your-production-secret-key-minimum-32-chars

# BioGPU Configuration
BIOGPU_EXECUTABLE=/usr/local/bin/bio_gpu_pipeline
GPU_DEVICE_ID=0
REFERENCE_DB_PATH=/data/microbial_ref_db
RESISTANCE_DB_PATH=/data/quinolone_resistance_db

# Performance Tuning
MAX_RETRIES=3
BASE_RETRY_DELAY=60
LOG_LEVEL=INFO
```

### Docker Production Deployment
```yaml
# docker-compose.prod.yml
version: '3.8'

services:
  biogpu-service:
    image: biogpu/streaming-service:latest
    environment:
      - REDIS_URL=redis://redis:6379
      - KAFKA_SERVERS=kafka:29092
      - JWT_SECRET=${JWT_SECRET}
      - LOG_LEVEL=INFO
    deploy:
      replicas: 2
      restart_policy:
        condition: on-failure
        delay: 5s
        max_attempts: 3
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
```

## Security Considerations

### Production Checklist
- [ ] Use WSS (WebSocket Secure) with valid SSL certificates
- [ ] Configure strong JWT secret (not the default!)
- [ ] Implement rate limiting at reverse proxy level
- [ ] Set up CORS policies if needed
- [ ] Use environment-specific configurations
- [ ] Enable authentication for all endpoints
- [ ] Monitor WebSocket connections and set connection limits
- [ ] Configure firewall rules for service ports
- [ ] Use separate Redis/Kafka instances for production
- [ ] Enable Redis password authentication
- [ ] Configure Kafka SASL/SSL authentication

### Connection Limits
```python
# In production, consider adding to config:
MAX_CONNECTIONS_PER_USER = 5
MAX_TOTAL_CONNECTIONS = 1000
CONNECTION_TIMEOUT = 3600  # 1 hour
```

## Monitoring

### Health Checks
```bash
# Check if service is running
curl http://localhost:8000/health

# Check WebSocket server
wscat -c ws://localhost:8765

# Monitor active connections
redis-cli
> PUBSUB CHANNELS job_updates:*
```

### Logs
```bash
# View service logs
docker-compose logs -f biogpu-service

# Or if running directly
tail -f biogpu_service.log
```

## Troubleshooting

### Port Already in Use
```bash
# Check what's using port 8765
ss -tulpn | grep :8765

# Kill the process if needed
kill -9 <PID>
```

### Connection Issues
1. Verify all services are running:
   ```bash
   docker-compose ps
   ```

2. Check service logs for errors:
   ```bash
   docker-compose logs redis kafka biogpu-service
   ```

3. Test WebSocket connection:
   ```bash
   # Install wscat if needed
   npm install -g wscat
   
   # Test connection
   wscat -c ws://localhost:8765
   ```

### Virtual Environment Issues
Always ensure the virtual environment is activated before running the service:
```bash
# Check if venv is activated
which python
# Should show: /path/to/project/venv/bin/python

# If not activated, activate it:
source venv/bin/activate
```

## Client Integration Examples

### Python Client
```python
import asyncio
import websockets
import json

async def connect_to_biogpu():
    uri = "ws://localhost:8765"  # or wss://api.biogpu.com/ws for production
    
    async with websockets.connect(uri) as websocket:
        # Authenticate
        await websocket.send(json.dumps({
            "type": "auth",
            "token": "your-jwt-token"
        }))
        
        # Subscribe to job
        await websocket.send(json.dumps({
            "type": "subscribe",
            "job_id": "job-123"
        }))
        
        # Listen for updates
        async for message in websocket:
            update = json.loads(message)
            print(f"Received update: {update}")

asyncio.run(connect_to_biogpu())
```

### React/JavaScript Client
```javascript
import { useEffect, useState } from 'react';

function useBioGPUStream(jobId, token) {
    const [updates, setUpdates] = useState([]);
    const [connected, setConnected] = useState(false);
    
    useEffect(() => {
        const ws = new WebSocket(
            process.env.NODE_ENV === 'production' 
                ? 'wss://api.biogpu.com/ws'
                : 'ws://localhost:8765'
        );
        
        ws.onopen = () => {
            setConnected(true);
            ws.send(JSON.stringify({ type: 'auth', token }));
            ws.send(JSON.stringify({ type: 'subscribe', job_id: jobId }));
        };
        
        ws.onmessage = (event) => {
            const update = JSON.parse(event.data);
            setUpdates(prev => [...prev, update]);
        };
        
        ws.onclose = () => setConnected(false);
        
        return () => ws.close();
    }, [jobId, token]);
    
    return { updates, connected };
}
```