# BioGPU Streaming Service

A production-ready streaming service for the BioGPU pipeline that provides real-time updates via WebSocket and Kafka integration.

## Features

- **Real-time Updates**: Stream analysis progress via WebSocket with JWT authentication
- **Job Queue Management**: Kafka-based job processing with retry logic
- **Error Handling**: Intelligent error classification and retry strategies
- **Dead Letter Queue**: Failed jobs are sent to DLQ for manual inspection
- **Scalable Architecture**: Supports multiple GPU workers and horizontal scaling
- **Monitoring**: Built-in health checks and Prometheus metrics

## Architecture

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Client    │────▶│  WebSocket  │────▶│    Redis    │
│   (React)   │     │   Handler   │     │  (PubSub)   │
└─────────────┘     └─────────────┘     └─────────────┘
                           │                     │
                           ▼                     │
┌─────────────┐     ┌─────────────┐            │
│    Kafka    │────▶│  BioGPU     │────────────┘
│  (Job Queue)│     │   Worker    │
└─────────────┘     └─────────────┘
       │                   │
       ▼                   ▼
┌─────────────┐     ┌─────────────┐
│     DLQ     │     │   Results   │
│  (Failed)   │     │  Storage    │
└─────────────┘     └─────────────┘
```

## Installation

### Requirements

- Python 3.12+
- Docker and Docker Compose
- CUDA-capable GPU (for BioGPU workers)
- Redis 7+
- Kafka 3+

### Quick Start

1. Clone the repository and navigate to the service directory:
```bash
cd runtime/common/io
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Start infrastructure services:
```bash
docker-compose up -d redis kafka
```

4. Run the service:
```bash
python run_service.py
```

### Docker Deployment

1. Build the Docker image:
```bash
docker build -t biogpu-service .
```

2. Run with Docker Compose:
```bash
docker-compose up -d
```

## Configuration

The service can be configured via environment variables or a JSON configuration file.

### Environment Variables

- `REDIS_URL`: Redis connection URL (default: `redis://localhost:6379`)
- `KAFKA_SERVERS`: Comma-separated list of Kafka brokers (default: `localhost:9092`)
- `BIOGPU_EXECUTABLE`: Path to BioGPU binary (default: `/usr/local/bin/bio_gpu_pipeline`)
- `JWT_SECRET`: Secret key for JWT authentication
- `GPU_DEVICE_ID`: GPU device ID to use (default: `0`)
- `REFERENCE_DB_PATH`: Path to reference database
- `RESISTANCE_DB_PATH`: Path to resistance database
- `MAX_RETRIES`: Maximum retry attempts for failed jobs (default: `3`)
- `BASE_RETRY_DELAY`: Base delay in seconds for retries (default: `60`)

### Configuration File

Create a `config.json` file based on `config.example.json`:

```bash
cp config.example.json config.json
# Edit config.json with your settings
python run_service.py --config config.json
```

## API Usage

### WebSocket Connection

1. Connect to WebSocket endpoint:
```javascript
const ws = new WebSocket('ws://localhost:8765');
```

2. Authenticate with JWT token:
```javascript
ws.send(JSON.stringify({
    type: 'auth',
    token: 'your-jwt-token'
}));
```

3. Subscribe to job updates:
```javascript
ws.send(JSON.stringify({
    type: 'subscribe',
    job_id: 'job-123'
}));
```

### Job Submission

Submit jobs to Kafka topic `new_jobs`:

```python
from kafka import KafkaProducer
import json

producer = KafkaProducer(
    bootstrap_servers=['localhost:9092'],
    value_serializer=lambda v: json.dumps(v).encode('utf-8')
)

job_data = {
    'job_id': 'job-123',
    'sample_id': 'Patient789_Day5',
    'path_to_r1': '/path/to/sample_R1.fastq.gz',
    'path_to_r2': '/path/to/sample_R2.fastq.gz',
    'output_path': '/path/to/output'
}

producer.send('new_jobs', job_data)
```

## Testing

Run the test suite:

```bash
pytest test_streamer.py -v
```

Run integration tests (requires Docker):

```bash
pytest test_streamer.py -v -m integration
```

## Monitoring

- **Health Check**: `http://localhost:8001/health`
- **Prometheus Metrics**: `http://localhost:8001/metrics`
- **Kafka UI**: `http://localhost:8080` (when using docker-compose)

## Error Handling

The service implements intelligent error classification:

- **Transient Errors**: Network issues, timeouts, GPU memory errors
  - Automatically retried with exponential backoff
  - Up to 3 retry attempts by default

- **Permanent Errors**: Invalid input, missing files, permission issues
  - Sent directly to Dead Letter Queue
  - No automatic retries

- **Unknown Errors**: Unclassified errors
  - Treated as transient by default
  - Logged for investigation

## Production Deployment

For production deployment:

1. Use proper SSL/TLS certificates for WebSocket
2. Configure authentication with your identity provider
3. Set up monitoring and alerting
4. Use persistent volumes for Kafka and Redis
5. Configure resource limits and scaling policies
6. Set up log aggregation

See the production docker-compose file in `requirements.txt` for a complete example.

## Troubleshooting

### Common Issues

1. **WebSocket connection fails**
   - Check JWT token is valid and not expired
   - Verify WebSocket port (8765) is accessible

2. **Jobs stuck in queue**
   - Check Kafka consumer group status
   - Verify GPU workers are healthy
   - Check for errors in worker logs

3. **High memory usage**
   - Adjust Redis memory limits
   - Configure Kafka retention policies
   - Monitor for memory leaks

### Debug Mode

Run with debug logging:

```bash
python run_service.py --log-level DEBUG
```

## License

This service is part of the BioGPU project.
