# BioGPU Production Testing Suite
# Comprehensive testing for the streaming service

import pytest
import asyncio
import json
import tempfile
import os
from unittest.mock import Mock, patch, AsyncMock
from datetime import datetime, timedelta
import jwt
import redis
from kafka import KafkaProducer
import docker
from pathlib import Path

# Import our modules
from biogpu_service import (
    ErrorClassifier, ErrorType, JobFailure, RetryManager,
    BioGPUWorkerIntegration, ProductionJobDispatcher,
    BioGPUConfig, AnalysisStage
)

# Test configuration
TEST_CONFIG = BioGPUConfig(
    redis_url='redis://localhost:6380',  # Test Redis instance
    kafka_servers=['localhost:9093'],    # Test Kafka instance
    biogpu_executable='/usr/local/bin/test_bio_gpu_pipeline',
    jwt_secret='test-secret-key-for-testing-only',
    reference_db_path='/tmp/test_ref_db',
    resistance_db_path='/tmp/test_resistance_db'
)

class TestErrorClassifier:
    """Test error classification logic"""
    
    def test_classify_permanent_errors(self):
        """Test permanent error classification"""
        permanent_errors = [
            "File not found: input.fastq",
            "Permission denied accessing file",
            "Invalid FASTQ format detected",
            "Corrupted input file",
            "Malformed sample manifest"
        ]
        
        for error in permanent_errors:
            assert ErrorClassifier.classify_error(error) == ErrorType.PERMANENT
    
    def test_classify_transient_errors(self):
        """Test transient error classification"""
        transient_errors = [
            "Connection refused to database",
            "Network timeout occurred",
            "CUDA out of memory error",
            "GPU memory allocation failed",
            "Temporary resource unavailable"
        ]
        
        for error in transient_errors:
            assert ErrorClassifier.classify_error(error) == ErrorType.TRANSIENT
    
    def test_classify_by_return_code(self):
        """Test classification by process return codes"""
        # Permanent failures
        assert ErrorClassifier.classify_error("Unknown error", 1) == ErrorType.PERMANENT
        assert ErrorClassifier.classify_error("Unknown error", 127) == ErrorType.PERMANENT
        
        # Transient failures
        assert ErrorClassifier.classify_error("Unknown error", 130) == ErrorType.TRANSIENT
        assert ErrorClassifier.classify_error("Unknown error", 143) == ErrorType.TRANSIENT
        
        # Unknown classification
        assert ErrorClassifier.classify_error("Unknown error", 99) == ErrorType.UNKNOWN

class TestRetryManager:
    """Test retry logic and backoff calculations"""
    
    def setup_method(self):
        self.retry_manager = RetryManager(max_retries=3, base_delay=60)
    
    def test_should_retry_transient_error(self):
        """Test retry decision for transient errors"""
        failure = JobFailure(
            job_id="test-job",
            error_type=ErrorType.TRANSIENT,
            error_message="Network timeout",
            retry_count=1,
            max_retries=3
        )
        
        assert self.retry_manager.should_retry(failure) is True
    
    def test_should_not_retry_permanent_error(self):
        """Test no retry for permanent errors"""
        failure = JobFailure(
            job_id="test-job",
            error_type=ErrorType.PERMANENT,
            error_message="File not found",
            retry_count=0,
            max_retries=3
        )
        
        assert self.retry_manager.should_retry(failure) is False
    
    def test_should_not_retry_max_attempts_reached(self):
        """Test no retry when max attempts reached"""
        failure = JobFailure(
            job_id="test-job",
            error_type=ErrorType.TRANSIENT,
            error_message="Network timeout",
            retry_count=3,
            max_retries=3
        )
        
        assert self.retry_manager.should_retry(failure) is False
    
    def test_exponential_backoff_calculation(self):
        """Test exponential backoff with jitter"""
        failure = JobFailure(
            job_id="test-job",
            error_type=ErrorType.TRANSIENT,
            error_message="Network timeout",
            retry_count=2,
            max_retries=3
        )
        
        start_time = datetime.utcnow()
        next_retry = self.retry_manager.calculate_next_retry(failure)
        
        # Should be between 240s (4*60) and 312s (4*60*1.3) from now
        min_delay = timedelta(seconds=240)
        max_delay = timedelta(seconds=312)
        
        actual_delay = next_retry - start_time
        assert min_delay <= actual_delay <= max_delay

@pytest.fixture
async def mock_streaming_service():
    """Mock streaming service for testing"""
    service = Mock()
    service.publish_update = AsyncMock()
    service.initialize = AsyncMock()
    return service

@pytest.fixture
def sample_job_data():
    """Sample job data for testing"""
    return {
        'job_id': 'test-job-123',
        'sample_id': 'Patient789_Day5',
        'path_to_r1': '/tmp/test_r1.fastq.gz',
        'path_to_r2': '/tmp/test_r2.fastq.gz',
        'output_path': '/tmp/test_output',
        'retry_count': 0
    }

class TestBioGPUWorkerIntegration:
    """Test the BioGPU worker integration"""
    
    @pytest.fixture
    def worker_integration(self, mock_streaming_service):
        retry_manager = RetryManager()
        return BioGPUWorkerIntegration(TEST_CONFIG, mock_streaming_service, retry_manager)
    
    def test_build_command(self, worker_integration, sample_job_data):
        """Test command line building"""
        cmd = worker_integration._build_command(sample_job_data)
        
        expected_args = [
            TEST_CONFIG.biogpu_executable,
            '--r1', sample_job_data['path_to_r1'],
            '--r2', sample_job_data['path_to_r2'],
            '--output-dir', sample_job_data['output_path'],
            '--reference-db', TEST_CONFIG.reference_db_path,
            '--resistance-db', TEST_CONFIG.resistance_db_path,
            '--gpu-device', str(TEST_CONFIG.gpu_device_id),
            '--sample-id', sample_job_data['sample_id'],
            '--progress-json',
            '--threads', '8'
        ]
        
        assert cmd == expected_args
    
    def test_is_valid_fastq_valid_file(self, worker_integration):
        """Test FASTQ validation with valid file"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIII\n")
            temp_path = f.name
        
        try:
            result = asyncio.run(worker_integration._is_valid_fastq(temp_path))
            assert result is True
        finally:
            os.unlink(temp_path)
    
    def test_is_valid_fastq_invalid_file(self, worker_integration):
        """Test FASTQ validation with invalid file"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write("This is not a FASTQ file\n")
            temp_path = f.name
        
        try:
            result = asyncio.run(worker_integration._is_valid_fastq(temp_path))
            assert result is False
        finally:
            os.unlink(temp_path)
    
    @patch('asyncio.create_subprocess_exec')
    async def test_run_analysis_success(self, mock_subprocess, worker_integration, sample_job_data, mock_streaming_service):
        """Test successful analysis execution"""
        # Mock successful process
        mock_process = Mock()
        mock_process.wait = AsyncMock(return_value=0)
        mock_process.stdout.readline = AsyncMock(side_effect=[
            b'{"stage": "minimizer_screening", "progress": 0.5, "message": "Processing..."}\n',
            b'{"stage": "completed", "progress": 1.0, "message": "Done", "results": {"total_reads": 1000}}\n',
            b''  # EOF
        ])
        mock_process.stderr.readline = AsyncMock(return_value=b'')
        mock_subprocess.return_value = mock_process
        
        # Create test output files
        os.makedirs(sample_job_data['output_path'], exist_ok=True)
        
        # Mock validation
        with patch.object(worker_integration, '_validate_inputs', AsyncMock()):
            with patch.object(worker_integration, '_load_final_results', AsyncMock(return_value={'test': 'results'})):
                results, failure = await worker_integration.run_biogpu_analysis(sample_job_data)
        
        assert failure is None
        assert results == {'test': 'results'}
        assert mock_streaming_service.publish_update.call_count >= 1

@pytest.mark.integration
class TestIntegrationWithDockerServices:
    """Integration tests that require Docker services"""
    
    @pytest.fixture(scope="class", autouse=True)
    def docker_services(self):
        """Start Docker services for integration testing"""
        client = docker.from_env()
        
        # Start Redis container
        redis_container = client.containers.run(
            "redis:7-alpine",
            ports={"6379/tcp": 6380},
            detach=True,
            remove=True,
            name="biogpu-test-redis"
        )
        
        # Start Kafka containers (simplified single-node setup)
        kafka_container = client.containers.run(
            "confluentinc/cp-kafka:latest",
            environment={
                "KAFKA_ZOOKEEPER_CONNECT": "localhost:2181",
                "KAFKA_ADVERTISED_LISTENERS": "PLAINTEXT://localhost:9093",
                "KAFKA_BROKER_ID": "1",
                "KAFKA_OFFSETS_TOPIC_REPLICATION_FACTOR": "1"
            },
            ports={"9092/tcp": 9093},
            detach=True,
            remove=True,
            name="biogpu-test-kafka"
        )
        
        # Wait for services to be ready
        import time
        time.sleep(10)
        
        yield
        
        # Cleanup
        redis_container.stop()
        kafka_container.stop()
    
    def test_redis_connection(self):
        """Test Redis connectivity"""
        r = redis.Redis(host='localhost', port=6380, decode_responses=True)
        r.set('test_key', 'test_value')
        assert r.get('test_key') == 'test_value'
    
    def test_kafka_producer(self):
        """Test Kafka message production"""
        producer = KafkaProducer(
            bootstrap_servers=['localhost:9093'],
            value_serializer=lambda v: json.dumps(v).encode('utf-8')
        )
        
        test_message = {'test': 'message'}
        future = producer.send('test_topic', test_message)
        
        # This should not raise an exception
        record_metadata = future.get(timeout=10)
        assert record_metadata.topic == 'test_topic'

class TestSecureWebSocketHandler:
    """Test WebSocket security and authentication"""
    
    def test_jwt_token_verification_valid(self):
        """Test valid JWT token verification"""
        from biogpu_streaming import SecureWebSocketHandler
        
        handler = SecureWebSocketHandler(TEST_CONFIG, None)
        
        # Create valid token
        payload = {
            'user_id': 'test-user-123',
            'exp': datetime.utcnow() + timedelta(hours=1)
        }
        token = jwt.encode(payload, TEST_CONFIG.jwt_secret, algorithm='HS256')
        
        # Verify token
        decoded = handler.verify_jwt_token(token)
        assert decoded['user_id'] == 'test-user-123'
    
    def test_jwt_token_verification_expired(self):
        """Test expired JWT token rejection"""
        from biogpu_streaming import SecureWebSocketHandler, AuthenticationError
        
        handler = SecureWebSocketHandler(TEST_CONFIG, None)
        
        # Create expired token
        payload = {
            'user_id': 'test-user-123',
            'exp': datetime.utcnow() - timedelta(hours=1)
        }
        token = jwt.encode(payload, TEST_CONFIG.jwt_secret, algorithm='HS256')
        
        # Should raise authentication error
        with pytest.raises(AuthenticationError):
            handler.verify_jwt_token(token)
    
    def test_jwt_token_verification_invalid_signature(self):
        """Test invalid JWT signature rejection"""
        from biogpu_streaming import SecureWebSocketHandler, AuthenticationError
        
        handler = SecureWebSocketHandler(TEST_CONFIG, None)
        
        # Create token with wrong secret
        payload = {
            'user_id': 'test-user-123',
            'exp': datetime.utcnow() + timedelta(hours=1)
        }
        token = jwt.encode(payload, 'wrong-secret', algorithm='HS256')
        
        # Should raise authentication error
        with pytest.raises(AuthenticationError):
            handler.verify_jwt_token(token)

class TestJobDispatcherLogic:
    """Test job dispatcher retry and DLQ logic"""
    
    @pytest.fixture
    def mock_dispatcher(self, mock_streaming_service):
        dispatcher = ProductionJobDispatcher(TEST_CONFIG)
        dispatcher.streaming_service = mock_streaming_service
        dispatcher.producer = Mock()
        dispatcher.producer.send = Mock()
        return dispatcher
    
    async def test_handle_success(self, mock_dispatcher, sample_job_data):
        """Test successful job completion handling"""
        results = {
            'taxonomy': {'species': ['E. coli', 'K. pneumoniae']},
            'resistance': {'mutations': [{'gene': 'gyrA', 'mutation': 'S83L'}]}
        }
        
        with patch.object(mock_dispatcher, '_store_results', AsyncMock()):
            with patch.object(mock_dispatcher, '_send_completion_notification', AsyncMock()):
                await mock_dispatcher._handle_success(sample_job_data, results)
        
        # Verify completion update was published
        mock_dispatcher.streaming_service.publish_update.assert_called()
        
        # Verify notification was sent
        mock_dispatcher.producer.send.assert_called_with('notifications', {
            'job_id': sample_job_data['job_id'],
            'sample_id': sample_job_data['sample_id'],
            'type': 'completion'
        })
    
    async def test_handle_failure_with_retry(self, mock_dispatcher, sample_job_data):
        """Test failure handling with retry scheduling"""
        failure = JobFailure(
            job_id=sample_job_data['job_id'],
            error_type=ErrorType.TRANSIENT,
            error_message="Network timeout",
            retry_count=1,
            max_retries=3
        )
        
        with patch.object(mock_dispatcher, '_schedule_retry', AsyncMock()):
            await mock_dispatcher._handle_failure(sample_job_data, failure)
        
        # Should publish retry pending update
        mock_dispatcher.streaming_service.publish_update.assert_called()
        
        # Should not send to DLQ
        dlq_calls = [call for call in mock_dispatcher.producer.send.call_args_list 
                    if call[0][0] == 'failed_jobs_dlq']
        assert len(dlq_calls) == 0
    
    async def test_handle_failure_send_to_dlq(self, mock_dispatcher, sample_job_data):
        """Test failure handling when sending to DLQ"""
        failure = JobFailure(
            job_id=sample_job_data['job_id'],
            error_type=ErrorType.PERMANENT,
            error_message="File not found",
            retry_count=0,
            max_retries=3
        )
        
        await mock_dispatcher._handle_failure(sample_job_data, failure)
        
        # Should send to DLQ
        mock_dispatcher.producer.send.assert_called_with('failed_jobs_dlq', {
            'job_data': sample_job_data,
            'failure_info': {
                'error_type': 'permanent',
                'error_message': 'File not found',
                'retry_count': 0,
                'failed_at': pytest.approx(datetime.utcnow().isoformat(), abs=5),
                'original_error': None
            }
        })

# Performance and load testing
class TestPerformanceCharacteristics:
    """Test performance characteristics of key components"""
    
    def test_error_classifier_performance(self):
        """Test error classification performance"""
        import time
        
        error_messages = [
            "Network timeout occurred during processing",
            "File not found: /path/to/input.fastq",
            "CUDA out of memory error",
            "Invalid FASTQ format detected"
        ] * 1000  # 4000 classifications
        
        start_time = time.time()
        
        for msg in error_messages:
            ErrorClassifier.classify_error(msg)
        
        end_time = time.time()
        classification_time = end_time - start_time
        
        # Should classify 4000 errors in less than 1 second
        assert classification_time < 1.0
        
        classifications_per_second = len(error_messages) / classification_time
        print(f"Error classifications per second: {classifications_per_second:.0f}")

# Test fixtures and utilities
@pytest.fixture
def test_fastq_files():
    """Create temporary test FASTQ files"""
    r1_content = """@read1/1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/1
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
    
    r2_content = """@read1/2
CGATCGATCGATCGATCGATCGATCGATCGATCGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/2
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='_R1.fastq', delete=False) as r1:
        r1.write(r1_content)
        r1_path = r1.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='_R2.fastq', delete=False) as r2:
        r2.write(r2_content)
        r2_path = r2.name
    
    yield r1_path, r2_path
    
    # Cleanup
    os.unlink(r1_path)
    os.unlink(r2_path)

# Run tests
if __name__ == "__main__":
    pytest.main([
        __file__,
        "-v",
        "--tb=short",
        "--cov=biogpu_production_final",
        "--cov-report=html",
        "--cov-report=term-missing"
    ])
