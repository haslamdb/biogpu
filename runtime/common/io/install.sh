#!/bin/bash
# Installation script for BioGPU Streaming Service

set -e

echo "🚀 BioGPU Streaming Service Installation"
echo "========================================"
echo ""

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)

echo "✓ Found Python $PYTHON_VERSION"

if [ "$PYTHON_MAJOR" -ne 3 ] || [ "$PYTHON_MINOR" -lt 12 ]; then
    echo "❌ Error: Python 3.12 or higher is required"
    echo "   Current version: $PYTHON_VERSION"
    exit 1
fi

# Check if we're in a virtual environment
if [ -z "$VIRTUAL_ENV" ]; then
    echo ""
    echo "⚠️  Warning: You're not in a virtual environment"
    echo "   It's recommended to use a virtual environment for installation"
    echo ""
    echo "   To create one:"
    echo "   python3 -m venv venv"
    echo "   source venv/bin/activate  # On Linux/Mac"
    echo "   # or"
    echo "   venv\\Scripts\\activate  # On Windows"
    echo ""
    read -p "Do you want to continue without a virtual environment? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Upgrade pip
echo ""
echo "📦 Upgrading pip..."
python3 -m pip install --upgrade pip

# Install dependencies
echo ""
echo "📦 Installing dependencies..."
echo "   Using requirements_updated.txt for Python 3.12 compatibility"

if [ -f "requirements_updated.txt" ]; then
    python3 -m pip install -r requirements_updated.txt
else
    echo "❌ Error: requirements_updated.txt not found"
    exit 1
fi

# Check for Docker
echo ""
echo "🐳 Checking Docker installation..."
if command -v docker &> /dev/null; then
    DOCKER_VERSION=$(docker --version | awk '{print $3}' | sed 's/,$//')
    echo "✓ Docker $DOCKER_VERSION is installed"
else
    echo "❌ Docker is not installed"
    echo "   Please install Docker from: https://docs.docker.com/get-docker/"
fi

# Check for Docker Compose
if command -v docker-compose &> /dev/null; then
    COMPOSE_VERSION=$(docker-compose --version | awk '{print $3}' | sed 's/,$//')
    echo "✓ Docker Compose $COMPOSE_VERSION is installed"
elif docker compose version &> /dev/null 2>&1; then
    COMPOSE_VERSION=$(docker compose version | awk '{print $4}')
    echo "✓ Docker Compose $COMPOSE_VERSION is installed (plugin)"
else
    echo "❌ Docker Compose is not installed"
    echo "   Please install Docker Compose from: https://docs.docker.com/compose/install/"
fi

# Create necessary directories
echo ""
echo "📁 Creating directories..."
mkdir -p logs data

# Make scripts executable
echo ""
echo "🔧 Setting executable permissions..."
chmod +x run_service.py

echo ""
echo "✅ Installation complete!"
echo ""
echo "📋 Next steps:"
echo "   1. Copy config.example.json to config.json and edit with your settings"
echo "   2. Start infrastructure services: docker-compose up -d redis kafka"
echo "   3. Run the service: python3 run_service.py"
echo ""
echo "   For full stack deployment: docker-compose up -d"
echo ""