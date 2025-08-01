#!/bin/bash
echo "Installing BioGPU Unified Pipeline (Process Integration)..."

# Check for required executables
echo "Checking dependencies..."
AMR_OK=false
RES_OK=false

if [ -f "/usr/local/bin/amr_detection" ] || find .. -name "amr_detection" -type f -executable 2>/dev/null | head -1; then
    AMR_OK=true
    echo "✓ AMR detection available"
else
    echo "⚠ AMR detection not found - install separately"
fi

if [ -f "/usr/local/bin/clean_resistance_pipeline" ] || find .. -name "*resistance*pipeline*" -type f -executable 2>/dev/null | head -1; then
    RES_OK=true
    echo "✓ Resistance pipeline available"
else
    echo "⚠ Resistance pipeline not found - install separately"
fi

# Install unified pipeline
sudo cp bio_gpu_pipeline /usr/local/bin/
sudo chmod +x /usr/local/bin/bio_gpu_pipeline

echo ""
echo "Installation complete!"
echo "Run: bio_gpu_pipeline --help"
echo ""
if [ "$AMR_OK" = false ] || [ "$RES_OK" = false ]; then
    echo "Note: Some pipeline components are missing. Install them separately."
fi
