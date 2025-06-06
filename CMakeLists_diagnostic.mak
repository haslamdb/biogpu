# Makefile for FQ Resistance Diagnostic Tool

NVCC = nvcc
NVCC_FLAGS = -arch=sm_61 -O3 -std=c++14
CUDA_LIBS = -lcudart

# Target executable
TARGET = fq_diagnostic

# Source files
SOURCES = diagnostic_fq_detection.cu

# Build rules
all: $(TARGET)

$(TARGET): $(SOURCES)
	$(NVCC) $(NVCC_FLAGS) -o $(TARGET) $(SOURCES) $(CUDA_LIBS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o test_fq_resistance.fastq

.PHONY: all run clean