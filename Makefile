# Makefile for CUDA Perishable Inventory VI
# Reproduces Ortega et al. (2018) two-product case

CC = gcc
NVCC = nvcc
CFLAGS = -O3 -Wall
NVCCFLAGS = -O3 -arch=sm_50  # sm_50 for Maxwell+, adjust for your GPU
BUILD_DIR = build
SRC_DIR = src

.PHONY: all clean sequential cuda single

# 'make all' builds sequential; cuda requires nvcc (use cloud GPU)
all: sequential
	@which nvcc >/dev/null 2>&1 && $(MAKE) cuda || echo "nvcc not found - skip cuda (build on cloud GPU)"

sequential: $(BUILD_DIR)/vi_sequential

cuda: $(BUILD_DIR)/vi_cuda

single: $(BUILD_DIR)/vi_single

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Sequential version
$(BUILD_DIR)/vi_sequential: $(SRC_DIR)/vi_sequential.c $(SRC_DIR)/common.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ $^

# CUDA version (nvcc handles .cu and .c)
$(BUILD_DIR)/vi_cuda: $(SRC_DIR)/vi_cuda.cu $(SRC_DIR)/common.c | $(BUILD_DIR)
	$(NVCC) $(NVCCFLAGS) -o $@ $^ -lcudart

# Single-product VI (Paper Section 2.3)
$(BUILD_DIR)/vi_single: $(SRC_DIR)/vi_single.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/vi_sequential $(BUILD_DIR)/vi_cuda $(BUILD_DIR)/vi_single
