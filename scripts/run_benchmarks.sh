#!/bin/bash
# Run P1-P4 benchmarks for sequential and CUDA versions
# Usage: ./scripts/run_benchmarks.sh [sequential|cuda|all]

BUILD_DIR="../build"
RESULTS_DIR="../results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "$RESULTS_DIR"

run_bench() {
    local prog=$1
    local outfile=$2
    echo "Running $prog, output: $outfile"
    "$BUILD_DIR/$prog" 2>&1 | tee "$outfile"
}

case "${1:-all}" in
    sequential)
        run_bench "vi_sequential" "$RESULTS_DIR/sequential_$TIMESTAMP.log"
        ;;
    cuda)
        run_bench "vi_cuda" "$RESULTS_DIR/cuda_$TIMESTAMP.log"
        ;;
    all)
        run_bench "vi_sequential" "$RESULTS_DIR/sequential_$TIMESTAMP.log"
        run_bench "vi_cuda" "$RESULTS_DIR/cuda_$TIMESTAMP.log"
        ;;
    *)
        echo "Usage: $0 [sequential|cuda|all]"
        exit 1
        ;;
esac
