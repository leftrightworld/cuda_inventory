#!/bin/bash
# Run vi_sequential (P1, P2, P3, P4; uses checkpoint to skip completed instances)
cd "$(dirname "$0")/.."
./build/vi_sequential
