# CUDA Perishable Inventory Control

Reproduction of *A CUDA approach to compute perishable inventory control policies using value iteration* (Ortega et al., The Journal of Supercomputing, 2018).

## Scope

- **Two-product** perishable inventory model with **substitution**
- Sequential Value Iteration (VI) implementation (Algorithm 3)
- CUDA GPU implementation (Algorithm 4/5)

## Requirements

- CUDA Toolkit 11+ (paper used CUDA 8)
- GCC/G++
- NVIDIA GPU (Compute Capability 3.7+)

## Project structure

```
cuda_inventory/
├── src/           # Source code
├── scripts/       # Run scripts
├── results/       # Run results and checkpoints
└── docs/          # Paper and references
```

## Build

```bash
make
```

## Run

```bash
# Sequential (P1 -> P2 -> P3 -> P4; skips instances with existing checkpoint)
./build/vi_sequential
```

**OpenMP** (faster on multi-core): on macOS run `brew install libomp`, then:
```bash
make sequential OPENMP=1
```

**CUDA** (requires GPU): `./build/vi_cuda`

## Reference results (paper Table 4)

| Instance | μa | μb | N     | Sequential (s) | GPU (s) | Speedup |
|----------|----|----|-------|----------------|---------|---------|
| P1       | 5  | 5  | 14461 | 357.89         | 91.29   | ~3.9×   |
| P2       | 5  | 6  | 20449 | 698.20         | 136.84  | ~5.1×   |
| P3       | 6  | 6  | 28561 | 1358.04        | 214.71  | ~6.3×   |
| P4       | 7  | 7  | 38416 | 4241.66        | 361.49  | ~11.7×  |

## GitHub

```bash
# After creating an empty repo on GitHub
git remote add origin https://github.com/YOUR_USERNAME/cuda_inventory.git
git push -u origin main
```

## License

Academic reproduction; follows the paper’s Creative Commons Attribution 4.0.
