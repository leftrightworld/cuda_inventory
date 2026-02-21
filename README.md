# CUDA Perishable Inventory Control

复现论文：*A CUDA approach to compute perishable inventory control policies using value iteration* (Ortega et al., The Journal of Supercomputing, 2018)

## 复现范围

- **双产品**易腐库存模型（含替代 substitution）
- Value Iteration (VI) 串行实现 (Algorithm 3)
- CUDA GPU 加速实现 (Algorithm 4/5)

## 环境要求

- CUDA Toolkit 11+ (论文原文 CUDA 8)
- GCC/G++
- NVIDIA GPU (Compute Capability 3.7+)

## 项目结构

```
cuda_inventory/
├── src/           # 源代码
├── scripts/       # 运行脚本
├── results/       # 运行结果
└── docs/          # 论文等资料
```

## 编译

```bash
make
```

## 运行

```bash
# 串行版本
./build/vi_sequential

# CUDA 版本
./build/vi_cuda
```

## 参考结果 (论文 Table 4)

| 实例 | μa | μb | N    | 串行(s) | GPU(s) | 加速比 |
|------|----|----|------|---------|--------|--------|
| P1   | 5  | 5  | 14461 | 357.89 | 91.29  | ~3.9×  |
| P2   | 5  | 6  | 20449 | 698.20 | 136.84 | ~5.1×  |
| P3   | 6  | 6  | 28561 | 1358.04| 214.71 | ~6.3×  |
| P4   | 7  | 7  | 38416 | 4241.66| 361.49 | ~11.7× |

## License

学术复现用途，遵循原论文的 Creative Commons Attribution 4.0。
