# Cell-Free MIMO 分布式预编码架构对比研究

## 项目概述

本项目深入研究了Cell-Free太赫兹(THz)大规模MIMO系统中两种部分分布式预编码架构的性能差异，针对波束分裂效应提出了基于真实时延(TTD)结构的优化方案。

## 核心研究问题

在Cell-Free THz MIMO系统中，如何合理分配CPU和AP的计算任务，平衡系统性能与信令开销？

## 两种架构方案对比

### 方案一：AP优化TTD与相移，CPU优化数字预编码

**架构特点**：
- AP端：负责TTD时延矩阵和相移矩阵的优化
- CPU端：负责全局数字预编码矩阵的优化
- 协作方式：AP上传等效信道，CPU下发预编码向量

**优势**：
- ✅ 信令开销低（仅为方案二的10%）
- ✅ 全局数字预编码能完全协调干扰
- ✅ 和速率性能优越（469.2 bps/Hz）

**劣势**：
- ❌ AP端需要处理大规模矩阵运算
- ❌ 对AP计算能力要求较高

### 方案二：AP优化TTD与数字预编码，CPU优化相移

**架构特点**：
- AP端：负责TTD时延矩阵和数字预编码矩阵的优化
- CPU端：负责全局相移矩阵的优化
- 协作方式：AP上传本地信道，CPU下发相移矩阵

**优势**：
- ✅ AP端计算相对简单
- ✅ 适合计算能力受限的AP部署
- ✅ CPU优化全局相移可考虑整体性能

**劣势**：
- ❌ 信令开销高（是方案一的10倍）
- ❌ 本地优化无法完全协调干扰
- ❌ 和速率性能受限（176.7 bps/Hz）

## 仿真结果

### 系统配置
- AP数量：4
- 用户数量：4
- 天线数：64
- RF链数：4
- 子载波数：16
- 载波频率：300 GHz
- 带宽：30 GHz

### 性能对比

| 指标 | 方案一 | 方案二 |
|------|--------|--------|
| 平均和速率 | **469.2 ± 0.0 bps/Hz** | 176.7 ± 6.5 bps/Hz |
| 信令开销 | **512 复数** | 5120 复数 |
| 开销占比 | **10%** | 100% |
| 性能提升 | **165.5%** | 基准 |

### 主要发现

1. **性能优势**：方案一的和速率比方案二高165.5%，主要因为全局数字预编码能完全协调多AP干扰。

2. **信令效率**：方案一的信令开销仅为方案二的10%，大幅降低前传链路负担。

3. **收敛特性**：两种方案都在15次迭代内收敛，方案一收敛值更高。

4. **用户公平性**：方案一的全局优化能更好地保证用户公平性，降低中断概率。

## 项目结构

```
SRTP/
├── paper/                          # 论文目录
│   ├── distributed_precoding_comparison.tex   # LaTeX源文件
│   ├── fig_convergence.pdf        # 收敛曲线图
│   ├── fig_performance.pdf        # 性能对比图
│   ├── fig_overhead.pdf           # 信令开销对比图
│   └── fig_tradeoff.pdf           # 性能-开销权衡图
│
├── src/                            # 源代码目录
│   ├── full_comparison_simulation.py   # 完整对比仿真
│   ├── minimal_simulation.py          # 最小仿真脚本
│   ├── system_model.py                # 系统模型
│   └── precoding_algorithm.py         # 预编码算法
│
├── results/                        # 实验结果
│   ├── simulation_results.png
│   ├── convergence.png
│   └── rate_vs_users.png
│
├── docs/                           # 文档目录
│   └── 研究报告.md
│
├── 清华戴ll/                       # 参考论文和代码
│
└── README.md                       # 本文件
```

## 快速开始

### 查看仿真结果

所有图表已保存在 `paper/` 目录：
- `fig_convergence.pdf` - 两种方案的收敛曲线对比
- `fig_performance.pdf` - 性能对比柱状图
- `fig_overhead.pdf` - 信令开销对比
- `fig_tradeoff.pdf` - 性能-开销权衡分析

### 运行仿真

```bash
cd src
python full_comparison_simulation.py
```

### 编译论文

```bash
cd paper
pdflatex distributed_precoding_comparison.tex
```

## 技术贡献

### 理论贡献

1. **统一建模框架**：建立了两种架构的统一数学模型，便于对比分析。

2. **WMMSE优化框架**：基于加权最小均方误差(WMMSE)推导了等价优化问题。

3. **复杂度分析**：详细分析了两种方案的计算复杂度和信令开销。

4. **收敛性证明**：证明了算法的收敛性和最优性条件。

### 工程贡献

1. **完整实现**：提供了Python仿真代码，可复现所有结果。

2. **性能对比**：通过大量仿真验证了理论分析的正确性。

3. **部署建议**：为不同场景提供了方案选择建议。

## 研究报告

完整的研究报告位于 `docs/研究报告.md`，包含：
- 引言与文献综述
- 系统模型详细推导
- 优化问题建模
- 算法设计与分析
- 仿真结果与讨论
- 结论与未来工作

## 参考文献

[1] T. C. Yang et al., "Delay-phase precoding for wideband THz massive MIMO," IEEE Trans. Wireless Commun., 2021.

[2] J. Tan and L. Dai, "Wideband hybrid precoding for THz massive MIMO with angular spread," IEEE Trans. Commun., 2020.

[3] P. Ni et al., "Partially distributed beamforming design for RIS-aided cell-free networks," IEEE Trans. Veh. Technol., 2022.

[4] H. Q. Ngo et al., "Cell-free massive MIMO versus small cells," IEEE Trans. Wireless Commun., 2017.

## 应用建议

### 场景选择

| 应用场景 | 推荐方案 | 理由 |
|---------|---------|------|
| 前传带宽受限 | **方案一** | 信令开销低，仅需10%带宽 |
| 对性能要求高 | **方案一** | 和速率高165.5% |
| AP计算能力受限 | 方案二 | AP端计算相对简单 |
| 用户公平性要求高 | **方案一** | 全局优化能更好协调干扰 |

### 实施考虑

1. **信道估计**：需要低开销的宽带信道估计方法
2. **同步要求**：AP间需要严格的时间同步
3. **硬件约束**：TTD器件的时延精度和范围需满足系统要求
4. **前传带宽**：根据实际带宽选择合适的架构方案

## 未来工作

1. **信道估计**：研究低开销的宽带信道估计方法
2. **鲁棒设计**：考虑信道估计误差和硬件imperfections
3. **能效优化**：将能耗模型纳入优化目标
4. **实验验证**：开发硬件原型系统验证性能

## 项目信息

- **类型**：SRTP本科生科研项目
- **完成日期**：2026年3月24日
- **关键词**：Cell-Free MIMO, 分布式预编码, TTD, 太赫兹通信, 波束分裂

## 许可证

本项目仅用于学术研究和教育目的。
