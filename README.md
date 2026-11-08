# Hard-core boson projected VMC（C++17）

这是依据上一级目录 `hardcore_boson_projected_VMC_notes.pdf` 搭建的可运行首版。实现二维正方晶格周期边界、固定偶数粒子数的 projected Gaussian / Hafnian 波函数，使用原始 hard-core Hamiltonian 测量能量，并可对三个参考态参数进行 stochastic reconfiguration（SR）优化。仅依赖 C++17 标准库。

当前版本已切换到 `(chi1, Delta, delta)`；详细公式、log(delta) 更新及固定 delta 扫描见 [DELTA_PARAMETERIZATION.md](docs/DELTA_PARAMETERIZATION.md)。

## 编译与运行

Linux / WSL / 安装 GNU make 和 GCC 的环境，在本目录执行：

```sh
make
make test
./hardcore_vmc model.inp
```

Windows 的 **x64 Native Tools Command Prompt for Visual Studio** 中：

```bat
nmake /f Makefile.msvc
nmake /f Makefile.msvc test
hardcore_vmc.exe model.inp
```

Windows 已编译的程序也可直接从 PowerShell 启动：

```powershell
.\hardcore_vmc.exe model.inp
.\hardcore_vmc.exe examples/small.inp
.\hardcore_vmc.exe examples/optimize.inp
```

`output` 相对于**启动程序时的工作目录**解析。为避免覆盖数据，输出目录必须尚不存在。交付时已保留 `results`、`results_small` 和 `results_optimized` 的验证结果；重新运行请先将对应 `.inp` 的 `output` 改为新目录，例如 `results_run2`。

GNU Makefile 支持 `CXX`、`CXXFLAGS`、`CPPFLAGS` 和 `LDFLAGS` 覆盖，例如 `make CXX=clang++`。本机验证使用 MSVC；GNU Makefile 尚未在本机实际执行。

## 输入参数

格式为 `key = value`，也接受空白分隔。`#` 后为注释。未知键、重复键和多余内容会报错；未给出的键使用 `include/vmc.hpp` 中默认值。大小写敏感，路径不支持空白或引号。默认模型是笔记的 `JB=0.06, VB=2.5, K=1, nu=1/8`。

| 键 | 默认值 | 含义 |
|---|---:|---|
| `Lx`, `Ly` | 8, 8 | 晶格长度；各自至少 3，避免短周期边界的键简并 |
| `Nb` | 8 | 固定粒子数；本版要求偶数、2–20 且小于格点数 |
| `JB`, `VB`, `K` | 0.06, 2.5, 1 | 原始模型的最近邻跃迁、最近邻密度作用、关联跃迁 |
| `chi1` | 0.05 | 正常最近邻参考关联参数 |
| `Delta` | 0.05 | 异常配对参考参数，不是物理凝聚比例 |
| `delta` | 0.66 | 到离散 Gaussian 稳定边界的正距离；lambda_eff=lambda_c+delta |
| `optimize_delta` | 1 | 1 优化三项；0 固定 delta、仅优化 chi1 和 Delta |
| `max_log_delta_step` | 1 | 内部 log(delta) 单步最大变化，超出则整体缩小更新 |
| `local_probability` | 0.5 | 最近邻移动概率；必须小于 1，保留非局域移动 |
| `warmup_sweeps` | 500 | 每次 SR 采样及最终测量前的热化 sweep 数 |
| `samples` | 2000 | 每次采样的测量数 |
| `sweeps_between` | 2 | 两次测量间的 sweep 数；1 sweep = Nb 次尝试 |
| `block_size` | 100 | 分块误差的块长度，以测量次数计；samples 必须整除且至少两块 |
| `seed` | 12345 | 随机种子；SR 第 t 轮使用 seed+t |
| `sr_iterations` | 0 | 0 表示仅测量；正数表示固定轮数的 SR 优化 |
| `sr_step` | 0.05 | SR 步长 tau |
| `sr_shift` | 0.001 | SR 正则化 epsilon，在内部 log(delta) 坐标使用 S + epsilon I |
| `correlations` | 1 | 最终测量全部单粒子/最近邻配对关联矩阵、平移平均与结构因子；0 可关闭 |
| `output` | results | 新输出目录 |

密度严格由 `nu = Nb/(Lx*Ly)` 得出，不另设可互相矛盾的密度输入。`lambda_eff` 由当前变分参数和稳定边界计算，固定粒子数的物理 Hamiltonian 不包含此化学势项。

默认参考参数只是满足稳定性条件的示例初值，没有宣称为自洽解或最优解。纯配对态在奇数粒子 sector 为零，所以明确拒绝奇数 `Nb`；`Delta=0` 或 `VB=4K` 导致 `Bk` 恒为零，也会报错。

## 实现与笔记的对应

1. `Kernel`：笔记式 (9)–(10)，`g=VB-4K`，`Ak=lambda_eff-2(JB-g*chi1)*gamma-8K*nu*eta`，`Bk=2g*Delta*gamma`。只接受 `Ak>|Bk|` 的稳定有隙分支。直接离散傅里叶求和得到完整长程 `f_ij`，并将对角元素及其导数置零。
2. `Hafnian`：式 (16)、(25)、(30)，以位掩码子集递归和缓存求精确 Hafnian。固定一个配置时复用 minor 计算单粒子移动振幅及三个解析对数导数。
3. `Sampler`：式 (18)，混合最近邻和全局 occupied-to-empty 移动。局域方案先均匀选粒子、再均匀选四个方向；遇到已占据目标直接拒绝，不重新挑选允许邻居。全局方案均匀选择空位。两者都保持对称 proposal，接受率为 `min(1,|Psi_new/Psi_old|^2)`，代码用对数避免平方溢出。
4. `local_energy`：式 (1)、(19)–(23)。最近邻键沿 +x、+y 各计一次；K 项直接遍历中心 i 的四个 `(delta_x, delta_y)`，再取唯一允许的 hard-core 跃迁方向，不额外乘 2。不用平均场 Hamiltonian 替代物理能量。
5. `log_derivative`：式 (26)–(27)，解析求 `df/dalpha` 并与 Hafnian minors 收缩。
6. `sr_direction`：式 (29)，求解 `(S+epsilon I) delta_alpha = -tau F`。注意 F 是能量梯度的一半，代码遵循笔记的步长约定。在线中心协方差减少数值抵消。每次更新通过 log(delta) 更新保持正值；对数步幅过大则整体缩小更新，计算出不可表示的参数时仍会回退。该更新不保证带噪声能量单调下降。
7. 每次 SR 重新初始化、热化并采样；最后一次参数更新后，单独重新采样并输出最终结果。相关函数只在最终采样测量。

## 输出和物理量约定

| 文件 | 内容 |
|---|---|
| `summary.txt` | 最终总能量、每格点能量、分项、分块误差、接受率、gap、初始/最终参数和运行控制 |
| `samples.csv` | 最终测量的逐样本总能量、E_J/E_V/E_K、三个对数导数 |
| `optimization.csv` | 每轮**更新前**的参数、能量、稳定性、步长缩放、F 和 S 上三角；零轮时只有表头 |
| `final_parameters.inp` | 最终模型/变分参数与采样控制，可重新读取；sr_iterations=0 |
| `one_body_correlation.dat` | 所有 i,j 的单粒子关联，三列 i j value |
| `nn_pair_correlation.dat` | 所有最近邻键对的配对关联，五列 i j k l value |
| `correlations.csv` | 平移平均的 G1(r)、G2_xx(r)、`<n_i n_(i+r)>` |
| `correlations.csv.structure.csv` | `S(q)=sum_r exp(i q.r) [<n_0 n_r>-nu^2]`，动量以整数索引输出 |

能量 CSV 中的 E 是**全系统总能量**，终端打印 E/Ns。误差来自相邻 block 均值的样本方差；应增加块长度检查误差是否饱和，并用不同 seed 复核。输出逐样本轨迹可用于后续自相关分析；程序没有自动估计积分自相关时间或判定收敛。

`G1(r)=<b_i^dagger b_(i+r)>`；`G2_xx(r)=<P_i^dagger P_(i+r)>`，`P_i=b_i b_(i+x)`。配对算符先消灭源键粒子再创建目标键粒子，允许源、目标键重合或共享端点。完整 .dat 文件覆盖 x-x、x-y、y-x、y-y 全部最近邻配对通道；CSV 仍保留平移平均 x-x 通道。格点编号从 0 开始，i=x+Lx*y。每条无向最近邻键沿 +x/+y 各列一次，不重复交换端点。文件第一行为 # 注释，其余为数字。采用实波函数的 Hermitian 对称估计量，逐格点结果仍有统计噪声。相关函数和结构因子当前只输出均值，不输出统计误差。没有实现完整配对密度矩阵的最大本征值或尺寸标度分析，不能单凭这些均值宣布纯 pair superfluid。

## 验证及适用范围

`make test` / `nmake /f Makefile.msvc test` 包含：

- 8 粒子全部配对枚举与缓存 Hafnian 的一致性；单粒子移动 minor 比值；秩一核 `(2M-1)!!` 恒等式。
- 三个解析对数导数与中心有限差分一致性。
- 4×4、Nb=2 的 120 个配置上构造显式 Hamiltonian，使用对角键及两个 spectator 的独立 K 项计数，逐配置比较局域能量并检查 Hermiticity。
- 全配置加权得到的精确**变分能量**与 Monte Carlo 结果比较（不是精确基态能量）；精确能量梯度验证 F 的符号和因子。
- SR 线性方程残差、Gaussian 不稳定分支和奇数粒子数拒绝。
- 3×3、Nb=4 全配置求和与全部单粒子及 NN 配对矩阵逐项比较，包含相同、共享端点及不相交键。

本机已用 MSVC C++17 `/O2 /W4` 编译运行上述测试，以及小体系相关函数、默认 8×8/Nb=8 和 10 轮 SR 示例。详见 `docs/validation.md`。

本版采用精确子集 Hafnian，时间/内存随粒子数指数增长；`Nb<=20` 是保护性上限，不代表 20 粒子已具有高效率。推荐从 Nb=2、4、8 起步。当前每次尝试重新建立缓存，没有实现跨配置的快速更新。直接 Fourier 求和为 O(Ns²)，核使用稠密矩阵。上限还要求各边不超过 128、Ns 不超过 4096；接近上限会消耗较多内存。

本版只实现实数、平移不变、笔记中的均匀参考态参数化。矩形有限尺寸允许使用同一 x/y 系数，但有限几何本身未必 C4 对称。未实现凝聚零模/coherent displacement、Jastrow、平均场自洽求解、复杂波函数及自动多链汇总。节点或极强振幅抵消可能带来遍历性和数值问题；使用全局移动不能代替多初值检查。内部使用 `long double`，其精度取决于编译器，MSVC 下与 double 相同。

## 源码组织

```text
include/vmc.hpp          配置和算法接口
src/vmc.cpp              模型、核、Hafnian、采样、SR、测量
src/main.cpp             命令行入口及错误处理
tests/test_vmc.cpp       独立物理/数值检查
model.inp                8×8、Nb=8 默认参数
examples/small.inp       4×4、Nb=2，开启相关函数
examples/optimize.inp    8×8、Nb=8，10 轮 SR
Makefile                GNU make 构建
Makefile.msvc           Windows nmake 构建
docs/                   参考 PDF 文本摘取、核公式页面和验证记录
```

服务器提交和目录布局见 `docs/SLURM_RUN.md`。新版结果目录仅保留 `logs/` 和 `results/`，不保留 `source/`。

每次更新后另存 `checkpoint.inp`；新日志包含物理 delta 和派生 lambda_c/lambda_eff。旧 lambda_eff、stability_margin 输入不再接受。

MPI 多链版本：`make mpi` 编译，服务器使用 `bash run_optimize_mpi.sh`（默认32进程），
或 `VMC_RANKS=16 bash run_optimize_mpi.sh examples/optimize.inp`。
`samples` 表示所有链的总样本数；各链分别热化、记录时间序列，统一合并SR和关联。
用法、统计归并和验证见 [MPI_CHAINS.md](docs/MPI_CHAINS.md)。

最终固定参数链现在默认输出 `diagnostics.csv` 和 `diagnostics_info.txt`，包括 warmup 与 measurement 两阶段的密度差、NA-NB、远距离单粒子和 pair 局域估计。输入 `diagnostics=1`、`diagnostic_stride=1` 控制开关和记录间隔；详见 [TIME_DIAGNOSTICS.md](docs/TIME_DIAGNOSTICS.md)。

独立采样预算：`samples_sr=500000` 控制每轮SR，`samples=5000000` 控制最终测量，均为所有链总数。省略samples_sr或设0则继承samples。


### Pair 测量与 Hafnian 缓存（2026-11-05）

`pair_measure_stride = 10`：每条链每 10 条正式样本测一次完整近邻 pair 矩阵（在第 10、20、… 条样本测量）。设为 1 可恢复逐样本测量。省略时默认 10；必须为正整数，且每条链至少有一个完整测量间隔。`samples` 仍是正式能量、单体矩阵和密度关联的样本预算；`samples_sr` 仍是每轮 SR 预算。这个间隔不影响 SR 更新或 Markov 链轨迹。

单体矩阵按全部正式样本数归一化，pair 矩阵按实际 pair 测量数归一化。MPI 汇总使用各链实际测量数（包含不能整除 stride 的情况）。`results/correlation_counts.txt` 给出两个总数；各链 `chain_info.txt` 记录本链 pair 样本数。原有 `.dat` 列格式不变。远距离关联时间序列仍由 `diagnostic_stride` 控制，独立于完整 pair 矩阵。

双粒子移动通过未移动粒子的子 Hafnian 展开，复用当前构型的缓存，避免对每个目标近邻键重新分配和求完整 Hafnian。Sampler 在拒绝移动后保留全部缓存；接受移动后仅失效包含该粒子的子集，其余子集继续复用。SR 参数改变后构建新的 Kernel 和 Sampler，不跨参数复用缓存。

降低 pair 测量频率减少测量成本，也减少 pair 的原始样本数；是否损失统计精度取决于其自相关时间。不能把更少测量当作误差自动缩小。
