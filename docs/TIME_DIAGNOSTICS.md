# 热化和自关联的时间记录

新增输入（默认启用）：

```ini
diagnostics = 1
diagnostic_stride = 1
```

正式结果目录增加 `results/diagnostics.csv` 和 `results/diagnostics_info.txt`。记录的是 **SR 后最终固定参数的采样链**，包括该链的初始配置、热化期间和正式测量；不把不同参数的 SR 轮拼接在一起。设 diagnostics=0 可关闭，且与 correlations=0/1 独立。

`diagnostic_stride=1` 每个正式测量样本记录一次；设为 m 则每 m 个样本记录一次。时间间隔是 `sweeps_between*m` 个 sweep。1 sweep 为 Nb 次移动尝试，不是一次接受的移动。诊断只读取配置，不消耗随机数，不修改采样链。热化记录间隔同样是 sweeps_between*m，另保存 sweep=0 的初始配置；非整除时不强行添加不等间隔的末行。

CSV 每行是当前配置的瞬时估计量，不是运行平均：

| 列 | 定义 |
|---|---|
| phase | warmup 或 measurement；正式统计和自关联脚本默认排除 warmup |
| sample | 正式测量序号；warmup 行为 0 |
| sweep | 从该链初始化起累计的 sweep 数 |
| E_total | 当前配置的总局域能量 |
| n0, nx, nfar | 格点 0、0+x、r 的占据数，均为 0 或 1 |
| density_diff_x | n0-nx |
| density_diff_far | n0-nfar |
| NA, NB | 棋盘染色 A/B 上的粒子数，A: (x+y)%2=0 |
| NA_minus_NB | NA-NB |
| NA_minus_NB_sq | (NA-NB)^2 |
| G1_far_0 | 固定起点 0，终点 r 的单粒子 Hermitian 局域估计 |
| G1_far_x | 两个端点均平移 +x 后的同一估计 |
| G1_far_diff | 上面两个瞬时估计之差，直接用于对称性残差及其误差分析 |
| G1_far_avg | 相同位移 r 的全起点空间平均估计 |
| G1_other_avg | 位移 r' 的全起点空间平均估计 |
| G2xx_far_0 | x 方向近邻 pair 从 0 到 r 的 Hermitian 局域估计 |
| G2xx_far_x | 源键与目标键均平移 +x 后的估计 |
| G2xx_far_diff | 上面两个 pair 瞬时估计之差 |
| G2xx_far_avg | 相同位移 r 的 x-x pair 全起点空间平均 |

这里 r=(floor(Lx/2),floor(Ly/2))，r'=(floor(Lx/2),floor(Ly/2)+1)，位移按周期边界取模。8×8 时分别为 (4,4)、(4,5)；后者在偶数晶格上与前者具有不同的子晶格奇偶性。格点编号 i=x+Lx*y，从 0 开始。pair 算符为 P_i=b_i b_(i+x)，r 表示两条键起点之间的位移。信息文件还保存具体编号、时间间隔和实际随机种子。

单粒子与 pair 均使用当前实波函数的 (O+O†)/2 估计，和完整关联矩阵采用相同约定。局域估计可超出 0–1 范围，这是振幅比估计的性质，不应裁剪。stride=1 时，诊断列的正式测量均值与完整矩阵相应元素/空间平均一致。

只有 Lx、Ly 都为偶数时，周期晶格的 A/B 染色才构成真正二分图；奇数尺寸时 NA、NB 仅是染色计数，不能当作 JB=0 的严格守恒 sector。

## 分析方式

```bash
python3 tools/analyze_diagnostics.py /路径/results/diagnostics.csv
# 额外丢弃最前面的 10000 个正式记录，再做独立报告：
python3 tools/analyze_diagnostics.py /路径/results/diagnostics.csv \
  --discard 10000 --output-dir /路径/results/diagnostics_discard10000
```

输出目录必须尚不存在，避免覆盖以前报告。脚本生成：

- summary.csv：每个观测量的均值、自关联因子 tau、对应 sweep 时间尺度、ESS、MCSE 和窗口敏感性。
- blocking.csv：不同块长度的误差估计、有效使用样本数及由 blocking 推算的相关因子，用于检查误差平台。
- segments.csv：按时间分成十段的均值；这是漂移检查，不是十条独立链。
- 每个非恒定量的 *_acf.dat：延迟记录数、延迟 sweeps、自关联函数。
- metadata.json：排除的 warmup 行数、额外丢弃数、时间单位和算法约定。

FFT 自关联需要 NumPy。当前服务器系统 Python 没有 NumPy，仍可运行脚本获得 blocking 和时间分段报告，脚本明确标记 FFT/ESS 未计算；可将数据下载到有 NumPy 的本机运行完整分析。本机 Codex bundled Python 已具备 NumPy。程序本身仍只依赖 C++17 标准库。

采用 tau=1+2 sum rho，ESS=N/tau，MCSE=sqrt(Var(X)*tau/N)。使用自洽截断窗口 W>=c*tau(W)，分别计算 c=5、10；窗口不稳定或链太短会标记。某列恒定时不报告 tau=0/ESS无穷大：它可能是真正守恒，也可能只是尚未访问其他 sector，需要结合模型和其他初态判断。

热化是否充分不能由单链某个 tau 自动认证。应固定最终参数（sr_iterations=0），比较不同 seed/初始配置、热化长度和误差分块平台。能量的 tau 不代表密度差或长程关联的 tau；优先检查 *_diff 和 NA_minus_NB 这些可能更慢的量。不要强行令空间对称性平均后的差为零，再用它诊断热化。

`examples/diagnostics_JB0.inp` 使用已完成 JB=0 任务的最终参数，设置 sr_iterations=0、20000 样本，用于验证记录和分析流程。它不是新的优化或充分热化结论。已有旧任务没有存下这些时间序列，无法从其平均矩阵补造历史记录；需要重新运行测量。


## 验证记录（2026-10-31）

服务器作业 74666 在指定 JB=0 最终参数下生成 251 条 warmup 记录和 20000 条正式记录，未重新执行 SR；作业 74667 使用相同参数和种子关闭诊断。两次运行的 samples.csv 完全相同，证明记录功能不改变能量/导数轨迹。原历史任务的参数经历了文本精度舍入，重新读入后不能保证历史轨迹永久逐位一致；新增记录不用于伪造旧任务的历史。

逐样本检查通过：NA+NB=8、粒子数差及平方、G1/G2 平移差、正式记录时刻和热化行数。时间均值与完整矩阵对应元素及全起点空间平均在 1e-11 内一致。额外测试通过 diagnostic_stride=3、非整除热化长度，以及 correlations=0 时独立开启 diagnostics 的情形。

分析器通过已知 AR(1) 系统（phi=0.8，理论 tau=9）、warmup 排除、追加 discard、时间单位转换以及恒定序列标记测试。服务器无 NumPy 时的 blocking 降级路径也已运行。

74666 的完整 FFT 分析在本机 NumPy 环境完成，并复制到该任务 results/diagnostics_analysis/；服务器标准 Python 另生成 results/diagnostics_blocking/。短跑结果仅展示诊断流程，不能据此认证热化；当前只有一条链。能量 tau 约 19 个记录间隔，density_diff_x 约 23、density_diff_far 约 30；量纲转换为 sweeps 时乘以本次间隔 2。
