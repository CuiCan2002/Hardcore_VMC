# 正 gap 的 delta 参数化

物理输入坐标现在为 `(chi1, Delta, delta)`。大写 `Delta` 是 anomalous pairing 参数，小写 `delta` 是到辅助 Gaussian 稳定边界的正距离，两者不同。

```ini
chi1 = 0.05
Delta = 0.05
delta = 0.66
optimize_delta = 1
max_log_delta_step = 1
```

对实际有限晶格的全部离散动量定义

\[
g=V_B-4K,\quad J_{\rm eff}=J_B-g\chi_1,\quad
B_k=2g\Delta\gamma_k,
\]
\[
c_k=2J_{\rm eff}\gamma_k+8K\nu\eta_k+|B_k|,
\qquad \lambda_c=\max_k c_k,
\qquad \lambda_{\rm eff}=\lambda_c+\delta.
\]

于是所有模式的稳定距离为

\[
m_k=A_k-|B_k|=(\lambda_c-c_k)+\delta>0,
\qquad \min_k m_k=\delta.
\]

注意 delta 不是 Bogoliubov 能隙本身；能隙为

\[
E_k=\sqrt{m_k}\sqrt{m_k+2|B_k|}.
\]

代码直接用这个正数表达式计算，避免先做 `A-|B|` 或 `A*A-B*B` 引起小量抵消。旧的 `stability_margin=1e-8` 硬阈值已移除；现在只要求 delta 严格正、参数和计算结果有限、gap 能以机器精度表示。当 delta 极小时，打印的 lambda_eff 与 lambda_c 可能舍入成相同数字，应读取单独保存的 delta/min_A_minus_absB，而不要相减恢复它。

## SR 和导数

因为 lambda_c 也随 chi1、Delta 改变，导数包含完整链式项：

\[
\partial_{\chi_1} A_k=\partial_{\chi_1}\lambda_c+2g\gamma_k,
\quad \partial_\Delta A_k=\partial_\Delta\lambda_c,
\quad \partial_\delta A_k=1.
\]

唯一最大模式 Q 上，`d lambda_c/d chi1 = -2g gamma_Q`，`d lambda_c/d Delta = sign(B_Q)*2g gamma_Q`。对数导数和解析能量梯度已经在两种不同最大动量分支上通过有限差分验证。

当多个模式同时达到最大值时，max 函数可能不可微。代码以机器精度识别相等的最大值，平均这些分支的梯度作为 subgradient；对称等价分支具有相同梯度，真正不同分支交叉时属于分段光滑优化。`boundary_modes` 输出被识别为最大值的模式数。

为了保证 delta 始终为正，SR 实际使用内部坐标 `(chi1, Delta, rho)`，其中 `rho=log(delta)`：

\[
O_\rho=\delta O_\delta,\qquad
S^{(\rho)}=J S^{(\delta)}J,\quad F^{(\rho)}=J F^{(\delta)},
\quad J=\mathrm{diag}(1,1,\delta).
\]

求解 `(S_rho + sr_shift*I) d_beta = -sr_step*F_rho` 后更新

\[
\delta_{\rm new}=\delta_{\rm old}\exp(s\,d\rho).
\]

`max_log_delta_step=1` 限制一次 log(delta) 的变化不超过 1，即 delta 的单步倍数在 exp(-1) 到 exp(1) 之间。如果超过限制，以同一比例 s 缩小整个三参数更新向量。该限制不是 delta 的下界。正则化作用在内部 log 坐标上，因此新旧坐标下相同的 sr_shift/sr_step 不代表相同更新轨迹。

`optimize_delta=0` 表示固定 delta，只求解 chi1、Delta 的 2×2 SR 子系统；不是先解三参数问题再丢掉第三个分量。它适用于 delta 扫描。其他参数变化时 lambda_c、lambda_eff 仍会重新计算，因此 delta 精确保持指定值。

输入和 checkpoint 使用物理 delta。`samples.csv` 中的 O_delta、`optimization.csv` 中的 F_delta/S02/S12/S22 也对应物理 delta 导数；求解器内部再做上面的坐标变换。优化日志额外包含派生 lambda_c、lambda_eff、实际 delta_log_change 和 boundary_modes。

每次参数更新后写 `checkpoint.inp`。SR 更新遇到无法表示的参数等异常时，保留最后有效参数并进入最终测量，summary 的 `sr_stop_reason` 明确标记提前停止，不能据退出码 0 认定 SR 已收敛。物理测量本身失败仍会报错。终端每轮立即刷新，并打印最终测量开始提示。

## 旧参数迁移

旧输入中的 `lambda_eff` 不再作为可优化输入接受，`stability_margin` 也不再接受，错误提示会要求迁移，避免将两个含义混用。迁移需先按同一模型和离散动量网格计算 lambda_c，再用 `delta=lambda_eff-lambda_c`。

默认 8×8、Nb=8、JB=0.06、VB=2.5、K=1、chi1=Delta=0.05 下，lambda_c=1.84，因此旧 lambda_eff=2.5 对应 delta=0.66。服务器当前 optimize.inp 已按此转换，保留 samples=1000000、sr_iterations=100；旧文件保存在项目 `.backups/delta_*` 中，已有结果不改写。

`final_parameters.inp` 与 checkpoint 均可用于后续测量，默认写 sr_iterations=0，保留采样控制。需要继续优化时自行设置轮数，并使用新的输出目录。

## 固定 delta 扫描

```bash
cd /DATA/cuic/Hardcore_boson_VMC
# 准备输入，不提交：
python3 tools/scan_delta.py --base examples/delta_scan.inp
# 准备并顺序提交，避免同时占据多个整节点：
python3 tools/scan_delta.py --base examples/delta_scan.inp --submit
# 自定义扫描：
python3 tools/scan_delta.py --base examples/delta_scan.inp \
  --deltas 1e-3 1e-5 1e-7 1e-9 --samples 100000 --sr-iterations 100 --submit
# 收集结果：
python3 tools/scan_delta.py --collect scan_inputs/时间戳/manifest.csv
```

默认扫描 1e-2、1e-4、1e-6、1e-8、1e-10，每个点从相同 chi1/Delta 开始、强制 optimize_delta=0，并测量全部关联。每点的结果分别保存在 data 中带 delta 的目录，仍仅包含 logs/ 和 results/。run_optimize.sh 现在接受可选输入路径；不传参数仍读取 examples/optimize.inp。

manifest.csv 记录输入、作业 ID 和结果目录，scan_summary.csv 汇总能量及其分块误差、最终参数、gap、G1(Lx/2,Ly/2)、G2_xx(Lx/2,Ly/2)。完整三列和五列关联矩阵仍在各自 results/ 下。相关均值当前没有误差列，单个远距离值也不能替代完整关联矩阵比较或尺寸标度。

扫描示例使用旧作业末段参数热启动、每点 20 轮/20000 样本，属于初步检查。应检查每点优化是否充分、不同 seed、块长度以及必要的其他初值；不能仅凭小 delta 的几次短跑就宣布极限收敛或凝聚。
