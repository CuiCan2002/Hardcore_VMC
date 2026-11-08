# MPI 多链版本

每个 MPI rank 一条独立 Markov 链；所有 rank 在每轮 SR 使用同一组参数。
只在一轮采样结束时汇总，rank 0 求解一次 SR 并广播更新后的参数。
保持有隙 Gaussian 参考态及 `(chi1, Delta, delta)` 参数化，内部仍优化 log(delta)。

## 服务器提交

```bash
cd /DATA/cuic/Hardcore_boson_VMC
bash run_optimize_mpi.sh
# 指定核数和输入：
VMC_RANKS=16 bash run_optimize_mpi.sh examples/optimize.inp
```

默认读取 examples/optimize.inp，在 norm2 单节点启动 32 个 MPI 进程，每进程一个 CPU。
提交器保留原有输入/源码快照、data 参数及时间命名、logs/results 布局和临时编译目录清理。
不要直接对这个提交器执行 sbatch；它会自行提交。
现有 run_optimize.sh 继续运行串行程序。MPI 脚本显式加载
`/opt/intel/oneapi/mpi/latest/env/vars.sh`，以同一套 MPI 编译和运行。

## samples 和 block_size

samples 是最终测量所有链合计的样本数；samples_sr 是每轮 SR 所有链合计的样本数。
省略 samples_sr 或设为0时继承 samples，兼容旧输入。两个阶段各自按下面的规则分配。
按完整 block 分配：令 B=samples/block_size，P=进程数，rank r 获得
`(B/P + (r < B%P))*block_size` 条记录（整数除法）。
总数精确不变；各链样本数最多相差一个 block；block 不跨链，不丢尾样本。
每条链至少需要两个完整 block，否则作业明确失败。统计上应远多于两个 block。
例如 samples=5000000、block_size=100、P=32 时，16条链各156300条，另16条各156200条。

warmup_sweeps、sweeps_between、block_size 和 diagnostic_stride 都是每条链的控制参数。
热化不能除以 P。仍然在每轮 SR 重新初始化和热化链，保持原采样语义。
不要为了满足 block 数约束而把每条链缩得过短；应结合自关联时间选择规模。

## 数值合并与随机数

保留单链在线中心化协方差，按各链实际样本数加权。
全局协方差包含链间均值偏移项，例如
`S = sum_r w_r [S_r + (Omean_r-Omean)(Omean_r-Omean)^T]`；
F 和能量方差采用同样的中心矩合并公式。
根进程汇总统计量并更新参数；不平均不同链的 SR 参数更新量。
通信使用 MPI_LONG_DOUBLE，相关矩阵累计和及权重用 MPI_Reduce 求和。

rank 0 保留原始 seed+iteration；其余 rank 使用该值与 rank 的 SplitMix64 混合获得种子。
各进程独立持有 mt19937_64、configuration、Hafnian 缓存；记录最终链的精确种子。
这些是不同伪随机种子，不是数学保证互不重叠的随机子流。
同输入、同 MPI 规模、同编译环境可复现；不同规模通常产生不同的 SR 随机轨迹。
参数日志和 checkpoint 使用 long double 的 max_digits10 精度，支持参数无损读回。

## 输出

P>1 时：

```text
results/
    summary.txt
    optimization.csv
    checkpoint.inp
    final_parameters.inp
    pooled_statistics.csv
    chain_summary.csv
    one_body_correlation.dat
    nn_pair_correlation.dat
    correlations.csv
    correlations.csv.structure.csv
    chains/
        chain_0000/
            samples.csv
            diagnostics.csv
            diagnostics_info.txt
            chain_info.txt
        chain_0001/ ...
```

diagnostics 和 correlation 文件分别受原有开关控制。
chain_summary.csv 记录 rank、seed、样本数、能量、block误差及接受率。
pooled_statistics.csv 保存最终合并的能量、导数均值、F 和 S，便于独立核验。
summary.txt 的 mpi_ranks 和 samples_scope 标明实际链数及总样本语义。
measurement_seed 为 rank 0 的种子；其他链种子见 chain_summary.csv。
P=1 保留原串行输出布局及采样轨迹。

不同链的时间记录不能首尾连接后计算 ACF。现有 analyze_diagnostics.py 可分别分析每个
chain_XXXX/diagnostics.csv。应比较链间均值和各自自关联，再进行多链诊断。
全局能量误差采用独立链公式 `SE^2=sum_r w_r^2 SE_r^2`，每个 SE_r 来自链内 block；
它仍取决于 block 是否足够长，也不是链已热化的证明。

任一进程出现不可恢复异常时 MPI_Abort 终止整个任务，避免其他 rank 挂在 collective 中。
SR 更新回退失败则广播正常停止状态，所有链共同测量最后有效参数。

## 编译和测试

```bash
source /opt/intel/oneapi/mpi/latest/env/vars.sh
make all mpi test LDFLAGS=-lstdc++fs  # 本服务器 GCC 7 需要 filesystem 链接库
# 在 Slurm 分配的节点中：
python3 tests/test_mpi.py "$PWD" /path/to/new_validation_directory
```

集成测试覆盖：MPI单进程与串行逐字节一致；多链不同样本数；每条链按相同seed串行重放；
用原始样本独立重算 pooled S/F；完整单体和pair矩阵的加权合并；SR逐轮重放；
非法采样预算的全体失败退出。tests/benchmark_mpi.py 用固定总样本量比较1/4/16/32进程。
tests/run_mpi_validation.sh 是服务器短时Slurm验证入口，与生产提交器独立。

## 2026-11-01 验证结果

- 本地 MSVC 串行物理与数值测试通过；服务器 GCC 串行和 MPI 构建无警告，物理测试通过。
- Slurm 74678：全部 MPI 集成测试通过，包括1进程串行一致、2/4进程统计合并、
  4进程不等长链 `[1200,1200,1100,1100]`、SR逐轮重放和异常退出。
- Slurm 74679：通过真实 run_optimize_mpi.sh 提交4进程任务，2轮SR和最终测量正常完成。
- 参数文件/优化日志升级为 max_digits10：避免 Linux long double 参数读回舍入改变轨迹。

节点 c28 短基准：8x8、Nb=8、每轮合计64000个样本、2轮SR、最终完整关联和诊断开启，
warmup=500 sweeps/链/轮。计时含 MPI 启动和输出，不含编译及排队。

|进程数|耗时/秒|相对1进程加速|
|---:|---:|---:|
|1|89.04|1.00|
|4|22.85|3.90|
|16|7.33|12.15|
|32|3.55|25.08|

这是一次短基准，参数轨迹随进程数改变；不保证大规模生产和不同节点具有完全相同的加速比。
服务器测试日志位于 validation_mpi_20261101/check-74678.log，
逐项结果位于 check_74678，基准数据位于 bench_74678/timings.csv。

2026-11-03：新增 samples_sr。生产示例设为500000，samples保持5000000。
两个阶段分别要求完整block且每条MPI链至少两个block；sr_iterations=0时不检查SR预算的block约束。
修改输入不会改变已运行作业；checkpoint记录两个预算，仍默认sr_iterations=0用于固定参数测量。
验证作业74703通过：串行/MPI独立预算、最终样本数、旧输入兼容及全部MPI回归检查。
