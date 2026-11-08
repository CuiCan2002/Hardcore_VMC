# 服务器 Slurm 运行

修改 `/DATA/cuic/Hardcore_boson_VMC/examples/optimize.inp` 后执行：

```bash
cd /DATA/cuic/Hardcore_boson_VMC
bash run_optimize.sh
```

脚本自动提交 sbatch，不需要再执行 sbatch run_optimize.sh。资源为 norm2、1 节点、1 task、1 CPU、4 GB、24 小时上限；可修改脚本开头的指令。

目录名包含尺寸、Nb、JB/VB/K、初始 chi1/Delta/delta、北京时间（含微秒）。最终变分参数见 results/summary.txt。

```text
data/L..._YYYYMMDD_HHMMSS_microseconds/
    logs/
        optimize.original.inp
        run.inp
        run_optimize.sh
        job_id.txt
        status.txt
        environment.txt
        checksums.sha256
        filesystem_probe.log
        slurm-<id>.out
        slurm-<id>.err
        run.log
    results/
        summary.txt
        optimization.csv
        samples.csv
        final_parameters.inp
        one_body_correlation.dat
        nn_pair_correlation.dat
        correlations.csv
        correlations.csv.structure.csv
```

logs/ 保存提交时的原始输入、实际输入、提交脚本及运行记录。实际输入中的 output 改为 results，相对于结果根目录解析。原项目输入不因此改变。文件系统日志中首次未加 -lstdc++fs 的链接失败属于正常兼容性探测。

结果目录不再保留 source/。源码快照临时放到项目共享目录 `.job_builds/run.*`，编译、测试和运行都使用独立快照，结束时自动清理。排队期间修改原项目不影响已提交任务。SIGKILL、节点故障或启动前取消可能留下临时目录，确认对应任务结束后再清理。日志中的源码校验值仅用于记录，源码不长期保留。

`correlations=1` 现在为默认值，示例输入也已开启。仅在 SR 结束后的最终采样测量关联；设为 0 可关闭。新增两个文件的数据列为：

```text
one_body_correlation.dat: i j <b_i^dagger b_j>
nn_pair_correlation.dat:  i j k l <b_i^dagger b_j^dagger b_k b_l>
```

每个文件第一行为 # 注释，后续为数字数据。编号从 0 开始，i=x+Lx*y，x=i%Lx、y=i/Lx（整数除法）。pair 文件中 (i,j) 和 (k,l) 都是无向最近邻键，每条键沿 +x/+y 方向列出一次，包含周期边界，覆盖 x-x、x-y、y-x、y-y 全部通道。包含重合键及共享端点的键；交换键内端点不改变 boson 算符，因此不重复输出反向端点。

Ns 个格点对应 Ns² 行单粒子数据和 (2Ns)² 行配对数据。8×8 时分别为 4096、16384 行（不含注释）。单粒子矩阵的迹等于 Nb；配对矩阵对角元为该键的 <n_i n_j>。当前实波函数采用 Hermitian 对称估计，逐格点矩阵未强制平移平均。文件只输出均值，不包含误差列。原 CSV 继续提供平移平均 G1、x-x G2、密度关联及结构因子。

```bash
squeue -u cuic
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed
```

2026-10-30 验证：作业 74647 完成编译、全部测试、10 轮 SR 及完整关联测量，退出码 0。3×3/Nb=4 完整关联矩阵与独立 Fock-space 精确求和逐项一致；服务器输出还验证了行列数、完整键覆盖、Hermiticity、trace=Nb，以及关联矩阵还原 E_J/E_V 的恒等式。

旧作业 74646 的指定目录已移除 source/，顶层文件移入 logs/，results/ 原文件保留。新旧运行的最终参数、optimization.csv 和 samples.csv 完全一致，因此将新增关联文件补入旧 results/；来源记录在旧 logs/correlation_measurement.txt。旧 summary 属于原始测量记录，没有为了补入关联而覆盖。

新版使用 delta 作为输入，lambda_eff 由稳定边界计算。固定 delta 扫描及参数迁移说明见 DELTA_PARAMETERIZATION.md。也可用 `bash run_optimize.sh 路径/input.inp` 提交指定输入。
