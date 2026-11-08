# Jastrow 使用说明

实现日期：2026-11-04。波函数为 exp[-v1 B1-v2 B2-v3 B3] Haf(F)。

```ini
Ifjastrow = 1
v1 = 0
v2 = 0
v3 = 0
samples_sr = 500000
samples = 5000000
```

Ifjastrow默认0；关闭时有效v归零。设1且sr_iterations>0时联合优化
chi1、Delta、delta及v1/v2/v3。sr_iterations=0时使用给定v直接测量。
v可正可负；delta继续严格为正，内部优化log(delta)。优化结束的checkpoint保留开关和v。
冻结delta时使用与其解耦的SR子空间；关闭Jastrow时新增方向置零、与原三维问题解耦。
SR存储统一使用六维；Gaussian核的df仍仅有三个分量。

B1是最近邻，B2是对角近邻，B3是轴向距离2的无向占据键数。
按最小周期距离分类并去重：L=4的正负2方向只计一次；L=3的距离2折回1，不算B3。
单粒子DeltaJ只遍历其他占据粒子，O(Nb)；pair使用完整配置占据键差，O(Nb^2)，
正确处理重合/共享端点。当前pair的Jastrow计算尚未进一步优化成只访问变化站点的局域更新。
所有非对角能量及关联均含exp(-DeltaJ)，Metropolis使用对数接受率中的-2DeltaJ。
Hafnian::amplitude和moved仍是纯Hafnian接口，调用者统一补充Jastrow修正。

samples.csv新增O_v1/O_v2/O_v3；diagnostics.csv新增B1/B2/B3/J；
optimization.csv新增v及六维F/S；summary和参数文件保存Ifjastrow及v。
关联.dat列格式不变。MPI传输55个long double统计值，保留链间均值修正；广播全部六个参数。
运行目录仍采用原Gaussian参数+时间命名；Jastrow开关与参数见logs输入快照及summary。
省略新输入仍可使用旧波函数；分析程序应按列名读取新增列。

小型4进程试跑：

```bash
cd /DATA/cuic/Hardcore_boson_VMC
VMC_RANKS=4 bash run_optimize_mpi.sh examples/jastrow.inp
```

生产运行可在自己的输入中加入上述四个新键。开启Jastrow不会修改已经提交任务的源码快照。
建议先做8x8、Nb=8及不同seed的试验，再开展12x12、Nb=18大任务。
Jastrow不降低Hafnian阶数，也不改变其零点支持；实际自关联和接受率需重新评估。

验证覆盖：Nb=2逐配置精确Hamiltonian能量、非零v的完整单体/pair矩阵Fock枚举、
局域DeltaJ与直接重算、Jastrow导数、MPI单进程/串行一致、六维统计池化、
各链独立重放、SR更新同步、独立SR/最终预算、零v与关闭开关对照。
