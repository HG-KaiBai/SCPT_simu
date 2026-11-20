将SurvPT更改为版本0.4.1

Simulation 针对 Exp / Weibull / Fiducials 的所有setting. 包含point estimate evaluation 和posterior uncertainty evaluation.

point estimate: 包含survPT方法. (survOPT for right censoring only). Fiducial, NPMLE 和两种 spbayessurv 方法. 对每一个setting, 做 1000 次 simulation replication. 每个模型都产生fitted survival distribution over survival grid, 为 0 到 tau (最大有限观测值) 的 equally spaced 200 grid (包括 0 点就是201个). 比较1000次average SSE, RMSE, RMSE_nolast. 和对应1000次的sd. 结果保存到文件夹point_estimate

posterior estimate: 包含survPT方法, 两种spbayessurv. 对每一个setting, 在同样1000次的simulation replication中, 针对survival time的underlying distribution S 的quantile point. [S^-1(0.75), S^-1(0.5), S^-1(0.25)]. 计算三个点的survival probability posterior mean 和 95% bayesian credible interval. 最后每个replication输出absolute_bias, credible interval length和 indicator of cover the true value. 表格最后一行输出average absolute bias, average length 和 coverage rate. 结果保存到文件夹posterior

posterior 的 bias 改为 absolute_bias = abs(posterior_mean - p)。Average 行和 summary 中的 average_absolute_bias 是逐次 absolute_bias 的平均值，不是 average signed bias 的绝对值。主分析和 sensitivity 均使用此定义。

posterior 同时包含 Fiducial 方法，在相同数据、相同方法种子和同一次拟合中计算三个真实 quantile 点的 fiducial mean、95% fiducial interval、absolute_bias、length、covered 和 Winkler Score。使用默认 mfid=1000、mburn=100；已有 point estimate 仍使用原 fiducial median。CDF 区间转为 survival 区间时上下端点取 [1-upper, 1-lower]。Fiducial 属于主分析比较方法，无 kappa sensitivity 变体。

posterior 的三个 quantile 均增加 Winkler Interval Score：对于 95% credible interval [L,U] 和真值 p（0.75、0.5、0.25），score = (U-L) + 40 * (max(L-p,0) + max(p-U,0))，数值越小越好。每个 replication 输出 score，Average 行及 posterior/summary.csv 输出各 replication score 的平均值；sensitivity 使用相同规则。

参数. survPT, OPT都使用 max_M = 8, max_EM_iter = 2000, EM_tol = 1e-6. survPT 的 posterior sample使用默认的sample 2000 burn 3000. 其他方法的参数都使用package或function的默认值.

survPT有多种方法. prior = 1(Unif), 2(Exp), 3(Weibull_PT) 三种 * flat = True, False 两种，共六个方法。新增方法命名为 Weibull_PT、Weibull_PT_Flat。prior = 3 使用同一 replication 的观测数据调用 SurvPT::estimate(L, R) 得到 Weibull MLE c(lambda, shape)，显式传入 prior_params；lambda 为 inverse scale。用于所有 setting 的 point_estimate、posterior 和 sensitivity。

此外再各自情况对 kappa 的 sensitivity analysis. 考虑 0.05, 0.1, 0.2, 0.3. 这部分的结果可以单独文件夹分开保存. 

六个 PT 方法各自交叉四个 kappa，共 24 个 sensitivity variants；kappa=0.1 精确复用主分析拟合和 posterior。新增方法不得改变旧方法的随机种子编号。所有方法共享同一 replication 的 dataset，data seed 仍为 replication 编号。


2026-11-06：重新安装更新后的 SurvPT 0.4.1，同一版本用 archive hash 区分。
新增 mixture 数据，使用 SurvPT::data_mixture(seed, n, type, q)，三种 censoring
(right/current/interval) × n=200,500 × q=0.3,0.5，共 12 个 setting，追加为 53–64，
每个仍为 1000 replications。mixture 的真分布为包内固定的 0.5 Weibull(5, scale=1)
+ 0.5 Weibull(5, scale=4.5)，不另设 lambda；真 quantile 由该 mixture survival 求根。
所有主方法以及 kappa sensitivity 均覆盖 mixture；point_estimate 和 posterior 格式不变。

2026-11-06 错误处理：MPT/TBP 任一方法在单次 replication 的拟合、预测或结果
校验中抛出 error，只跳过该方法该次结果，同数据的其他方法继续。失败结果保留
原 replication 编号，估计值、误差、区间、coverage indicator、Winkler score 填 NA。
均值、SD、coverage rate 排除失败结果，summary 的 replicates 为该方法实际成功数；
全部失败时有效数为 0、汇总指标为 NA。warnings 不单独触发跳过，不改模型参数。
汇总同时写 fit_failures.csv（逐条原始错误）、fit_failure_summary.csv（按方法及
错误原因计数）与 fit_status_summary.csv（各 setting/method 的尝试、成功、失败数量）。
历史冻结运行使用单独保存并校验 hash 的恢复入口，保留原 source/manifest。

gamma_sensitivity：仅对非 flat 的 Unif、Exp、Weibull_PT，考虑 gamma=0.5,1,2，
固定 kappa=0.1，共九种方法，覆盖全部 64 个 setting。flat 覆盖 gamma，不纳入。
kappa sensitivity 仍固定 gamma=0.5，不做 gamma × kappa 交叉。gamma=0.5 精确复用
主分析结果及 posterior；新增方法编号追加，旧 dataset 和 method seed 编号不变。
gamma_sensitivity/point_estimate 和 gamma_sensitivity/posterior 保存相同指标。
已有运行按独立 supplement 补充数据/结果，保留原始冻结记录；下一次 Submit.sh 完整运行一次输出全部。


2026-11-06 代码整理：新实验统一使用项目根目录 Submit_Simulation.sh，默认一次完成
64 个 setting × 1000 replications，以及全部主方法、kappa/gamma sensitivity、
point/posterior 和失败报告。新实验不需要后补 Weibull、Fiducial、mixture/gamma
或手工 backfill。运行时只打包统一流程及必要验证代码。历史补算/恢复工具保留在
已有 run 的冻结 source 中；正在排队的恢复作业继续使用其冻结副本。
保持数据生成、setting/method 编号、随机种子、模型参数、采样参数及统计定义，
因此在相同 R/package/build/RNG 环境下，科学结果与既有恢复版本一致。
运行耗时、创建时间及 manifest ID 不要求逐字一致。


2026-11-07：全量从头运行 SurvPT 0.4.2。mixture 真分布同步改为等权
Weibull(shape=3,scale=3) + Weibull(shape=4,scale=10)，真 quantile 同步求根。
使用新版 data_ic 修复后的生成流程，不复用 0.4.1 数据或拟合结果。
所有 survPT 变体及 OPT 的 error 均记录并隔离；PT posterior 超过 tau 的 NA
保留有效 point curve 和其余 quantile。失败覆盖指标为 NA，不当作未覆盖。
汇总每个 setting 的 package tau 和 Exact/Left/Interval/Right 比例向量均值。
true_S_at_average_tau 是将 average_tau 代入真实生存函数，绝非 mean(S(tau))。
posterior 同时输出 estimate 与 mean，分别计算 estimate_absolute_bias 和 absolute_bias；
后者继续表示 posterior mean 的 absolute bias。Average 行和 summary 均分别汇总。
PT estimate 直接取 predict 的 MAP-EM estimate；Fiducial 取原 median；MPT/TBP 取 Shat。
旧 run_20261104_0.4.1 原封保留，新运行使用独立 library、source、manifest 和 archive。


2026-11-08：posterior estimate/mean/lower/upper 的概率边界允许 1e-8 浮点误差；
只把容差内的有限越界值截回 [0,1]，超过容差仍报错，NA 处理保持不变。
合法值逐位保留，所有 posterior 指标使用纠正后的概率。point grid 保留既有
1e-8 校验及原始数值，不额外截断。新 manifest v6 冻结此策略，模型、数据与方法
seed、采样次数不变。完整运行使用统一 source、manifest、chunks 和 COMPLETE 结构；续跑和汇总均使用
该运行 source/Pipeline.R，无需额外补全入口。
