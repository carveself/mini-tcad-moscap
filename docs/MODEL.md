# 模型、方程与数值方法

## 1. 坐标、符号和单位

氧化层/硅界面为 $`x=0`$，正方向进入 P 型硅；深处体区接地。$`\psi`$ 为相对于体区的静电势，$`V_G`$ 为栅对体电压，$`V_T=k_BT/q`$ 为热电压。阈值电压写成 $`V_{th}`$，避免和热电压混淆。

内部采用 SI（International System of Units，国际单位制）。受主浓度输入 cm⁻³，乘 $`10^6`$ 得到 m⁻³；厚度输入 nm，乘 $`10^{-9}`$ 得到 m。电荷按单位面积记 C/m²，电容按单位面积记 F/m²。电场为 V/m。

器件参数为 $`N_A`$、$`t_{ox}`$、温度、金属/半导体功函数差 $`\phi_{ms}`$ 和固定氧化层片电荷 $`Q_{ox}`$。输入 `qox_cm2` 表示每 cm² 的正元电荷数，可为负；转为 C/m² 时乘 $`q\times10^4`$。这里将固定电荷集中在界面等效位置。

```math
C_{ox}=\varepsilon_{ox}/t_{ox},\qquad
V_{FB}=\phi_{ms}-Q_{ox}/C_{ox}.
```

Flat-band voltage（平带电压，$`V_{FB}`$）是硅中没有带弯曲时的栅压。正固定电荷使平带电压向负方向移动。

## 2. 电荷与非线性泊松方程

完全电离且无施主时，体区满足电中性与质量作用定律：

```math
p_0=\frac{N_A+\sqrt{N_A^2+4n_i^2}}{2},\qquad
n_0=n_i^2/p_0,\qquad p_0-n_0=N_A.
```

玻尔兹曼统计给出：

```math
u=\psi/V_T,\quad p=p_0 e^{-u},\quad n=n_0e^u,\quad
\rho=q(p-n-N_A),\quad -\varepsilon_{Si}\psi''=\rho.
```

代码用 `expm1` 计算

```math
r(u)=\rho/(qp_0)=\operatorname{expm1}(-u)-(n_0/p_0)\operatorname{expm1}(u)
```

以避免体区附近相近数相减造成精度损失。牛顿迭代试探值超出安全指数范围会被拒绝，代码没有通过裁剪载流子密度来偷偷改变方程。

硅禁带使用 Varshni 型近似；本征浓度在 300 K 标定为 $`10^{10}`$ cm⁻³，温度关系包含 $`T^{3/2}`$ 与禁带指数项。它是教学参数化，不是材料精密标定。

## 3. 边界条件与栅电荷

Robin boundary condition（罗宾边界条件）将表面电势与法向导数联系起来，体现栅电压在氧化层与硅中的分配。远端采用 Dirichlet boundary condition（狄利克雷边界条件），指定体区电势为零。

```math
-\varepsilon_{Si}\psi'(0)=C_{ox}(V_G-V_{FB}-\psi_s),\qquad \psi(L)=0.
```

实际金属栅片电荷与半导体片电荷分别为：

```math
Q_G=C_{ox}(V_G-\phi_{ms}-\psi_s),\qquad Q_S=\int_0^L\rho\,dx.
```

注意有固定氧化层电荷时，表面硅中的位移通量等于 $`Q_G+Q_{ox}`$。有限深度域的电荷守恒检查为：

```math
Q_G+Q_{ox}+Q_S-D(L)=0.
```

这里 $`D(L)`$ 用最后一段网格的电势通量计算；随深度增加应接近零。曲线图中的电场由 `gradient` 得到，边界守恒诊断使用离散通量，两者不要混为一谈。

## 4. 离散与阻尼牛顿法

Finite Volume Method（有限体积法，FVM）在节点控制区间上积分泊松方程。网格靠界面更密，进入体区后逐渐变疏。取 Debye length（德拜长度，电势被载流子屏蔽的典型尺度）归一化位置：

```math
L_D=\sqrt{\varepsilon_{Si}V_T/(qp_0)},\quad y=x/L_D,\quad
\beta=C_{ox}L_D/\varepsilon_{Si}.
```

令左右网格间距为 $`h_l,h_r`$，内部节点的残差为：

```math
F_i=\frac{u_i-u_{i-1}}{h_l}+\frac{u_i-u_{i+1}}{h_r}
-r(u_i)\frac{h_l+h_r}{2}.
```

表面节点包括半个控制区间的电荷：

```math
F_0=\frac{u_0-u_1}{h_0}+\beta\left(u_0-\frac{V_G-V_{FB}}{V_T}\right)
-r(u_0)h_0/2.
```

若遗漏表面半格电荷，界面密度变化很快时会产生显著偏差。远端残差为 $`u_{N-1}`$。解析 Jacobian matrix（雅可比矩阵，即残差对未知变量的导数矩阵）是三对角矩阵，使用带状线性求解器；导数包含电子和空穴的指数项。

```math
J\Delta u=-F,\qquad u_{new}=u+\alpha\Delta u.
```

每步先限制无量纲电势增量，再回溯减小 $`\alpha`$，直到残差下降。归一化残差无穷范数小于 $`2\times10^{-10}`$ 时收敛。偏压扫描使用前一解作初值；独立三状态图的牛顿历史从零初值开始。不存在用解析势数组替代数值迭代的路径。

## 5. 电容、独立参考与阈值

通过微分离散方程得到灵敏度 $`s=du/dV_G`$：

```math
Js=b,\quad b_0=\beta/V_T,\quad b_{i>0}=0,\quad
C_{QS}=C_{ox}(1-V_Ts_0).
```

它与对每个偏压求解后的 $`Q_G`$ 做中心差分相互验证。高频近似保留同一个直流解，在交流雅可比中将电子响应冻结，即把电荷导数替换成仅包含空穴的项。该曲线是明确的冻结电子近似，不代表完整频率响应。

独立参考解通过半无限硅的第一积分求片电荷，再对栅压做一维标量求根：

```math
H(u)=e^{-u}+u-1+(n_0/p_0)(e^u-u-1),\qquad
Q_S=-\operatorname{sgn}(u)\sqrt{2\varepsilon_{Si}qp_0V_TH(u)},
```

```math
V_G=V_{FB}+\psi_s-Q_S/C_{ox}.
```

在平带附近使用级数避免相消。参考解不使用空间网格，不提供牛顿初值；它用于检查空间离散误差。

采用 $`\psi_s=2\phi_F`$ 定义阈值，$`\phi_F=V_T\ln(p_0/n_i)`$。解析估计为：

```math
V_{th}=V_{FB}+2\phi_F+\frac{\sqrt{4\varepsilon_{Si}qN_A\phi_F}}{C_{ox}}.
```

数值阈值通过标量求根寻找有限体积解达到该表面势时的栅压，不是直接返回上述解析公式。

## 6. 参数提取与误差含义

数据：独立第一积分模型，61 个偏压点，固定种子高斯噪声，已知单位面积归一化。真实厚度为 12 nm，真实掺杂为 $`1.2\times10^{16}`$ cm⁻³。Inverse problem（逆问题）使用有限体积模型，在对数参数空间拟合厚度与掺杂，目标是噪声标准差加权残差的平方和。

拟合函数只使用已知的温度、功函数差和氧化层固定电荷；真实厚度与掺杂只在生成数据及最后评分时使用。已知参数由 `known` 对象传入，但其中的厚度和掺杂没有作为残差函数的初值或约束中心。局部协方差来自最优点雅可比矩阵，不包含模型失配、面积误差或未知功函数误差。

数值误差、参数提取误差和物理模型误差必须分开：小残差说明方程求解充分；网格误差说明离散逼近程度；它们都不能证明玻尔兹曼统计或高频假设对真实器件完全有效。`validation.json` 还记录浓度/有效态密度筛查值以及有限域端的电位移。

参考背景：[TU Graz MOS 电容课程](https://lampz.tugraz.at/~hadley/psd/L10/moscap_poisson.php)、[MIT MOS 电容课程](https://ocw.mit.edu/courses/6-012-microelectronic-devices-and-circuits-fall-2009/2a94c4436b68993a04ad47bda22eb067_MIT6_012F09_lec10.pdf)。本项目的离散、验证与合成数据实验以仓库代码和输出记录为准。
