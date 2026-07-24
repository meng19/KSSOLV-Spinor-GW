# ISDF 插值基构造与 `c1 / c2` 推导

本文说明 `src/GW/ISDF/isdf_product_gram.m` 和
`src/GW/ISDF/isdf_zeta_g_from_product_gram.m` 中

```matlab
c1 = (phi * phi_mu') .* (psi * psi_mu');
c2 = (phi_mu * phi_mu') .* (psi_mu * psi_mu');
zeta_real = c1 / c2;
```

的数学含义，以及它和论文第三章 ISDF 公式、柯西积分降标度公式之间的关系。

## 1. 乘积态矩阵

GW 极化函数和自能中反复出现成对波函数乘积。以标量形式为例，定义

$$
P(r, ij) = \phi_i(r)\psi_j(r).
$$

其中：

- $r$ 是实空间网格点；
- $i$ 和 $j$ 是两组能带指标；
- $\phi$ 可以是已经取过共轭的波函数，例如 epsilon 中 `phi = conj(valence_real)`；
- $\psi$ 是另一组波函数，例如 conduction states。

如果共有 $N_\mathrm{grid}$ 个实空间点，$N_\phi$ 个 $\phi$ 轨道，$N_\psi$ 个 $\psi$ 轨道，则

$$
P \in \mathbb{C}^{N_\mathrm{grid} \times N_\mathrm{pair}},
\qquad
N_\mathrm{pair} = N_\phi N_\psi.
$$

$P$ 的每一列是一个乘积态 $\phi_i(r)\psi_j(r)$。

这里的乘积是同一个实空间网格点上的逐点乘法。若 $r_g$ 表示第 $g$ 个实空间网格点，则

$$
P(g,ij)=\phi_i(r_g)\psi_j(r_g).
$$

在 MATLAB 中，如果 `phi(:, i)` 和 `psi(:, j)` 都是 $N_\mathrm{grid}\times 1$ 列向量，则单个乘积态应写成

```matlab
P(:, ipair) = phi(:, i) .* psi(:, j);
```

也就是说，单个 $P(:,ij)$ 的构造使用的是 `.*` 点乘，而不是矩阵乘法。

## 2. ISDF 近似

ISDF 选取 $N_\mu$ 个插值点 $r_\mu$，并用这些插值点上的乘积态值表示所有乘积态：

$$
P(r, ij) \approx \sum_{\mu=1}^{N_\mu} Z(r,\mu) P(r_\mu, ij).
$$

矩阵形式为

$$
P \approx Z P_\mu,
$$

其中

$$
Z \in \mathbb{C}^{N_\mathrm{grid} \times N_\mu},
\qquad
P_\mu \in \mathbb{C}^{N_\mu \times N_\mathrm{pair}}.
$$

$Z(r,\mu)$ 就是 ISDF 辅助插值基 $\zeta_\mu(r)$。

注意，辅助基 $Z$ 不带 $ij$ 指标。ISDF 的目标不是为每一个带对 $ij$ 单独构造一套 $\zeta_\mu^{ij}(r)$，而是为所有乘积态张成的整体空间构造一套公共辅助基：

$$
\{\phi_i(r)\psi_j(r)\}_{ij}
\quad\Longrightarrow\quad
\{\zeta_\mu(r)\}_{\mu=1}^{N_\mu}.
$$

不同带对 $ij$ 的差别放在插值点上的系数

$$
P(r_\mu,ij)=\phi_i(r_\mu)\psi_j(r_\mu)
$$

里，而不是放在不同的辅助基里。因此同一组 $\zeta_\mu(r)$ 可以服务于所有 $ij$。

如果为每个 $ij$ 都单独求一套辅助基，

$$
P(r,ij)\approx \sum_\mu \zeta_\mu^{ij}(r)c_\mu^{ij},
$$

那么辅助基数量会随带对数 $N_\phi N_\psi$ 增长，基本失去 ISDF 降标度的意义。ISDF 的压缩来自这样一个事实：大量乘积态虽然列数很多，但它们通常近似落在一个低维公共子空间中。

这对应论文第三章中的 ISDF 乘积态分解：

$$
\rho_{ij}(r)
= \psi_i^*(r)\psi_j(r)
\approx
\sum_{\mu=1}^{N_\mu}
\zeta_\mu(r)\psi_i^*(r_\mu)\psi_j(r_\mu).
$$

## 3. 最小二乘解

$Z$ 通过最小二乘拟合得到：

$$
\min_Z \left\| P - ZP_\mu \right\|_F.
$$

这个 Frobenius 范数包含对所有带对列的误差：

$$
\left\| P - ZP_\mu \right\|_F^2
=
\sum_{ij}
\left\|
P(:,ij)-ZP_\mu(:,ij)
\right\|_2^2.
$$

所以这里的 $ij$ 求和不是说单个乘积态 $P(r,ij)$ 本身含有求和；单个乘积态没有求和，就是逐点乘积 $\phi_i(r)\psi_j(r)$。求和出现在求公共辅助基 $Z$ 的最小二乘问题中，因为 $Z$ 要同时拟合所有带对。

如果只取一个带对，例如 $i=1,j=1$，则得到的问题是

$$
\min_Z
\left\|
P(:,11)-ZP_\mu(:,11)
\right\|_2^2.
$$

这只能保证一个乘积函数被拟合好，不能保证同一组 $Z$ 能表示 $P(:,12)$、$P(:,21)$ 或其他带对。因此，GW 中需要的公共辅助基必须由所有相关的 occupied-unoccupied 或 band-state 乘积共同决定。实际实现中也可以用随机抽样的一部分带对来近似这个整体拟合，但其目标仍然是近似整体乘积空间，而不是只拟合某一个带对。

对 $Z$ 求正规方程：

$$
(P - ZP_\mu)P_\mu^* = 0.
$$

于是

$$
PP_\mu^* = ZP_\mu P_\mu^*.
$$

如果 $P_\mu P_\mu^*$ 可逆，则

$$
Z = PP_\mu^*(P_\mu P_\mu^*)^{-1}.
$$

这就是代码里 `zeta_real = c1 / c2` 的来源。

MATLAB 中：

```matlab
Z = c1 / c2
```

表示求右线性方程

$$
Zc_2 = c_1,
$$

也就是

$$
Z = c_1 c_2^{-1}.
$$

## 4. 为什么不是 `Z = P * inv(P_mu^*)`

一般不能写成

$$
Z = P(P_\mu^*)^{-1}.
$$

原因是 $P_\mu^*$ 通常不是方阵。维度为

$$
P_\mu \in \mathbb{C}^{N_\mu \times N_\mathrm{pair}},
\qquad
P_\mu^* \in \mathbb{C}^{N_\mathrm{pair} \times N_\mu}.
$$

做降标度时通常有

$$
N_\mu \ll N_\mathrm{pair}.
$$

所以 $P_\mu^*$ 是一个高矩阵，没有普通逆。

ISDF 实际使用的是 $P_\mu$ 的右伪逆：

$$
P_\mu^+ = P_\mu^*(P_\mu P_\mu^*)^{-1}.
$$

于是

$$
Z = PP_\mu^+
  = PP_\mu^*(P_\mu P_\mu^*)^{-1}.
$$

注意：

$$
P_\mu^*P_\mu \ne I,
\qquad
P_\mu P_\mu^* \ne I.
$$

一般都不是单位矩阵。只有在 $P_\mu$ 行满秩时，

$$
P_\mu P_\mu^+ = I_{N_\mu}
$$

成立。

## 5. 矩阵乘、点乘和共轭转置

先解释代码中三个 MATLAB 操作：

```matlab
phi * phi_mu'
psi * psi_mu'
(phi * phi_mu') .* (psi * psi_mu')
```

### 5.1 矩阵乘 `*`

设

$$
\Phi(r,i)=\phi_i(r),
\qquad
\Phi_\mu(\mu,i)=\phi_i(r_\mu).
$$

MATLAB 中的 `'` 是共轭转置，因此

$$
\Phi_\mu' = \Phi_\mu^*.
$$

矩阵乘

```matlab
phi * phi_mu'
```

对应

$$
A = \Phi\Phi_\mu^*,
$$

其矩阵元为

$$
A(r,\mu)
=
\sum_i
\Phi(r,i)\overline{\Phi_\mu(\mu,i)}
=
\sum_i
\phi_i(r)\overline{\phi_i(r_\mu)}.
$$

同理，设

$$
\Psi(r,j)=\psi_j(r),
\qquad
\Psi_\mu(\mu,j)=\psi_j(r_\mu),
$$

则

```matlab
psi * psi_mu'
```

对应

$$
B = \Psi\Psi_\mu^*,
$$

其矩阵元为

$$
B(r,\mu)
=
\sum_j
\psi_j(r)\overline{\psi_j(r_\mu)}.
$$

### 5.2 点乘 `.*`

MATLAB 中的 `.*` 是逐元素乘法。也就是说

```matlab
c1 = A .* B;
```

对应

$$
c_1(r,\mu)=A(r,\mu)B(r,\mu).
$$

代入上面的 $A$ 和 $B$：

$$
c_1(r,\mu)
=
\left[
\sum_i
\phi_i(r)\overline{\phi_i(r_\mu)}
\right]
\left[
\sum_j
\psi_j(r)\overline{\psi_j(r_\mu)}
\right].
$$

### 5.3 为什么可以变成双重求和

这里用的是普通乘法的分配律，而不是矩阵乘法交换律。

对于两个独立求和指标 $i$ 和 $j$，

$$
\left(\sum_i a_i\right)
\left(\sum_j b_j\right)
=
\sum_i\sum_j a_i b_j.
$$

令

$$
a_i=\phi_i(r)\overline{\phi_i(r_\mu)},
\qquad
b_j=\psi_j(r)\overline{\psi_j(r_\mu)},
$$

则

$$
\begin{aligned}
c_1(r,\mu)
&=
\sum_i\sum_j
\phi_i(r)\overline{\phi_i(r_\mu)}
\psi_j(r)\overline{\psi_j(r_\mu)} \\
&=
\sum_{ij}
\phi_i(r)\psi_j(r)
\overline{\phi_i(r_\mu)}
\overline{\psi_j(r_\mu)}.
\end{aligned}
$$

又因为复数共轭满足

$$
\overline{ab}=\overline{a}\,\overline{b},
$$

所以

$$
\overline{\phi_i(r_\mu)}
\overline{\psi_j(r_\mu)}
=
\overline{\phi_i(r_\mu)\psi_j(r_\mu)}.
$$

于是

$$
c_1(r,\mu)
=
\sum_{ij}
\phi_i(r)\psi_j(r)
\overline{\phi_i(r_\mu)\psi_j(r_\mu)}.
$$

这正是乘积态矩阵 $P$ 与插值点乘积态矩阵 $P_\mu$ 的相关矩阵：

$$
c_1(r,\mu)=[PP_\mu^*](r,\mu).
$$

### 5.4 “交换次序”到底是什么

这里没有使用

$$
AB=BA
$$

这样的矩阵乘法交换。矩阵乘法一般不可交换。

真正使用的是以下事实：

$$
\sum_i\sum_j a_i b_j
=
\sum_j\sum_i a_i b_j,
$$

以及 $a_i$ 只依赖 $i$、$b_j$ 只依赖 $j$，因此双重求和可以因式分解为两个单独求和的乘积：

$$
\sum_{ij} a_i b_j
=
\left(\sum_i a_i\right)
\left(\sum_j b_j\right).
$$

所以代码并不是把矩阵乘的顺序交换了，而是利用了乘积态

$$
P(r,ij)=\phi_i(r)\psi_j(r)
$$

的可分离结构，把一个大矩阵乘法 $PP_\mu^*$ 写成两个小矩阵乘法的逐元素乘积：

$$
PP_\mu^*
=
(\Phi\Phi_\mu^*) \odot (\Psi\Psi_\mu^*),
$$

其中 $\odot$ 表示逐元素乘法，也就是 MATLAB 的 `.*`。

## 6. 代码中的 `c1`

直接构造 $P$ 会非常大，所以代码利用乘积态的结构避免显式形成 $P$。

定义

$$
\phi_\mu(\mu,i) = \phi_i(r_\mu),
\qquad
\psi_\mu(\mu,j) = \psi_j(r_\mu).
$$

则

$$
\begin{aligned}
[PP_\mu^*](r,\mu)
&= \sum_{ij} P(r,ij)\overline{P_\mu(\mu,ij)} \\
&= \sum_{ij}
\phi_i(r)\psi_j(r)
\overline{\phi_i(r_\mu)\psi_j(r_\mu)} \\
&= \sum_{ij}
\phi_i(r)\psi_j(r)
\overline{\phi_i(r_\mu)}
\overline{\psi_j(r_\mu)}.
\end{aligned}
$$

这个双重求和可以拆成两个单独求和的乘积：

$$
\begin{aligned}
[PP_\mu^*](r,\mu)
&=
\left[
\sum_i \phi_i(r)\overline{\phi_i(r_\mu)}
\right]
\left[
\sum_j \psi_j(r)\overline{\psi_j(r_\mu)}
\right].
\end{aligned}
$$

对应 MATLAB：

```matlab
phi * phi_mu'
psi * psi_mu'
```

因此

```matlab
c1 = (phi * phi_mu') .* (psi * psi_mu');
```

正好等于

$$
c_1 = PP_\mu^*.
$$

这里没有交换矩阵乘法次序；只是使用乘积态的可分离结构，把大矩阵乘法拆成两个小矩阵乘法和逐元素乘法。

## 7. 代码中的 `c2`

同理：

$$
\begin{aligned}
[P_\mu P_\mu^*](\mu,\nu)
&= \sum_{ij}
P_\mu(\mu,ij)
\overline{P_\mu(\nu,ij)} \\
&=
\left[
\sum_i
\phi_i(r_\mu)\overline{\phi_i(r_\nu)}
\right]
\left[
\sum_j
\psi_j(r_\mu)\overline{\psi_j(r_\nu)}
\right].
\end{aligned}
$$

对应 MATLAB：

```matlab
c2 = (phi_mu * phi_mu') .* (psi_mu * psi_mu');
```

所以

$$
c_2 = P_\mu P_\mu^*.
$$

## 8. 最终插值基

因此

```matlab
zeta_real = c1 / c2;
```

就是

$$
Z = PP_\mu^*(P_\mu P_\mu^*)^{-1}.
$$

也就是 ISDF 辅助基函数

$$
Z(r,\mu)=\zeta_\mu(r).
$$

代码随后做 FFT：

```matlab
zeta_grid = reshape(zeta_real(:, imu), fftgrid);
zeta_fft = fftn(zeta_grid) / ngrid;
zetag(:, imu) = zeta_fft(idx_q);
```

得到倒空间中的辅助基

$$
\zeta_\mu(G),
$$

用于构造 epsilon 或 sigma 的矩阵元。

## 9. 与柯西积分的关系

上面的 `c1 / c2` 只是在构造 ISDF 辅助基 $Z$，还不是柯西积分。

柯西积分用于计算低秩极化函数中的系数矩阵：

$$
C\Omega^{-1}C^*.
$$

其中

$$
C(\mu,vc)=\psi_v^*(r_\mu)\psi_c(r_\mu),
\qquad
\Omega_{vc,vc} = \epsilon_v-\epsilon_c.
$$

直接求和为

$$
\sum_{vc}
\frac{C(:,vc)C(:,vc)^*}{\epsilon_v-\epsilon_c}.
$$

论文第三章中用柯西积分把 occupied 和 unoccupied 的能量依赖解耦，从而避免显式高标度的 $v-c$ 双重耦合求和。

代码中对应：

```matlab
isdf_comega_cstar(...)
isdf_comega_cstar({left_component_1, ...}, {right_component_1, ...}, ...)
```

## 10. 输入参数怎么选

ISDF 相关输入分成三层，不应把它们混在一起理解。

第一层是总开关：

```matlab
eps.isdf.enable = true;
sig.isdf.enable = true;
```

只有打开后才启用 ISDF。

第二层是外层工作流：

```matlab
isdf.algorithm = 'matrix_elements';
```

表示在 epsilon 或 sigma 的共享主循环中，只用 ISDF 加速每一个波函数矩阵元。这个路径仍然显式遍历能带对。

```matlab
isdf.algorithm = 'reduced_basis';
```

表示使用 reduced ISDF 路径。静态 epsilon 中直接在 ISDF 辅助基空间构造低秩极化函数；静态 sigma 中使用对应的低秩屏蔽相互作用和 COHSEX 收缩。这个路径目前用于静态计算。

当前 reduced 路径同时支持多 K 点和多分量 spinor。epsilon 对每个 q 遍历不可约 k 点，并把同一星中各 k 点的辅助基按 G 空间对称映射后合并成块低秩极化表示；sigma 再把不可约 q 上的 screened W 映射到 full BZ，并按 q 星权重累加。spinor 各分量先组成

$$
P(r,ij)=\sum_s \psi_{i,s}^*(r)\psi_{j,s}(r),
$$

其分量间交叉项保留在 Gram 矩阵和 reduced 极化系数中。

epsilon 和 sigma 不再调用独立的 reduced-basis 工作流函数。两者分别只有一套共享主循环：

```matlab
% epsilon.m: q -> spin -> irreducible k
ops = epsilon_ops(ctx);

% sigma.m: spin -> band -> target k -> irreducible q
ops = sigma_ops(ctx);
```

`ctx.method` 在进入主循环前解析一次，`ops` 随后预选矩阵元构造、收缩、累加和 finalize 算子。循环内部不再解析 `isdf.enable` 或 `isdf.algorithm`，也没有 reduced-basis 提前 `return`。direct、matrix-elements 和 reduced-basis 因此共用循环、对称性映射、进度和结果存储逻辑，只在数值算子处不同。

可复用的 ISDF 数值接口位于 `+isdf` package：

```matlab
space = isdf.build_space(left_components, right_components, ...
    idx_q, fftgrid, isdf_options);
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver_options);
screened = isdf.screened_w(space, coulg, polar);
```

矩阵元、实空间分量和低秩屏蔽收缩分别使用 `isdf.matrix_elements`、`isdf.real_component`、`isdf.screened_kernel` 和 `isdf.screened_contract`。`isdf_comega_cstar` 仍只是 `isdf.polarizability` 内部的数值核；`cauchy` 仅表示 `reduced_solver` 的一种选择，不再出现在 epsilon/sigma 工作流函数名中。

第三层只在 `reduced_basis` 中生效：

```matlab
isdf.reduced_solver = 'cauchy';
```

用柯西积分近似计算 reduced 系数矩阵 $C\Omega^{-1}C^*$，目标是避免显式的 occupied-unoccupied 双重求和。

```matlab
isdf.reduced_solver = 'direct';
```

用显式 occupied-unoccupied 双重求和计算同一个 reduced 系数矩阵，主要用于验证和小体系对照。这个选项不使用柯西积分，因此不应该理解成一种 Cauchy method。

推荐写法为：

```matlab
eps.isdf.enable = true;
eps.isdf.algorithm = 'reduced_basis';
eps.isdf.reduced_solver = 'cauchy';

sig.isdf.enable = true;
sig.isdf.algorithm = 'reduced_basis';
sig.isdf.reduced_solver = 'cauchy';
```

插值点采样由 `isdf.sample_method` 控制。`sample_method = 'qrcp'` 是精确 QRCP 参考路径，需要完整乘积态矩阵 \(P\)，因此多分量 spinor 情况会显式构造 `products`。`sample_method = 'qrcp_randomized'`、`'default'` 和 `'kmeans'` 不显式构造完整 \(P\)：随机 QRCP 只形成随机压缩后的乘积矩阵，`kmeans` 只使用乘积态范数权重，随后都通过 Gram 矩阵和插值点上的 `product_mu` 构造 \(Z\)。因此正式计算中建议优先使用非 `qrcp` 采样，`qrcp` 主要用于小体系验证。

`eps.isdf.output` 决定 reduced-basis epsilon 返回哪种对象：

```matlab
eps.isdf.output = 'screened_w';    % 默认，只返回低秩 screened W
eps.isdf.output = 'full_inverse'; % 只返回完整 eps.inv，用于对照
eps.isdf.output = 'both';         % 同时返回二者，用于验证
```

`eps.isdf_screened_w` 是 reduced-basis epsilon 自动生成的低秩屏蔽相互作用，不是用户输入参数。它保存的是 ISDF 辅助基空间中的 screened \(W\) 对象，用于后续 reduced-basis sigma 收缩。默认 `screened_w` 模式不保存完整的 `eps.inv`，因为完整 \(G\)-space 逆矩阵会抵消低秩路径的内存优势。

自能计算会根据传入的 epsilon 结果自动选择路径：

- 如果 `eps.isdf_screened_w` 存在，则使用低秩 screened \(W\) 收缩；
- 如果没有 `eps.isdf_screened_w` 但有普通 `eps.inv`，则回退到完整 \(G\)-space 屏蔽矩阵收缩。

因此，`reduced_basis` 下可以明确区分三种用途：

```matlab
eps.isdf.output = 'screened_w';
```

这是默认正式低秩路线。epsilon 不展开完整的 `chi0_sum`，也不保存完整 `eps.inv`，而是在 ISDF reduced space 中构造 `eps.isdf_screened_w`。后续 reduced-basis sigma 会直接使用这个低秩 \(W\) 做收缩。

```matlab
eps.isdf.output = 'full_inverse';
```

这是对照路线。epsilon 先用 ISDF reduced basis 得到低秩极化表示，再展开为完整 \(G\)-space 的 `chi0_sum` 并直接求 `eps.inv`。它保留了 ISDF 低秩构造极化函数的部分，但最后的求逆仍是完整矩阵求逆，因此主要用于和普通 epsilon 或 matrix-elements 路径比较，不是推荐的加速路线。

```matlab
eps.isdf.output = 'both';
```

这是验证路线。同时保存 `eps.inv` 和 `eps.isdf_screened_w`，用于确认 full inverse 对照和低秩 screened \(W\) 路径的一致性。

## 11. epsilon 中的组合

ISDF+Cauchy 后，静态极化函数可以写成低秩形式：

$$
\chi_0(G,G')
\approx
\sum_{\mu\nu}
Z(G,\mu) S_{\mu\nu} Z(G',\nu),
$$

其中

$$
S = C\Omega^{-1}C^*.
$$

在代码的多分量 reduced-basis epsilon 路径中，`coeff` 不再默认展开成完整的 \(G\)-space `eps.inv`。它会和 ISDF 辅助基、Coulomb 核一起保存为低秩 screened \(W\)：

```matlab
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver_options);
screened_polar.coeff = conj(coeff) * fact;
eps.isdf_screened_w{iq, 1} = isdf.screened_w( ...
    screened_space, coulg(:), screened_polar);
```

其中：

- `vc_space.zeta_g` 是 $Z(G,\mu)$；
- `coeff` 是 $C\Omega^{-1}C^*$；
- `eps.isdf_screened_w` 保存 reduced \(W\)，供 reduced-basis sigma 使用；
- 外面的共轭来自当前代码中矩阵元约定：

```matlab
chi0_tmp = chi0_tmp + conj(gme_temp) * gme_temp.' * eden
```

## 12. 公式 (3.28) 的 SMW 推导

公式 (3.28) 的目的，是在 ISDF 低秩表示下避免显式构造和反演完整的介电矩阵。

先把原始 occupied-unoccupied 乘积态写成 ISDF 形式：

$$
P_{vc}(G,vc)
\approx
\sum_\mu Z(G,\mu) C(\mu,vc).
$$

其中 \(Z(G,\mu)\) 是 ISDF 辅助基组，\(C(\mu,vc)=\psi_v^*(r_\mu)\psi_c(r_\mu)\) 是乘积态在插值点上的系数。于是静态极化函数可以写成：

$$
\begin{aligned}
\chi_0
&\approx
P_{vc}\Omega^{-1}P_{vc}^* \\
&\approx
Z C\Omega^{-1}C^* Z^* \\
&=
Z S Z^*,
\end{aligned}
$$

这里：

- $Z$ 表示 ISDF 辅助基组，对应代码中的 `vc_space.zeta_g`；
- $C$ 表示插值点上的 occupied-unoccupied 乘积系数；
- $S$ 表示 reduced 系数矩阵，对应代码中的 `polar.coeff`；
- $V$ 是对角的 Coulomb 核，对应代码中的 `vcoul` 或 `V_g`；
- $*$ 表示共轭转置。

于是介电矩阵为

$$
\epsilon
= I - V\chi_0
\approx I - V Z S Z^*.
$$

直接求逆需要形成完整的 \(G\)-space 矩阵：

$$
\epsilon^{-1}
=
\left(I - V Z S Z^*\right)^{-1}.
$$

对 Sherman-Morrison-Woodbury 公式

$$
\left(A + U C V\right)^{-1}
=
A^{-1}
- A^{-1}U
\left(C^{-1} + V A^{-1}U\right)^{-1}
V A^{-1}
$$

取

$$
A = I,\qquad
U = -VZ,\qquad
C = S,\qquad
V = Z^*,
$$

则有

$$
\begin{aligned}
\epsilon^{-1}
&=
\left(I - V Z S Z^*\right)^{-1} \\
&=
I
- (-VZ)
\left(S^{-1} + Z^*(-VZ)\right)^{-1}
Z^* \\
&=
I
+ VZ
\left(S^{-1} - Z^*VZ\right)^{-1}
Z^*.
\end{aligned}
$$

这就是公式 (3.28) 的核心结构：

$$
\boxed{
\epsilon^{-1}
=
I
+ VZ
\left(S^{-1} - Z^*VZ\right)^{-1}
Z^*
}.
$$

它把完整 \(G\)-space 逆矩阵的问题，转化成 reduced space 中的小矩阵逆：

$$
K =
\left(S^{-1} - Z^*VZ\right)^{-1}.
$$

因此不需要显式保存完整的 \(\chi_0(G,G')\)、\(\epsilon(G,G')\) 和 \(\epsilon^{-1}(G,G')\)。

代码中的对应关系在 `isdf.screened_w`：

```matlab
vmat = vc_space.zeta_g' * (vcoul .* vc_space.zeta_g);
smw_denominator = inv(polar.coeff) - vmat;
k_mu = (eye(size(polar.coeff)) - polar.coeff * vmat) \ polar.coeff;
screened.k_mu = k_mu;
```

其中：

- `vc_space.zeta_g` 对应 \(Z\)，也就是 ISDF 辅助基组；
- `polar.coeff` 对应 \(S\)；
- `vmat` 对应 \(Z^*VZ\)；
- `smw_denominator` 对应 \(S^{-1} - Z^*VZ\)；
- `screened.k_mu` 对应公式 (3.30) 中的 \(K\) 矩阵。

代码没有直接写成

```matlab
inv(polar.coeff) - vmat
```

再求逆，而是使用等价的右解形式：

```matlab
k_mu = (I - S * vmat) \ S;
```

因为

$$
\left(S^{-1} - B\right)^{-1}
=
\left(I - SB\right)^{-1}S.
$$

这里 \(B=Z^*VZ\)。这种写法避免显式形成部分逆矩阵，更符合数值计算习惯。

在后续自能计算中，如果 `eps.isdf_screened_w` 存在，`sigma` 不再使用完整的 `eps.inv`，而是通过

```matlab
screened_kernel = isdf.screened_kernel(...);
aqs_eps_coul = isdf.screened_contract(...);
```

完成 reduced screened \(W\) 的投影和收缩。

## 13. 后续自能各部分的计算

由上一节得到的 reduced screened \(W\) 系数矩阵记为

$$
K =
\left(S^{-1} - Z_{vc}^* V Z_{vc}\right)^{-1}.
$$

这里 \(Z_{vc}\) 是 occupied-unoccupied 乘积空间的 ISDF 辅助基。根据公式 (3.28)，完整 \(G\)-space 中的屏蔽修正部分可以写成

$$
W_1
=
V Z_{vc} K Z_{vc}^* V.
$$

这里的 \(W_1\) 对应 \((\epsilon^{-1}-I)V\)，也就是 screened Coulomb 相互作用中相对裸 Coulomb \(V\) 的修正部分。因此静态 COHSEX 中通常分成三类量：

- 裸交换 \(E_X\)，只含 \(V\)；
- 屏蔽交换修正 \(E_{\mathrm{SEX}-X}\)，含 \(W_1\)；
- Coulomb-hole 项 \(E_{\mathrm{COH}}\)，含 \(W_1\) 并带 \(1/2\) 因子。

### 13.1 将 \(W_1\) 投影到自能乘积空间

自能计算并不直接需要完整的 \(W_1(G,G')\)。对给定目标态 \(n\)，需要的是目标乘积态和 \(W_1\) 的二次型。定义目标乘积空间

$$
P_{tn}(G,j) = \psi_n^*(G) \psi_j(G),
$$

其中 \(j\) 可以是占据态 \(v\)，也可以是求和用的全部能带 \(n''\)。对这个目标乘积空间再做 ISDF：

$$
P_{tn} \approx Z_{tn} C_{tn}.
$$

于是公式 (3.29) 中的 reduced kernel 可以写成

$$
\begin{aligned}
W_{1,tn}
&=
P_{tn}^* V Z_{vc} K Z_{vc}^* V P_{tn} \\
&\approx
C_{tn}^*
\left[
Z_{tn}^* V Z_{vc} K Z_{vc}^* V Z_{tn}
\right]
C_{tn}.
\end{aligned}
$$

方括号中的小矩阵就是代码中的 `screened_kernel`。它由

```matlab
screened_kernel = isdf.screened_kernel( ...
    screened, target_zeta, coulg_cutoff);
```

计算。其中：

- `screened.zeta_g` 对应 \(Z_{vc}\)；
- `screened.k_mu` 对应 \(K\)；
- `target_zeta` 对应 \(Z_{tn}\)；
- `coulg_cutoff` 对应 \(V\)。

`isdf.screened_kernel` 中的三步对应为：

```matlab
left_projector = target_zeta_g.' * (left_vcoul .* screened.zeta_g);
right_projector = screened.zeta_g' * (right_vcoul .* conj(target_zeta_g));
kernel = left_projector * screened.k_mu * right_projector;
```

也就是

$$
\mathrm{kernel}
=
Z_{tn}^T V Z_{vc}
K
Z_{vc}^* V \overline{Z_{tn}}.
$$

共轭和转置的具体位置来自当前代码对矩阵元 `aqs` 的约定；核心结构是把完整 \(G\)-space 的 \(W_1\) 投影成目标 ISDF 空间中的小矩阵。

### 13.2 对单个能带的 reduced 收缩

目标乘积态在插值点上的系数为

$$
c_j = C_{tn}(:,j).
$$

代码中对应

```matlab
coeff = sigma_space.product_mu(:, nn);
aqs_eps_coul = isdf.screened_contract(screened_kernel, coeff);
```

而 `isdf.screened_contract` 做的是

```matlab
value = coeff_left(:).' * kernel * conj(coeff_right(:));
```

数学上就是

$$
\langle P_{tn,j} | W_1 | P_{tn,j} \rangle
\approx
c_j^T \mathrm{kernel}\,\overline{c_j}.
$$

这个值在代码中记为 `aqs_eps_coul`。

### 13.3 裸交换 \(E_X\)

裸交换项不需要 screened \(W_1\)，只用裸 Coulomb \(V\)。代码先由 ISDF 还原目标乘积态在 \(G\)-space 的矩阵元：

```matlab
aqs_isdf = sigma_space.zeta_g * sigma_space.product_mu;
aqs_nocut = aqs_isdf(:, nn);
```

其中 `aqs_isdf(:, nn)` 对应 \(P_{tn,j}(G)\)。裸交换收缩为

```matlab
ax_loc = ax_loc - occ_kq(nn) * sum(abs(aqs_nocut).^2 .* coulg_nocut);
```

即

$$
E_X(n)
=
-
\sum_{v}^{\mathrm{occ}}
f_v
\langle P_{vn} | V | P_{vn} \rangle.
$$

这里 \(f_v\) 对应 `occ_kq(nn)`。只有占据态贡献裸交换。

### 13.4 屏蔽交换修正 \(E_{\mathrm{SEX}-X}\)

当 `eps.isdf_screened_w` 存在时，屏蔽交换修正使用 reduced \(W_1=(\epsilon^{-1}-I)V\)：

```matlab
if occ_kq(nn) > 0
    asx_loc = asx_loc - occ_kq(nn) * aqs_eps_coul;
end
```

对应

$$
E_{\mathrm{SEX}-X}(n)
=
-
\sum_{v}^{\mathrm{occ}}
f_v
\langle P_{vn} | W_1 | P_{vn} \rangle.
$$

最终代码把裸交换和屏蔽交换修正相加：

$$
E_{\mathrm{SEX}}
=
E_X + E_{\mathrm{SEX}-X}.
$$

这就是为什么代码同时保存 `ax_loc` 和 `asx_loc`：

```matlab
sig.sig = real(asx + ax + ach) * ryd;
```

### 13.5 Coulomb-hole 项 \(E_{\mathrm{COH}}\)

Coulomb-hole 项对所有求和能带贡献：

```matlab
ach_loc = ach_loc + aqs_eps_coul;
```

循环结束后代码写入

```matlab
ach(n_index, ik, ispin) = 0.5 * ach_loc;
```

对应

$$
E_{\mathrm{COH}}(n)
=
\frac{1}{2}
\sum_{n''}
\langle P_{n''n} | W_1 | P_{n''n} \rangle.
$$

这里的 \(1/2\) 是静态 COHSEX 中 Coulomb-hole 项的标准因子。

### 13.6 当前实现和论文公式的关系

论文公式 (3.29) 分别写出 \(P_{vn}\) 和 \(P_{nn}\) 两个目标空间：

$$
W_{1,vn}
=
P_{vn}^* V Z_{vc} K Z_{vc}^* V P_{vn},
$$

$$
W_{1,nn}
=
P_{nn}^* V Z_{vc} K Z_{vc}^* V P_{nn}.
$$

当前代码为了实现简洁，对给定目标态 \(n\) 构造一个包含 `1:nbands` 的目标乘积空间：

```matlab
right_components{ispinor} = isdf.real_component( ...
    wfnkq, fft.Nfft2, idx.kq, ispin, ispinor, 1:nbands);
sigma_space = isdf.build_space(left_components, right_components, ...
    idx.q, grid_size, isdf_options);
```

然后：

- 当 `occ_kq(nn) > 0` 时，这一列用于 \(P_{vn}\) 的裸交换和屏蔽交换；
- 对所有 `nn = 1:nbands`，这一列用于 \(P_{nn}\) 的 COH 求和。

也就是说，代码用一个统一的 `sigma_space` 同时覆盖 \(vn\) 和 \(nn\) 两类目标乘积态。数学结构仍然是公式 (3.29) 和 (3.31) 的 reduced 收缩，只是实现上没有拆成两个独立的 `vn_space` 和 `nn_space`。

如果没有 `eps.isdf_screened_w`，代码会回退到完整 \(G\)-space 的 `eps.inv`：

```matlab
eps_inv_i_coul = (eps_inv - eye(n_cutoff)) .* coulg_cutoff';
sigma_cohsex(...)
```

这时仍可以使用 ISDF 构造 `aqs_isdf`，但 screened \(W_1\) 不再通过公式 (3.29) 的 reduced kernel 计算。

## 14. 病态矩阵警告

如果 $P_\mu P_\mu^*$ 接近奇异，MATLAB 会在：

```matlab
zeta_real = c1 / c2;
```

处提示矩阵接近奇异。

这通常表示：

- 插值 rank 太高；
- 乘积态之间线性相关强；
- full-rank 验证时尤其容易出现；
- spinor 或多 k 点体系中也更常见。

当前代码使用 `isdf_stable_right_solve` 包装这个求解。默认仍保持直接右除，以维持和原始验证路径一致，但会静默 MATLAB 的 near-singular warning。若需要诊断，可设置：

```matlab
eps.isdf.warn_ill_conditioned = true;
sig.isdf.warn_ill_conditioned = true;
```

这样会显示更明确的 ISDF 条件数警告。
