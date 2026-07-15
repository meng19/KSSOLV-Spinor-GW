# ISDF 插值基构造与 `c1 / c2` 推导

本文说明 `src/GW/ISDF/isdf_kernelg_current_fft.m` 中

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
isdf_spinor_comega_cstar(...)
```

## 10. epsilon 中的组合

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

在代码的 spinor Cauchy epsilon 路径中对应：

```matlab
[coeff, info] = isdf_spinor_comega_cstar(...);
chi0_sum = chi0_sum + conj(vc_space.zeta_g) * conj(coeff) * vc_space.zeta_g.';
```

其中：

- `vc_space.zeta_g` 是 $Z(G,\mu)$；
- `coeff` 是 $C\Omega^{-1}C^*$；
- 外面的共轭来自当前代码中 `chi0` 的矩阵元约定：

```matlab
chi0_tmp = chi0_tmp + conj(gme_temp) * gme_temp.' * eden
```

## 11. 病态矩阵警告

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
