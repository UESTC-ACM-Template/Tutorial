# 最优化

$\def\brac#1{\left(#1\right)}$前半部分为连续优化，后半部分为离散优化。

## 向量微积分

$\def\bs{\boldsymbol}\def\par#1#2{\frac{\part #1}{\part #2}}$以下用$\bs A, \bs B$指代矩阵，$\bs x, \bs y$指代列向量，$f, g, k,l$指代标量。

定义（向量对标量求导）：
$$
\par{\bs x}{k}=\brac{\par{x_1}{k},\par{x_2}{k}, \cdots \par{x_n}{k}}^T
$$
定义（矩阵对标量求导）：
$$
\par{\bs A}{k}=\left[
\begin{array}{cl}
\par{a_{11}}k & \par{a_{12}}k & \cdots & \par{a_{1m}}k\\
\par{a_{21}}k & \par{a_{22}}k & \cdots & \par{a_{2m}}k\\
\vdots & \vdots & \ddots & \vdots\\
\par{a_{n1}}k & \par{a_{n2}}k & \cdots & \par{a_{1m}}k\\
\end{array}
\right]
$$
定义（标量对向量求导）：
$$
\par{f}{\bs x}=\brac{\par{f}{x_1},\par{f}{x_2}, \cdots \par{f}{x_n}}^T
$$
例：$\displaystyle \par{}{\bs x}{\bs x}^T\bs x= \brac{\par{{\bs x}^T\bs x}{x_1},\par{{\bs x}^T\bs x}{x_2}, \cdots \par{{\bs x}^T\bs x}{x_n}}^T=\brac{2x_1,2x_2, \cdots, 2x_n}^T=2 \bs x$

例：$\displaystyle \par {}{\bs x} \bs y^T\bs x=\brac{y_i}^T=\bs y$

例：
$$
\displaystyle \par{}{\bs x} {\bs x}^T \bs A \bs x=\par{}{\bs x} \brac{\sum_{i=1}^n\sum_{j=1}^na_{ij}x_ix_j}^T=\brac{\par{}{x_k}\sum_{i=1}^n\sum_{j=1}^na_{ij}x_ix_j}^T\\
=\brac{2a_{kk}x_k+\sum_{i=1 \wedge i \neq k}^na_{ik}x_i+\sum_{j=1 \wedge j \neq k}^na_{kj}x_j}^T=\brac{\sum_{i=1}^na_{ik}x_k+\sum_{j=1}^nx_ka_{kj}}^T\\=\bs A \bs x+(\bs x^T \bs A)^T=(\bs A+\bs A^T)\bs x
$$
当$\bs A = {\bs A}^T$时，有$\par{}{\bs x} {\bs x}^T \bs A \bs x=2\bs A \bs x$

定义（向量对向量求导）：
$$
\par {\bs y}{\bs x}=\left[
\begin{array}{cl}
\par{y_1}{x_1} & \par{y_2}{x_1} & \cdots & \par{y_m}{x_1}\\
\par{y_1}{x_2} & \par{y_2}{x_2} & \cdots & \par{y_m}{x_2}\\
\vdots & \vdots & \ddots & \vdots\\
\par{y_1}{x_n} & \par{y_2}{x_n} & \cdots & \par{y_m}{x_n}\\
\end{array}
\right]
$$
例：$\displaystyle \par{\bs A^T\bs x}{\bs x}=\left[
\begin{array}{cl}
\par{}{x_1}\sum_{k=1}^na_{k1}x_k & \par{}{x_1}\sum_{k=1}^na_{k2}x_k & \cdots & \par{}{x_1}\sum_{k=1}^na_{km}x_k\\
\par{}{x_2}\sum_{k=1}^na_{k1}x_k & \par{}{x_2}\sum_{k=1}^na_{k2}x_k & \cdots & \par{}{x_2}\sum_{k=1}^na_{km}x_k\\
\vdots & \vdots & \ddots & \vdots\\
\par{}{x_n}\sum_{k=1}^na_{k1}x_k & \par{}{x_n}\sum_{k=1}^na_{k2}x_k & \cdots & \par{}{x_n}\sum_{k=1}^na_{km}x_k\\
\end{array}
\right]=\bs A$

当自变量确定时（如为$\bs x$），一般将$\displaystyle \par {f}{\bs x}, \par{\bs g}{\bs x}$记作$\nabla f, \nabla \bs g$。即$f,\bs g$的梯度。

## 一维单峰函数

## 一阶必要条件与拉格朗日对偶

给定一阶可微的函数$f : \mathbb{R} ^n \rightarrow \mathbb{R}$，若$\bs x_0$使得$f$取得极值，则必有$\nabla f(\bs x_0)=0$。



给定一阶可微的函数$f : \mathbb{R} ^n \rightarrow \mathbb{R}$和约束函数$\boldsymbol g : \mathbb{R}^n\rightarrow\mathbb{R}^{k},\boldsymbol h : \mathbb{R}^n\rightarrow\mathbb{R}^{m}$，对于带约束优化问题
$$
\begin{array}{lcl}
\min & &f(\boldsymbol x)\\
s.t.& &\boldsymbol h(\boldsymbol x) = 0\\
& &\boldsymbol g(\boldsymbol x) \leq 0
\end{array}
$$

可以构造一个奥妙重重的函数来把约束塞进目标函数里，即转化为无约束问题。

定义（带约束优化问题的拉格朗日函数）：
$$
\mathscr L(\bs x, \bs \lambda, \bs \mu)=f(\bs x)-\bs \lambda^T \bs h(\bs x)-\bs \mu^T\bs g(\bs x)\\
$$
里面的$\bs \lambda$和$\bs \mu$被称为拉格朗日乘子（后者为KKT乘子）。

先考虑只有等式约束的情况，观察在可行域内的极值点的特殊性质。

注意到任意一个满足约束的极小值点$\bs x_0$附近，朝着可行域内的任意一个方向走都不能让$f$减小（否则就不是极小点了）。

可行域内的方向即$\bs h$的每个分量变化率均为$0$的方向，即属于每个等式约束在$\bs x_0$处切空间的交。

让$f$保持不变的方向即$f(\bs x)=f(\bs x_0)$在点$\bs x_0$处切空间中的任意一个方向。

超曲面在某点的切空间可视作垂直于其在某点梯度的向量组成的空间，即由梯度的倍数形成的线性空间的补。

注：这里可以想象$\R^3$中有一个球形约束$h_1(\bs x)=||x||^2-r=0$，和一个平面约束$h_2(\bs x)=\bs w^T\bs x=0$ 。这两个约束的交（即可行域）是一个圆环，且对其上任意一点，$h_1$的切空间是球面上该点的切点，$h_2$的切空间即为该平面，其交即是圆环过此点的切线。

因为$\nabla \bs h(\bs x_0)$的列空间的补即每个等式约束在该点切空间的交，所以若$\nabla f(\bs x_0)$可表示成$\nabla \bs h(\bs x_0)$的列的线性组合则$\nabla f(\bs x_0)$的补包含$\nabla \bs h(\bs x_0)$的补。即每个等式约束在$\bs x_0$处切空间的交属于$f(\bs x)=f(\bs x_0)$在$\bs x_0$处的的切空间。

所以$\nabla f(\bs x_0)=\bs \lambda^T\nabla \bs h(\bs x_0)$等价于在$\bs x_0$附近，朝着可行域内的任意一个方向走都不能使$f$减小。

对于不等式约束$\bs g$，在最小值点取到等号的约束转化为等式约束。但其梯度方向应与$\nabla f(\bs x_0)$相反，即$\bs \mu \geq 0$，朝着可行域内方向走必须使$f$增大。考虑到在最小值点可能只会有$\bs g$的某个子集取到等号，而未取到等号的不等式约束在该点梯度向量的系数应该是0，所以在该点有
$$
\nabla f(\bs x_0)=\bs \lambda^T\nabla \bs h(\bs x_0)+\bs \mu^T\nabla \bs g(\bs x_0)\\
\bs \mu^T \bs g(\bs x)=0
$$
注意到第一个式子等价于
$$
\par{\mathscr L}{\bs x}\bigg|_{\bs x=\bs x_0}=0
$$
最终形态如下（Karush–Kuhn–Tucker条件）：

给定$f,\bs g, \bs h$，$f$在可行域中的极值点满足
$$
\par{\mathscr L}{\bs x}=\nabla f(\bs x)-\bs \lambda^T\nabla \bs h(\bs x)-\bs \mu^T\bs g(\bs x)=\bs 0\\
\bs h(\bs x)=\bs 0\\
\bs g(\bs x) \leq \bs 0\\
\bs \mu \geq 0\\
\bs \mu^T\bs g(\bs x)= \bs 0\\
$$
现在已经证明拉格朗日函数的梯度能给出和前面一样的必要条件，接下来就可以利用拉格朗日函数构造一个对应原问题的无约束优化问题了。

设
$$
\displaystyle \theta_{primal}(\bs x)=\max_{\lambda,\mu \geq 0} \mathscr L(\bs x, \bs \lambda, \bs \mu)
$$
则
$$
\min_{\bs x}\theta_{primal}(\bs x)=\min_{\bs x\in D}f(\bs x)
$$
注：$D$为可行域。

证明：对于不在可行域中的任意一个点$\bs x$，存在一个$h_i$或$g_i$使得$h_i(\bs x) \neq 0$或$g_i(\bs x) \geq 0$，因而改变对应的乘子分量即可使$\mathscr L$增大，且没有上限。所以$\theta_{primal}$在可行域外的取值是正无穷，因而其最小值等于$f$在可行域中的最小值。

观察到交换$\min$和$\max$的顺序不会改变答案，所以定义
$$
\theta_{dual}(\bs \lambda, \bs \mu)=\min_{\bs x} \mathscr L(\bs x, \bs \lambda, \bs \mu)
$$
有对偶问题
$$
\max_{\bs \lambda, \bs \mu \geq 0} \theta_{dual}(\bs \lambda, \bs \mu)=\min_{\bs x}\theta_{primal}(\bs x)=\min_{\bs x\in D}f(\bs x)
$$
对偶问题有时候会比原问题好解。

例：

支持向量机(SVM)是一个解决二分类问题的模型。其基本思想是使用一个超平面将两类数据点分开。

这里考虑SVM的训练过程，训练的目的即是使两类点到超平面的距离尽可能大。

第$i$个数据点为$\bs x_i \in \R^n$，标签为$y_i \in \{-1, 1\}$。超平面方程为${\bs w}^T \bs x+b=0$，$\bs w$相当于超平面的法向量。

预测函数为$y=\bs w^T\bs x+b$。约束即使得训练集中不存在使得$|\bs w^T \bs x+b| <1$的点。

恰好满足$|\bs w^T \bs x_i+b|=1$的点即是支持向量，即离超平面最近的那些点。

即若第$i$个数据点是支持向量，则有$y_i(\bs w^T \bs x_i+b)=1$。

接下来考虑解出支持向量到超平面的距离。

$\bs x_i$在超平面上的投影$\bs p$满足：$\bs w^T(\bs x_i+k \bs w)+b=0$。

展开：$\bs w^T \bs x_i+b+k\bs w^T\bs w=0$

代入：$y_i+k \bs w^T\bs w=0$

解得$\displaystyle k=\frac{y_i}{\bs w^T\bs w}$。即$\bs x_i$到超平面的距离的平方为$||k\bs w||^2=\displaystyle k^2\bs w^T\bs w=\frac{y_i^2}{\bs w^T\bs w}=\frac{1}{||\bs w||^2}$。

得到优化问题
$$
\max \frac{1}{||w||^2}\\
s.t. \quad \forall i, y_i(\bs w^T\bs x_i+b) \geq 1
$$
等价于解
$$
\min \frac 12 \bs w^T\bs w\\
s.t. \quad \forall i, 1-y_i(\bs w^T\bs x_i+b) \leq 0
$$

设
$$
g_i(\bs w)=1-y_i(\bs w^T\bs x_i+b)\\
$$


由KKT条件有
$$
\mathscr L= \frac 12 \bs w^T\bs w+ \sum_{i=1}^k\mu_ig_i=\frac 12 \bs w^T\bs w+ \sum_{i=1}^k\mu_i\brac{1-y_i(\bs w^T\bs x_i+b)}\\
s.t. \\
\par {\mathscr L}{\bs w}=\bs w-\par{}{\bs w}\brac{\sum_{i=1}^k\mu_i-\sum_{i=1}^k\mu_iy_i\bs w^T\bs x_i-\sum_{i=1}^k\mu_iy_ib}=\bs w-\sum_{i=1}^k\mu_iy_i\bs x_i=0\\
\par {\mathscr L}{b}=-\sum_{i=1}^k{\mu_iy_i}=0\\
\quad \mu_i \geq 0\\
1-y_i(\bs w^T\bs x_i+b) \leq 0
$$
代回$\mathscr L$有
$$
\mathscr L=\frac 12 \sum_{i=1}^m\sum_{j=1}^m\mu_i\mu_jy_iy_j\bs x_i^T\bs x_j+\sum_{i=1}^m\mu_i-\sum_{i=1}^m\mu_iy_i\brac{\sum_{j=1}^m\mu_jy_j\bs x_j^T}\bs x^i\\
=\sum_{i=1}^m\mu_i-\frac 12 \sum_{i=1}^m\sum_{j=1}^m\mu_i\mu_jy_iy_j\bs x_i^T\bs x_j
$$

由拉格朗日对偶可以丢掉不等式约束，原问题答案为
$$
\min _{\bs w} \max_{\bs \mu \geq 0} \mathscr L
$$
前式中$\bs w$已被消去，因此
$$
\max_{\bs \mu \geq 0} \mathscr L=\max_{\bs \mu \geq 0}\sum_{i=1}^m\mu_i-\frac 12 \sum_{i=1}^m\sum_{j=1}^m\mu_i\mu_jy_iy_j\bs x_i^T\bs x_j
$$


## 线性规划

### 单纯形法

### 对偶

## 最短路

## 最大流

## 费用流

## 拟阵

