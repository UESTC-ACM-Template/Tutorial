# 线性代数

## 行列式与线性方程组

### Cramer法则

对于非齐次线性方程组$A\vec x= \vec b$，其中$A=\{\vec{a}_1, \vec{a}_2, \cdots, \vec{a}_n\};A \in \R^{n\times n};\vec a_i, b \in \R^{n \times 1}; |A| \neq 0$。

设$A_i=\{\vec a_1, \cdots, \vec a_{i-1}, \vec b, \vec a_{i+1}, \cdots, \vec a_n\}$，则$\vec x$的第$i$个分量等于$\frac{|A_i|}{|A|}$。

## 高斯消元与逆矩阵

### 解线性方程组

### 求逆矩阵

### 求伴随阵

### 实数行列式

### 模质数意义下的行列式

### 模任意整数意义下的行列式

## 线性空间

## 特征值与特征向量

对域$\mathbb F$上的矩阵$A \in \mathbb F^{n \times n}$，若存在$\lambda \in \mathbb F,\vec x \in \mathbb F^{n \times 1}$使得$A\vec x=\lambda \vec x$，则$\lambda$被称为$A$的特征值，而$\vec x$为$A$的特征向量。

### 特征多项式

矩阵$A \in \mathbb F^{n \times n}$的特征多项式为：$f(\lambda) = |\lambda I-A|$

定理(Hamilton-Caylay)：$f(A)=0$

注：若$A$为$n$阶方阵，则其特征多项式$f(x)$的次数至多为$n$。

## 线性递推

定义：数列$\{a_i\}$为线性递推数列当且仅当存在$k$与数列$c_1, c_2, \cdots, c_k$使得

$$\forall i \geq k,a_i=\sum_{j=1}^kc_ja_{i-j}$$

构造向量$\vec x$与矩阵$C$如下

$$\vec x_i = \left[
\begin{array}{cl}
a_{i+k-1}\\
a_{i+k-2}\\
\vdots\\
a_{i+1}\\
a_{i}
\end{array}
\right]
,
C=\left[
\begin{array}{lcl}
c_1 & c_2 & c_3 & \cdots & c_{k-1} & c_k\\
1 & 0 & 0 & \cdots & 0 & 0\\
0 & 1 & 0 & \cdots & 0 & 0\\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots\\
0 & 0 & 0 & \cdots & 0 & 0\\
0 & 0 & 0 & \cdots & 1 & 0
\end{array}
\right]$$

易得$\vec x_i=C\vec x_{i-1}=C^i\vec x_0$。

利用快速幂可在$O(k^3 \log n)$内求出$a_n$。

### 快速线性递推

观察到对于正整数$n$，若存在序列$\{b_j\}$使得$C^n=\sum_{j=0}^{k-1}b_jC^j$

将两边乘上$\vec x_0$：$C^n \vec x_0=\sum_{j=0}^{k-1}b_j C^j \vec x_0=\sum_{j=0}^{k-1}b_j \vec x_j$

取出第一行：$a_n=\sum_{j=0}^{k-1}b_ja_j$

则求出序列$\{b_j\}$后可在$O(k)$内求出$a_n$。

接下来考虑如何求序列$\{b_j\}$。

考虑$C$的特征多项式$f(x)$与多项式$g(x)=x^n$。

使用多项式除法将$g$除以$f$求出$q(x),r(x)$满足$g(x)=q(x)f(x)+r(x)$，则由Hamilton-Caylay定理有$C^n=g(C)=q(C)f(C)+r(C)=r(C)$，即$x^n$除以$f(x)$所得余多项式的各项系数即为所求的序列$\{b_j\}$。

于是使用多项式快速幂计算$C^n$的时候对$f$取模可在$O(k \log k \log n)$内求出$\{b_i\}$。

注：这里的快速幂不是用对数和指数直接对$x^n$取模的那个。

接下来考虑如何求$C$的特征多项式。

计算矩阵$C$的特征多项式$f(x)$，即行列式

$$|xI-C|=\left|
\begin{array}{cl}
x-c_1 & -c_2 & -c_3 & \cdots & -c_{k-1} & -c_k\\
-1 & x & 0 & \cdots & 0 & 0\\
0 & -1 & x & \cdots & 0 & 0\\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots\\
0 & 0 & 0 & \cdots & x & 0\\
0 & 0 & 0 & \cdots & -1 & x
\end{array}
\right|$$

对$i=k\cdots 2$，将第$i$列乘上$\frac 1 x$后加到第$i-1$列上后可得

$$f(x)=|xI-C|=x^{k-1}\left(x-c_1-\frac {c_2}x-\frac{c_3}{x^2}-\cdots\right)=x^k-\sum_{j=1}^kc_jx^{k-j}$$

### 最短线性递推式(Berlekamp-Massey)


