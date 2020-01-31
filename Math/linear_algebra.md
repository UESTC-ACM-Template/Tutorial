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

## 线性空间与线性算子

定义：线性空间是一个四元组$<\mathbb F, V, +, \cdot>$，其中$+:V \times V \rightarrow V, \cdot: \mathbb F \times V \rightarrow V$且满足8条公理：

第一组（群公理）：

结合律：$\forall \vec x, \vec y, \vec z \in V,\vec x + (\vec y + \vec z) = (\vec x + \vec y) + \vec z$

单位元存在性：$\exist \vec 0 \in V \forall \vec x \in V,\vec x + \vec 0 = \vec x$

逆元存在性：$\forall \vec x \in V \exists \vec y \in V,\vec x + \vec y = \vec 0$

交换律：$\forall \vec x, \vec y \in V, \vec x + \vec y = \vec y + \vec x$

即$<V,+>$是一个交换群。

第二组（线性公理）：

结合律：$\forall k_1, k_2 \in \mathbb F\forall \vec x \in V,k_1(k_2 \vec x)=(k_1 k_2)\vec x$

幺元存在性：$\forall \vec x \in V,1 \vec x=\vec x$

向量分配律：$\forall k_1, k_2 \in \mathbb F \forall \vec x \in V,(k_1 + k_2) \vec x =k_1 \vec x + k_2 \vec x$

标量分配律：$\forall k \in \mathbb F \forall \vec x_1, \vec x_2 \in V,k(\vec x_1 + \vec x_2)=k\vec x_1 + k\vec x_2$

定义：线性算子是线性空间之间的函数

## 特征值与特征向量

对域$\mathbb F$上的矩阵$A \in \mathbb F^{n \times n}$，若存在$\lambda \in \mathbb F,\vec x \in \mathbb F^{n \times 1}$使得$A\vec x=\lambda \vec x$，则$\lambda$被称为$A$的特征值，而$\vec x$为$A$的特征向量。

### 特征多项式

矩阵$A \in \mathbb F^{n \times n}$的特征多项式为：$f(\lambda) = |\lambda I-A|$

定理(Hamilton-Caylay)：$f(A)=0$

注：若$A$为$n$阶方阵，则其特征多项式$f(x)$的次数至多为$n$。

