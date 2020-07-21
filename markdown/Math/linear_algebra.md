$$
\def\brac#1{\left(#1\right)}
\def\mat#1{\left[\begin{array}{cl}#1\end{array}\right]}
$$

# 线性代数

## 行列式与线性方程组

### Cramer法则

对于非齐次线性方程组$A\vec x= \vec b$，其中$A=\{\vec{a}_1, \vec{a}_2, \cdots, \vec{a}_n\};A \in \R^{n\times n};\vec a_i, b \in \R^{n \times 1}; |A| \neq 0$。

设$A_i=\{\vec a_1, \cdots, \vec a_{i-1}, \vec b, \vec a_{i+1}, \cdots, \vec a_n\}$，则$\vec x$的第$i$个分量等于$\frac{|A_i|}{|A|}$。

## 高斯消元与逆矩阵

### 解线性方程组

### 求逆矩阵

### 求伴随阵

### 模质数意义下的行列式

高斯消元，$O(n^3)$。

```cpp
typedef vector<int> vec;
typedef vector<vec> mat;

int det(mat a) {
    const int n = a.size(); int res = 1;
    for (int i = 0; i != n; ++i) {
        for (int j = i + 1; j != n && !a[i][i]; ++j)
            if (a[j][i]) swap(a[i], a[j]), res = sub(0, res);
        if (!a[i][i]) return 0; else res = mul(res, a[i][i]);
        for (int j = i + 1, t = inv(a[i][i]); j != n; ++j)
            for (int k = i, w = mul(a[j][i], t); k != n; ++k)
                a[j][k] = sub(a[j][k], mul(w, a[i][k]));
    }
    return res;
}
```

### 模任意正整数意义下的行列式

高斯消元+辗转相除，$O(n^3 \log U)$。其中$U$是模数。

```cpp
typedef vector<int> vec;
typedef vector<vec> mat;

int det(mat a) {
    const int n = a.size(); int res = 1, f = 0;
    for (int i = 0; i != n; ++i) {
        for (int j = i + 1; j != n && !a[i][i]; ++j)
            if (a[j][i]) swap(a[i], a[j]), f ^= 1;
        if (!a[i][i]) return 0;
        for (int j = i + 1; j != n; swap(a[i], a[j]), f ^= 1, ++j)
            for (; a[i][i]; swap(a[i], a[j]), f ^= 1)
                for (int k = i, q = a[j][i] / a[i][i]; k != n; ++k)
                    a[j][k] = sub(a[j][k], mul(q, a[i][k]));
        res = mul(res, a[i][i]);
    }
    return f ? sub(0, res) : res;
}
```

### 实数行列式

高斯消元，$O(n^3)$

```cpp
typedef double dbl;
const dbl eps = 1e-10;
int sgn(dbl f) { return f < -eps ? -1 : f > eps; }

typedef vector<dbl> vec;
typedef vector<vec> mat;

dbl det(mat a) {
    const int n = a.size(); dbl res = 1;
    for (int i = 0; i != n; ++i) {
        for (int j = i + 1; j != n; ++j)
            if (abs(a[j][i]) > abs(a[i][i]))
                swap(a[i], a[j]), res = -res;
        if (!sgn(a[i][i])) return 0; else res *= a[i][i];
        for (int j = i + 1; j != n; ++j) {
            if (!sgn(a[j][i])) continue;
            double w = a[j][i] / a[i][i];
            for (int k = i; k != n; ++k)
                a[j][k] -= w * a[i][k];
        }
    }
    return res;
}
```

## 线性空间与线性算子

### 基

线性无关的极大向量组被称为线性空间的一组基。

命题：所有基大小相同。

证明：



命题：设$V$，$W$分别为$\mathbb F$上的$n,m$维线性空间，则对于$V$的一组基$B$，$W$的一组基$C$来说$T:V\to W$有唯一的矩阵表示$A \in \mathbb F^{m \times n}$。

证明：考虑域$F$上的$n$维线性空间$V$，$m$维线性空间$W$和线性算子$T:V \to W$

设$B=\mat{v_1,v_2,\cdots,v_n}$为$V$的基，$C=\mat{w_1,w_2,\cdots,w_m}$为$W$的一组基

对于$V$中的向量$v$，其在$B$下的坐标为$x$，$W$中向量$w$，其在$C$下的坐标为$y$，$T$在基$B,C$下的矩阵表示为$A$，设
$$
T(v_i)=\sum_{j=1}^mw_ja_{ji}
$$
则有
$$
T(\vec v)=T(B\vec x)=T\brac{\sum_{i=1}^n\vec v_ix_i}=\sum_{i=1}^nT(\vec v_ix_i)=\sum_{i=1}^nT(\vec v_i)x_i\\
=\sum_{i=1}^n\sum_{j=1}^ma_{ji}\vec w_jx_i=\sum_{j=1}^m\vec w_j\sum_{i=1}^na_{ji}x_i=\sum_{j=1}^m\vec w_jy_j
$$
由此有
$$
A=\{a_{ij}\},\vec y=A\vec x
$$

$$
\begin{CD}V@>T>>W\\@A{B}AA @AA{C}A\\\mathbb F^n @>A>>\mathbb F^m\end{CD}
$$

### 基变换

设$B=\{\vec v_1,\cdots,\vec v_n\}$为$V$的一组基，$B'=\{\vec v'_1,\cdots,\vec v'_n\}$为另一组基。

由基的定义，可将$v'_i$用$B$展开。即
$$
\vec v'_i=\sum_{j=1}^n\vec v_jp_{ji}
$$
写作矩阵形式即
$$
B'=\mat{\vec v'_1,\cdots,\vec v'_n}
=\mat{\displaystyle\sum_{j=1}^np_{ij}\vec v_j}
=\mat{\vec v_1, \cdots, \vec v_n}
\mat{
p_{11} & p_{12} & \cdots & p_{1n}\\
p_{21} & p_{22} & \cdots & p_{2n}\\
\vdots & \vdots & \ddots & \vdots \\
p_{n1} & p_{n2} & \cdots & p_{nn}
}
=BP
$$
考虑基变换对向量的影响：设向量$\vec v$在$B$中坐标为$\vec x$，在$B'$中坐标为$\vec x'$，则
$$
\vec v=B\vec x=B'\vec x'
$$
所以有
$$
\vec x'=P^{-1}\vec x
$$
设$C=\{\vec w_1, \cdots , \vec w_m\}$为$W$的一组基，$C'=\{w'_1, \cdots, w'_m\}$为另一组基，同上有$C'=CQ$。

考虑基变换对线性变换$T:v \to w$在基$B,C$下的矩阵的影响。

设$T(v)=w$，则因$v=Bx=B'x'.w=Cy=C'y'$有

$T(B'x')=C'y'$



### 子空间投影

定理：设子空间$W$的一组基为$A=\mat{a_1,a_2,\cdots,a_m} \in \R^{n\times m}$，则矩阵$P=A(A^TA)^{-1}A^T \in \R^{n \times n}$将任何$\R^n$中的向量投影至该子空间中。

证明：设$x$的投影为$x'$，$x'$在$A$下的坐标为$y$，则有
$$
A^T(Ay)=A^Tx
$$
由$x'=Ay$，有
$$
x'=Ay=A(A^TA)^{-1}A^Tx=Px
$$
命题：$A^TA$正定，因此$(A^TA)^{-1}$存在。

证明：对于任意$x \in \R^m \setminus \{0\}$，有
$$
x^T(A^TA)x=(Ax)^T(Ax)=||Ax||_2>0
$$

## 特征值与特征向量

对域$\mathbb F$上的矩阵$A \in \mathbb F^{n \times n}$，若存在$\lambda \in \mathbb F,\vec x \in \mathbb F^{n \times 1}$使得$A\vec x=\lambda \vec x$，则$\lambda$被称为$A$的特征值，而$\vec x$为$A$的特征向量。

### 特征多项式

矩阵$A \in \mathbb F^{n \times n}$的特征多项式为：$f(\lambda) = |\lambda I-A|$

定理(Hamilton-Caylay)：$f(A)=0$

注：若$A$为$n$阶方阵，则其特征多项式$f(x)$的次数至多为$n$。
