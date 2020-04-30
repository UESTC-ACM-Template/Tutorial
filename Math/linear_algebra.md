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

### 基变换

考虑域$F$上的$n$维线性空间$V$，$m$维线性空间$W$和线性算子$T:V \to W$

设$B=\mat{v_1,v_2,\cdots,v_n}$为$V$的基，$C=\mat{w_1,w_2,\cdots,w_m}$为$W$的一组基

对于$V$中的向量$v$，其在$B$下的坐标为$x$，$W$中向量$w$,其在$C$下的坐标为$y$

$T$的矩阵表示为$A$，则
$$
 \begin{CD}...@>\partial_*>>H_{q}(A)@>i_*>> H_{q}(X) @>j_*>> H_{q}(X,A)@>\partial_*>>H_{q-1}(A)@>i_*>>... \\  & @V{f|_A}_* VV & @VVf_* V @VVf_* V@V{f|_A}_* VV\\  ...@>\partial_*>>H_{q}(B)@>i_*>> H_{q}(Y) @>j_*>> H_{q}(Y,B)@>\partial_*>>H_{q-1}(B)@>i_*>>...  \end{CD}\\
 \begin{CD}
 ...@> \partial_ 2>>\\
 \end{CD}
$$


### 子空间投影

## 特征值与特征向量

对域$\mathbb F$上的矩阵$A \in \mathbb F^{n \times n}$，若存在$\lambda \in \mathbb F,\vec x \in \mathbb F^{n \times 1}$使得$A\vec x=\lambda \vec x$，则$\lambda$被称为$A$的特征值，而$\vec x$为$A$的特征向量。

### 特征多项式

矩阵$A \in \mathbb F^{n \times n}$的特征多项式为：$f(\lambda) = |\lambda I-A|$

定理(Hamilton-Caylay)：$f(A)=0$

注：若$A$为$n$阶方阵，则其特征多项式$f(x)$的次数至多为$n$。
