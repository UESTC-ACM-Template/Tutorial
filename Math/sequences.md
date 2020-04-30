$$
\def\brac#1#2#3{\left#1#2\right#3}
\def\floor#1{\left\lfloor #1 \right\rfloor}
\def\stirlingf#1#2{\genfrac{[}{]}{0pt}{0}{#1}{#2}}
\def\stirlings#1#2{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
\def\brac#1{\left(#1\right)}\def\mat#1{\left[\begin{array}{cl}#1\end{array}\right]}
$$

## 计数序列

### 等比数列与求和

$$
S_0(n)=\sum_{i=0}^nq^i=\frac{q^{n+1}-1}{q-1}\\
S_k(n)=\sum_{i=0}^ni^kq^i=\frac{1}{q-1}\brac{[}{(n+1)^kq^{n+1}-q\sum_{j=0}^{k-1}\binom kjS_j(n)}{]}
$$
证明：扰动法。
$$
S_k(n)=\sum_{i=0}^ni^kq^i=\sum_{i=1}^ni^kq^i=\sum_{i=0}^{n-1}i^kq^i+n^kq^n\\
q\sum_{i=0}^{n-1}(i+1)^kq^i=\sum_{i=1}^{n-1}i^kq^i+n^kq^n\\
\sum_{i=0}^{n-1}q^i\brac{[}{q(i+1)^k-i^k}{]}=n^kq^n\\
\sum_{i=0}^{n-1}q^i\brac{[}{q\sum_{j=0}^k\binom kji^j-i^k}{]}=n^kq^n\\
\sum_{i=0}^{n-1}q^i\brac{[}{q\sum_{j=0}^{k-1}\binom kji^j+(q-1)i^k}{]}=n^kq^n\\
q\sum_{j=0}^{k-1}\binom kj\sum_{i=0}^{n-1}i^jq^i+(q-1)\sum_{i=0}^{n-1}i^kq^i=n^kq^n\\
q\sum_{j=0}^{k-1}\binom kjS_j(n-1)+(q-1)S_k(n-1)=n^kq^n\\
S_k(n-1)=\frac{1}{q-1}\brac{[}{n^kq^n-q\sum_{j=0}^{k-1}\binom kjS_j(n-1)}{]}
$$

### 斐波那契数

定义：$F_0=0,F_1=1,F_n=F_{n-1}+F_{n-2}$

设转移矩阵$A=\mat{1 & 1\\1 & 0}$，则$\mat{F_n\\F_{n-1}}=A^{n-1}\mat{F_1\\F_0}$。

且因$A=\mat{F_2 & F_1 \\ F_1 & F_0}$有$A^{n}=\mat{F_{n+1} & F_n\\F_n & F_{n-1}}$。

---

性质：$\displaystyle \sum_{k=1}^nF_k=F_{n+2}-1 $

证明：$\displaystyle \sum_{k=1}^{n+1}F_k=\sum_{k=1}^{n}F_k+F_{n+1}=F_{n+2}-1+F_{n+1}=F_{n+3}-1$

---

性质：$\displaystyle \sum_{k=1}^nF_k^2=F_nF_{n+1}$

证明：$\displaystyle \sum_{k=1}^{n+1}=\sum_{k=1}^nF_k^2+F_{n+1}^2=F_nF_{n+1}+F_{n+1}^2=F_{n+1}(F_{n}+F_{n+1})=F_{n+1}F_{n+2}$

---

性质：$\displaystyle \sum_{k=0}^{n-2}F_{k}F_{k+3}=F_n^2-1$

证明：扰动平方和

$$
\sum_{k=1}^nF_k^2=1+\sum_{k=2}^nF_k^2=\sum_{k=1}^{n-1}F_k^2+F_n^2\\
\sum_{k=1}^{n-1}(F_{k+1}^2-F_k^2)=F_n^2-1\\
\sum_{k=1}^{n-1}F_{k+2}F_{k-1}=F_n^2-1\\
$$

---

$$
\mat{F_n\\F_{n-1}}=A^{n-k}\mat{F_k\\F_{k-1}}=\mat{F_{n-k+1}&F_{n-k}\\F_{n-k}&F_{n-k-1}}\mat{F_k\\F_{k-1}}\\
F_n=F_kF_{n-k+1}+F_{n-k}F_{k-1}
$$
换元即得
$$
F_{a+b}=F_aF_{b+1}+F_{b}F_{a-1}=F_{a}F_{b-1}+F_{b}F_{a+1}
$$

---

考虑$|A^n|=|A|^n$，于是有
$$
|A^n|=|A|^n=F_{n+1}F_{n-1}-F_n^2=(-1)^n|A^n|=|A|^n=F_{n+1}F_{n-1}-F_n^2=(-1)^n
$$

---

性质：$\displaystyle \sum_{k=1}^nF_kF_{k+1}$

---

---

（见线性递推-特征值分解部分）

$A$的特征多项式为$p(\lambda)=(1-\lambda)(-\lambda)-1$。令$p(\lambda)=0$解得$\displaystyle \lambda=\frac{\pm\sqrt 5+1}{2}$。

设其两个特征值为$\lambda_1, \lambda_2$，则有$\lambda_1\lambda_2=-1,\lambda_1+\lambda_2=1$
$$
P=\mat{\lambda_1 & \lambda_2\\1 & 1};P^{-1}=\frac 1{\sqrt 5}\mat{-1 & \lambda_2\\1 & -\lambda_1}
$$

$$
A^k
=\frac 1{\sqrt 5}\mat{\lambda_1 & \lambda_2\\1 & 1}\mat{\lambda_1 ^k& 0\\ 0 & \lambda_2^k}\mat{-1 & \lambda_2\\1 & -\lambda_1}
=\frac 1{\sqrt 5}\mat{\lambda_1^{k+1} &\lambda_2^{k+1}\\\lambda_1^k & \lambda_2^k}\mat{-1 & \lambda_2\\1 & -\lambda_1}\\
=\frac 1{\sqrt 5}\mat{\lambda_2^{k+1}-\lambda_1^{k+1} & \lambda_1^{k+1}\lambda_2-\lambda_2^{k+1}\lambda_1\\\lambda_2^k-\lambda_1^k & \lambda_1^{k}\lambda_2-\lambda_2^{k}\lambda_1}
=\frac 1{\sqrt 5}\mat{\lambda_2^{k+1}-\lambda_1^{k+1}&\lambda_2^k-\lambda_1^k\\\lambda_2^k-\lambda_1^k&\lambda_2^{k-1}-\lambda_1^{k-1}}
$$


$$
F_n=\frac 15(\lambda_1^n-\lambda_2^n)=\frac{1}{\sqrt 5}\left[\left(\frac{1+\sqrt 5}{2}\right)^n-\left(\frac{1-\sqrt 5}{2}\right)^n\right]\\
F_n^2=\frac 15\brac{\lambda_1^{2n}+\lambda_2^{2n}-2(-1)^n}
$$

---


$$
\sum_{k=1}^{n}{F_{2k-1}}=F_{2n}\\
\sum_{k=1}^{n}{F_{2k}}=F_{2n+1}-1\\
F_{n-1}F_{n+1}=F_n^2+(-1)^n\\
gcd(F_n,F_m)=F_{gcd(n,m)}\\
n|m \Leftrightarrow F(n)|F(m)
$$


OGF:
$$
F(x)=xF(x)+x^2F(x)+x\\
F(x)=\frac{x}{1-x-x^2}
$$

### 二项式系数

定义：

$$
\binom{n}{k}=\frac{n!}{k!(n-k)!}
$$

常见恒等式：

$$
\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}\\
\binom{n}{k}=\frac{n}{k}\binom{n-1}{k-1}=\frac{n-k+1}{k}\binom {n}{k-1}\\
\sum_{k=0}^n{\binom{m+k}{k}}=\binom{m+n+1}{n}\\
\sum_{k=0}^n{\binom{k}{m}}=\binom{n+1}{m+1}\\
\sum_{k=0}^{n}{\binom{m_1}{k}\binom{m_2}{n-k}}=\binom{m_1+m_2}{n}
$$

### 卡特兰数

性质：

$$
C_n=\frac{1}{n+1}\binom {2n}{n}=\prod_{k=2}^n{\left(1+\frac nk\right)}=\sum_{k=0}^{n-1}{C_{k}C_{n-k-1}}
$$

| $n=$  | 0    | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    | 9    | 10    |
| ----- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ----- |
| $C_n$ | 1    | 1    | 2    | 5    | 14   | 42   | 132  | 429  | 1430 | 4862 | 16796 |

OGF：
$$
F(x)=1+xF(x)^2=\frac{1-\sqrt{1-4x}}{2x}=\frac{2}{1+\sqrt{1-4x}}
$$

### 默慈金数

定义：$M_n$表示从$(0,0)$开始每次向右上或正右或右下走一格且不走到第四象限的情况下走到$(n,0)$的方案总数。

$M_n$表示在$n$个点的圆上画出数条不相交弦的全部方法的总数。

性质：

$$
M_n=\frac{(2(n-1)+3)M_{n-1}+(3(n-2)+3)M_{n-2}}{n+2}\\
M_n=M_{n-1}+\sum_{k=0}^{n-2}{M_kM_{n-2-k}}=\sum_{k=0}^{\left\lfloor\frac{n}{2}\right\rfloor}{\binom{n}{2k}C_k}
$$

OGF:
$$
F(x)=xF(x)+x^2F^2(x)=\frac{1-x-\sqrt{1-2x-3x^2}}{2x^2}=\frac{2}{1-x+\sqrt{1-2x-3x^2}}
$$

注：$C_k$为卡塔兰数。

| n=    | 0    | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    | 9    | 10   |
| ----- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| $M_n$ | 1    | 1    | 2    | 4    | 9    | 21   | 51   | 127  | 323  | 835  | 2188 |

### 那罗延数

定义：$N(n,k)$表示长度为$2n$的合法括号序列中有$k$对直接相邻的左右括号的方案数。

性质：

$$
N(n,k)=\frac{1}{n}\binom{n}{k-1}\binom{n}{k}
$$

$$
\sum_{k=1}^{n}{N(n,k)}=C_{n}
$$

注：$C_k$为卡塔兰数。

| n\k  | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| 1    | 1    |      |      |      |      |      |      |      |
| 2    | 1    | 1    |      |      |      |      |      |      |
| 3    | 1    | 3    | 1    |      |      |      |      |      |
| 4    | 1    | 6    | 6    | 1    |      |      |      |      |
| 5    | 1    | 10   | 20   | 10   | 1    |      |      |      |
| 6    | 1    | 15   | 50   | 50   | 15   | 1    |      |      |
| 7    | 1    | 21   | 105  | 175  | 105  | 21   | 1    |      |
| 8    | 1    | 28   | 196  | 490  | 490  | 196  | 28   | 1    |

### 第一类斯特林数

定义：$s(n,m)$表示将大小为$n$的集合划分成$m$个圆排列的方案数

考虑$n$所在排列（加入某个圆排列或单独作为一个新的圆排列）可得到递推式：

$$
s(n,m)=s(n-1,m-1)+(n-1)s(n-1,m)
$$

特别的，$s(n,0)=[n==0]$

注：每个圆排列可被视作一个置换，因而

$$
\sum_{k=0}^n{s(n,k)}=n!
$$

行的OGF：

由递推式有
$$
s_n(x)=xs_{n-1}(x)+(n-1)s_{n-1}=(x+n-1)s_{n-1}(x)=\prod_{i=0}^{n-1}{(x+i)}
$$

列的EGF:

单一圆排列的EGF为

$$
\sum_{k=1}^{\infty}{\frac{(k-1)!x^k}{k!}}=\sum_{k=1}^{\infty}{\frac{x^k}{k}}=-\ln(1-x)
$$

卷$m$次后除去$m!$消除圆排列之间的先后顺序即可得到第一类斯特林数的EGF

$$
s_m(x)=\frac{\left(-\ln(1-x)\right)^m}{m!}
$$

注：上式是无符号第一类斯特林数的EGF，带符号的为

$$
s_m(x)=\frac{\left(\ln(1+x)\right)^m}{m!}
$$

| n/m  | 0    | 1    | 2     | 3     | 4    | 5    | 6    | 7    | 8    |
| ---- | ---- | ---- | ----- | ----- | ---- | ---- | ---- | ---- | ---- |
| 0    | 1    |      |       |       |      |      |      |      |      |
| 1    | 0    | 1    |       |       |      |      |      |      |      |
| 2    | 0    | 1    | 1     |       |      |      |      |      |      |
| 3    | 0    | 2    | 3     | 1     |      |      |      |      |      |
| 4    | 0    | 6    | 11    | 6     | 1    |      |      |      |      |
| 5    | 0    | 24   | 50    | 35    | 10   | 1    |      |      |      |
| 6    | 0    | 120  | 274   | 225   | 85   | 15   | 1    |      |      |
| 7    | 0    | 720  | 1764  | 1624  | 735  | 175  | 21   | 1    |      |
| 8    | 0    | 5040 | 13068 | 13132 | 6769 | 1960 | 322  | 28   | 1    |

### 第二类斯特林数

定义：$S(n,m)$表示将大小为$n$的集合划分成$m$个非空集合的方案数

考虑$n$所在集合（加入$m$个集合中的一个或单独作为一个新的集合）可得到递推式：

$$
S(n,m)=S(n-1,m-1)+mS(n-1,m)
$$

考虑容斥，给每个集合编号后设$A_i$为第$i$个集合为空的放法，则所求为$\bar{A_i}$的交，即

$$
m!S(n,m)=m^n-\sum_{k=1}^{m}{(-1)^{k-1}\binom mk(m-k)^n}
$$

得

$$
S(n,m)=\frac{1}{m!}\sum_{k=0}^m{(-1)^{k-1}\binom mk (m-k)^n}
$$

性质：

考虑将$m$个物品分入$n$个盒子的方案数，枚举有$k$个有物品的盒子，可得

$$
m^n=\sum_{k=0}^mS(m,k)\binom nk k!
$$

二项式反演可得

$$
k!S(m,k)=\sum_{n=0}^k(-1)^{k-n}\binom kn n^m
$$
此式与上式等价。

行的OGF:

由通项公式得

$$
S_n(x)=\sum_{k=0}^{n}{\frac {x^k}{k!}\sum_{i=0}^k{(-1)^{i}\binom ki (k-i)^n}}=\sum_{k=0}^{n}{x^k\left(\sum_{i=0}^{k}{\frac{(-1)^{i}}{i!}\frac{(k-i)^n}{(k-i)!}}\right)}
$$

列的OGF:由递推式有：

$$
S_m(x)=xS_{m-1}(x)+mxS_{m}(x)
$$

得

$$
S_m(x)=\frac{x}{1-mx}S_{m-1}(x)=\frac{x^m}{\prod_{i=1}^{m}{1-ix}}
$$

| n/m  | 0    | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| 0    | 1    |      |      |      |      |      |      |      |      |
| 1    | 0    | 1    |      |      |      |      |      |      |      |
| 2    | 0    | 1    | 1    |      |      |      |      |      |      |
| 3    | 0    | 1    | 3    | 1    |      |      |      |      |      |
| 4    | 0    | 1    | 7    | 6    | 1    |      |      |      |      |
| 5    | 0    | 1    | 15   | 25   | 10   | 1    |      |      |      |
| 6    | 0    | 1    | 31   | 90   | 65   | 15   | 1    |      |      |
| 7    | 0    | 1    | 63   | 301  | 350  | 140  | 21   | 1    |      |
| 8    | 0    | 1    | 127  | 966  | 1701 | 1050 | 266  | 28   | 1    |

### 贝尔数

定义：贝尔数$B_n$表示将$n$化为至少一个非空集合的方案数。

性质：

$$Bell_n=\sum_{k=1}^nS(n,k)$$

注：$S(n,k)$为第二类斯特林数。

EGF: $Bell(x)=e^{e^x-1}$

| n=       | 0    | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    | 9     | 10     |
| -------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ----- | ------ |
| $Bell_n$ | 1    | 1    | 2    | 5    | 15   | 52   | 203  | 877  | 4140 | 21147 | 115975 |

### 伯努利数与等幂求和

定义（伯努利数）：
$$
B_0=1,\sum_{i=0}^{n}\binom{n+1}{i}B_i=0 \quad (n\geq 1)
$$
由其递归定义推导EGF$ B(x)$：
$$
\sum_{i=0}^{n-1}\binom ni B_i+B_n=B_n \quad (n \geq 2)\\
\sum_{i=0}^{n}\binom niB_i=B_n \quad (n \geq 2)\\
$$
当$n=1$时有
$$
\binom 10B_0+\binom 11B_1=\frac 12=-\frac 12 +1
$$
因此有
$$
e^xB(x)=B(x)+x\\
B(x)=\frac{x}{e^x-1}
$$

定义（等幂求和）：
$$
S_k(n)=\sum_{i=0}^{n-1}i^k
$$
设$F_n(x)[\frac{x^k}{k!}]=S_k(n)$
$$
F_n(x)=\sum_{i=0}^{\infty}S_k(n)\frac{x^k}{k!}=\sum_{k=0}^{\infty}\sum_{i=0}^{n-1}i^k\frac{x^k}{k!}=\sum_{i=0}^{n-1}\sum_{k=0}^{\infty}\frac{(xi)^k}{k!}=\sum_{i=0}^{n-1}e^{xi}=\frac{e^{xn}-1}{e^x-1}
$$
观察到
$$
F_n(x)=B(x)\frac{e^{xn}-1}{x}
$$
因为
$$
\frac{e^{xn}-1}{x}=\frac 1x\left(-1+\sum_{i=0}^\infty\frac{(nx)^i}{i!}\right)=\frac 1x\sum_{i=1}^\infty\frac{n^ix^i}{i!}
=\sum_{i=0}^\infty\frac{n^{i+1}}{i+1}\frac{x^i}{i!}
$$
所以
$$
S_k(n)=F_n(x)[x^k]=\sum_{i=0}^k\binom ki\frac{n^{k-i+1}}{k-i+1}B_i=\frac 1{k+1}\sum_{i=0}^k\binom{k+1}{i}n^{k-i+1}B_i
$$
是一个关于$n$的$k+1$次多项式。

| n=    | 0    | 1           | 2          | 3    | 4               | 5    | 6              | 7    | 8               | 9    | 10            |
| ----- | ---- | ----------- | ---------- | ---- | --------------- | ---- | -------------- | ---- | --------------- | ---- | ------------- |
| $B_n$ | $1$  | $-\frac 12$ | $\frac 16$ | $0$  | $-\frac{1}{30}$ | $0$  | $\frac{1}{42}$ | $0$  | $-\frac{1}{30}$ | $0$  | $\frac 5{66}$ |


### 五边形数与拆分数

拆分数：$a_n=|\{S|\forall b \in S, b > 0 \wedge b = \sum_{b \in S} b\}|$，即将$n$拆成数个正整数的方案。

不难列出拆分数的OGF:
$$
F(x)=\prod_{i=1}^{\infty}\sum_{j=0}^\infty{x^{ij}}=\prod_{i=1}^{\infty}{\frac{1}{1-x^i}}
$$

和五边形数相关的神奇函数：

$$
P(x)=\prod_{i=1}^\infty(1-x^i)=\sum_{i=0}^\infty{(-1)^ix^{\frac{1}{2}i(3i \pm 1)}}
$$

易得

$$
F(x)P(x)=1
$$

若限定每个数至多用$k$次，则有OGF:

$$
F_k(x)=\prod_{i=1}^{\infty}\sum_{j=0}^k{x^{ij}}=\prod_{i=1}^{\infty}{\frac{1-x^{(k+1)i}}{1-x^i}}=P(x^{k+1})F(x)
$$

注：$P(x)$前$n$项的系数只有$O(\sqrt{n})$个非零项，因此在预处理出$F(x)$的条件下可在$O(\sqrt{n})$内算出$F_k(x)$的第$n$项系数。

| n=    | 0    | 1           | 2          | 3    | 4               | 5    | 6              | 7    | 8               | 9    | 10            |
| ----- | ---- | ----------- | ---------- | ---- | --------------- | ---- | -------------- | ---- | --------------- | ---- | ------------- |
| $a_n$ | $1$  | $-\frac 12$ | $\frac 16$ | $0$  | $-\frac{1}{30}$ | $0$  | $\frac{1}{42}$ | $0$  | $-\frac{1}{30}$ | $0$  | $\frac 5{66}$ |

### 质数与质数幂