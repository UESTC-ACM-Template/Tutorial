$$
\def\brac#1{\left(#1\right)}
\def\floor#1{\left\lfloor #1 \right\rfloor}
\def\stirlingf#1#2{\genfrac{[}{]}{0pt}{0}{#1}{#2}}
\def\stirlings#1#2{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
$$

# 生成函数

## 幂级数

$$
\frac 1 {1-x}=\sum_{i=0}^{+\infty}x^i=1+x+x^2+\cdots
$$

$$
\frac x{(1-x)^2}=\sum_{i=0}^{+\infty}nx^n=x+2x^2+3x^3+\cdots
$$

$$
\frac {x(x+1)}{(1-x)^3}=\sum_{i=0}^{+\infty}n^2x^n=x+4x^2+9x^3+\cdots
$$

$$
\exp x=\sum_{i=0}^{+\infty} \frac {x^i}{i!}=1+x+\frac {x^2}{2}+\frac {x^3}{6}+\cdots
$$

$$
-\log {(1-x)}=\sum_{i=1}^{+\infty}\frac {x^n}n=x+\frac x2 + \frac {x^3} 3+\cdots
$$

$$
(1+x)^n=\sum_{i=0}^{n}\binom ni x^i
$$

### 多元幂级数

$$
\frac {1}{1-(1+x)y}=\sum_{i=0}^{+\infty}(1+x)^iy^i=\sum_{i=0}^{+\infty}\sum_{j=0}^{i} \binom ij x^jy^i
$$

## 常见变换

$$
f'(x)=\sum_{i=0}^{+\infty}(i+1)a_{i+1}x^i=a_1+2a_2x+3a_3x+\cdots
$$

第$i$项乘$i$：
$$
(Tf)(x)=xf'(x)=\sum_{i=0}^{\infty}ia_ix^i
$$
设$\displaystyle (T^{(n)}f)(x)=\sum_{k=0}^{n}c_{n,k}x^kf^{(k)}(x)$，则 
$$
(T^{(n+1)}f)(x)=x(T^{(n)}f)'(x)=x\left[\sum_{k=0}^nc_{n,k}(x^kf^{(k+1)}(x)+kx^{k-1}f^{(k)}(x)\right]\\
=\sum_{k=0}^nc_{n,k}(x^{k+1}f^{(k+1)}(x)+kx^{k}f^{(k)}(x))\\

=\sum_{k=0}^nc_{n,k}kx^{k}f^{(k)}(x)+\sum_{k=1}^{n+1}c_{n,k-1}x^{k}f^{(k)}(x)
$$
考虑到$c_{n,n+1}=0$，因此有
$$
=\sum_{k=0}^{n+1}(kc_{n,k}+c_{n,{k-1}})x^{k}f^{(k)}(x)=\sum_{k=0}^{n+1}c_{n+1,k}x^{k}f^{(k)}(x)
$$
即$c_{n,k}=kc_{n-1,k}+c_{n-1,k-1}$。结合初始条件有$\displaystyle c_{n,k}=\stirlings nk$。即
$$
(T^{(n)}f)(x)=\sum_{k=0}^{n}\stirlings nk x^kf^{(k)}(x)
$$
例：令$f(x)=(1-x)^{-1}$，即得
$$
f^{(k)}(x)=\sum_{i=k}^{\infty}i^{\underline k}x^{i-k}\\
(T^{(k)}f)(x)=\sum_{i=0}^{\infty}i^kx^i=\sum_{i=0}^k \stirlings ki  x^i\sum_{j=i}^{\infty}j^{\underline i}x^{j-i}=\sum_{i=0}^k \stirlings ki \sum_{j=i}^{\infty}j^{\underline i}x^{j}\\
=\sum_{j=0}^{\infty}x^j\sum_{i=0}^k\stirlings kij^{\underline i}
=\sum_{j=0}^{\infty}x^j\sum_{i=0}^k\stirlings ki \binom jii!
$$
因此有
$$
j^k=\sum_{i=0}^k\stirlings ki \binom jii!
$$
卷上一个$(1+x)^{-1}$即得
$$
\sum_{j=0}^nj^k=\sum_{j=0}^n\sum_{i=0}^k\stirlings ki \binom jii!=\sum_{i=0}^k\stirlings kii!\sum_{j=0}^n\binom ji=\sum_{i=0}^k\stirlings kii!\binom{n+1}{i+1}
$$


## 一般生成函数(OGF)

$$
f(x)=\sum_{i=0}^{+\infty}a_i x^i
$$

一般生成函数一般用于组合的计数（位置不可区分，或无标号）

### OGF卷积的组合意义

考虑$f$和$g$的卷积：

$$
h(x)=f(x)g(x)=\sum_{i=0}^{+\infty}\left(\sum_{j=0}^ia_jb_{i-j}\right)x^i=\sum_{i=0}^{+\infty}c_i x^i
$$

可以理解为：有两类物品组合。第一类中大小为$i$的有$a_i$种，第二类中大小为$i$的有$b_i$种。则分别从第一类和第二类中分别选出两个组合将其合并，大小为$i$的有$c_i$种。

### 欧拉变换

大小为$i$的组合有$a_i$种，现在计算总大小为$i$的由一个或多个组合组成的组合方案数量。

因为组合之间不可区分，所以可以枚举大小为$i$的组合中每一个用了多少。

即

$$
\mathscr E(f)(x)=\left(1+x+x^2+\cdots\right)^{a_1}\left(1+x^2+x^4+\cdots\right)^{a_2}\cdots=\prod_{i=1}^{+\infty}\left(\frac{1}{1-x^i}\right)^{a_i}
$$
$$
= \exp\left(\sum_{i=0}^{+\infty} a_i \ln \frac{1}{1-x^i}\right)=\exp\left(\sum_{i=0}^{+\infty} -a_i \ln{(1-x^i)}\right)
$$
$$
= \exp\left(\sum_{i=0}^{+\infty} a_i \sum_{j=1}^{+\infty}\frac{x^{ij}}{j}\right)
= \exp\left(\sum_{j=1}^{+\infty}\frac{1}{j}\sum_{i=0}^{+\infty}a_ix^{ij}\right)=\exp\left(\sum_{j=1}^{+\infty}\frac {f(x^j)}{j}\right)
$$

## 指数生成函数(EGF)

$$
f(x)=\sum_{i=0}^{+\infty}a_i \frac{x^i}{i!}
$$

指数生成函数一般用于排列的计数（位置可区分，或有标号）

### EGF卷积的组合意义

考虑$f$和$g$的卷积：

$$
h(x)=f(x)g(x)=\sum_{i=0}^{+\infty}\left(\sum_{j=0}^i\frac{a_j}{j!}\frac{b_{i-j}}{(i-j)!}\right)x^i=\sum_{i=0}^{+\infty}\left(\sum_{j=0}^i\binom ija_jb_{i-j}\right)\frac{x^i}{i!}=\sum_{i=0}^{+\infty}c_i x^i
$$

可以理解为：有两类物品排列。第一类中大小为$i$的有$a_i$种，第二类中大小为$i$的有$b_i$种。则分别从第一类和第二类中分别选出两个排列将其合并，大小为$i$的有$c_i$种。

注：先钦定来自第一类的物品放在哪些位置，来自第二类的物品放在剩下那些位置，因此有个二项式系数。

### $\exp f(x)$的组合意义

大小为$i$的排列有$a_i$种，现在计算大小为$i$的由一个或多个排列组成的排列方案。

$$
g(x)=\sum_{j=1}^{+\infty}\frac{f(x)^j}{j!}
$$

