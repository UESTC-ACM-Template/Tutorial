# 生成函数

## 一般生成函数(OGF)

$$f(x)=\sum_{i=0}^{+\infty}a_i x^i$$

一般生成函数一般用于组合的计数（位置不可区分，或无标号）

### OGF卷积的组合意义

考虑$f$和$g$的卷积：

$$h(x)=f(x)g(x)=\sum_{i=0}^{+\infty}\left(\sum_{j=0}^ia_jb_{i-j}\right)x^i=\sum_{i=0}^{+\infty}c_i x^i$$

可以理解为：有两类物品组合。第一类中大小为$i$的有$a_i$种，第二类中大小为$i$的有$b_i$种。则分别从第一类和第二类中分别选出两个组合将其合并，大小为$i$的有$c_i$种。

### 欧拉变换

大小为$i$的组合有$a_i$种，现在计算总大小为$i$的由一个或多个组合组成的组合方案数量。

因为组合之间不可区分，所以可以枚举大小为$i$的组合中每一个用了多少。

即

$$\mathscr E(f)(x)=\left(1+x+x^2+\cdots\right)^{a_1}\left(1+x^2+x^4+\cdots\right)^{a_2}\cdots=\prod_{i=1}^{+\infty}\left(\frac{1}{1-x^i}\right)^{a_i}
$$

$$= \exp\left(\sum_{i=0}^{+\infty} a_i \ln \frac{1}{1-x^i}\right)=\exp\left(\sum_{i=0}^{+\infty} -a_i \ln{(1-x^i)}\right)$$
$$ = \exp\left(\sum_{i=0}^{+\infty} a_i \sum_{j=1}^{+\infty}\frac{x^{ij}}{j}\right)
= \exp\left(\sum_{j=1}^{+\infty}\frac{1}{j}\sum_{i=0}^{+\infty}a_ix^{ij}\right)=\exp\left(\sum_{j=1}^{+\infty}\frac {f(x^j)}{j}\right)
$$

## 指数生成函数(EGF)

$$f(x)=\sum_{i=0}^{+\infty}a_i \frac{x^i}{i!}$$

指数生成函数一般用于排列的计数（位置可区分，或有标号）

### EGF卷积的组合意义

考虑$f$和$g$的卷积：

$$h(x)=f(x)g(x)=\sum_{i=0}^{+\infty}\left(\sum_{j=0}^i\frac{a_j}{j!}\frac{b_{i-j}}{(i-j)!}\right)x^i=\sum_{i=0}^{+\infty}\left(\sum_{j=0}^i\binom ija_jb_{i-j}\right)\frac{x^i}{i!}=\sum_{i=0}^{+\infty}c_i x^i$$

可以理解为：有两类物品排列。第一类中大小为$i$的有$a_i$种，第二类中大小为$i$的有$b_i$种。则分别从第一类和第二类中分别选出两个排列将其合并，大小为$i$的有$c_i$种。

注：先钦定来自第一类的物品放在哪些位置，来自第二类的物品放在剩下那些位置，因此有个二项式系数。

### $\exp f(x)$的组合意义

大小为$i$的排列有$a_i$种，现在计算大小为$i$的由一个或多个排列组成的排列方案。

$$g(x)=\sum_{j=1}^{+\infty}\frac{f(x)^j}{j!}$$
