## 积性函数

定义（数论函数）：定义域在正整数上，且值域中元素能互相做加法和乘法（即在某个交换环中）的函数称为数论函数。

定义（积性函数）：给定数论函数$f$，若对任意互质的$a,b$有$f(ab)=f(a)f(b)$，则称$f$为积性函数。

命题：$f(1)=1$。证明显然。

命题：给定任意一个积性函数在所有质数的幂次$p^e$上的取值，则可以唯一确定这个积性函数。

证明：对于任意正整数$n$，考虑$n$的质因子分解$\displaystyle n=\prod p_i^{q_i}$，则$\displaystyle f(n)=\prod f(p_i^{q_i})$。

涉及到积性函数的计算中经常出现对某个函数$f$在某个数$n$的所有因子位置上的值求和，这个符号记为$\displaystyle \sum_{d|n}f(d)$。

即对所有能够整除$n$的数$d$统计$f(d)$的和。

### 除数函数

定义（除数函数）：$d(n)$为$n$的因子数量。

命题：设$n$的质因子分解为$\displaystyle n=\prod p_i^{q_i}$，则$\displaystyle d(n)=\prod(q_i+1)$。

证明：对于$n$的任意因子$x$，$x$的每个质因子的在$x$中幂次必然小于等于其在$n$中的幂次。因此对于在$n$中幂次为$q_i$的质因子，其在$n$的因子中有$q_i+1$种可能，且与其他质因子互相独立。

命题：除数函数是积性函数。

证明：由上式显然。

命题：$d(n)$的前缀和是$O(n \log n)$级别的。

证明：
$$
\sum_{i=1}^nd(n)=\sum_{i=1}^n\sum_{j | i}1=\sum_{j=1}^n\sum_{j | i}1=\sum_{j=1}^n\left\lfloor\frac nj \right\rfloor \leq n\sum_{j=1}^n\frac1j =O(n\log n)
$$


### 欧拉函数

定义（欧拉函数）：$\varphi(n)$为$1$到$n-1$中与$n$互质的数的数量。

命题：设$n$的质因子分解为$\displaystyle n=\prod p_i^{q_i}$，则$\displaystyle \varphi(n)=n\prod(1-\frac{1}{p_i})$。

证明：考虑容斥。设$n$的质因子分别是$p_1,p_2,\cdots,p_k$，则$1$到$n$能被$p_i$整除的数的数量是$n/p_i$，能被$p_ip_j$整除的数量是$n/p_ip_j$，由此可以写出
$$
\varphi(n)=n-\sum_{1 \leq i \leq k}\frac n{p_i}+\sum_{1 \leq i < j \leq k}\frac n{p_ip_j}- \cdots=n\left(1-\frac{1}{p_1}\right)\left(1-\frac{1}{p_2}\right)\cdots\left(1-\frac{1}{p_k}\right)
$$
注：将右边展开，每带上一个质因子都会乘上一个$-1$，因此每一项的符号为$(-1)^{质因子数量}$

命题：欧拉函数是积性函数。

证明：由上式显然。

命题：$\sum_{d|n}\varphi(d)=n$

证明：
$$
\sum_{d|n}\varphi(d)=\sum_{d|n}\varphi\left(\frac nd\right)=\sum_{d|n}\sum_{i=1}^{\frac nd}[\gcd\left(i,\frac nd\right)=1]\\
=\sum_{d|n}\sum_{k=1}^n[\gcd(k,n)=d]=n
$$
注：$\varphi(n/d)$等于$1,2,\cdots,n$内和$n$的$\gcd$为$d$的数的数量。

命题：$n$以内的欧拉函数值可在$O(n \log n)$的时间复杂度内计算出来。

由上一个命题可得$\varphi(n)=n-\sum_{d|n,d \neq n}\varphi(d)$。

```cpp
int phi[N];
void get_phi(int n) {
    for (int i = 1; i <= n; ++i) {
        phi[i] = i;
        for (int j = 2 * i; j <= n; j += i)
            phi[j] -= phi[i];
    }
}
```

### 默比乌斯函数

定义（默比乌斯函数）：
$$
\mu(n)=\begin{cases}(-1)^{k} & n=p_1p_2\cdots p_k\\1&n=1\\0&else\end{cases}
$$

命题：若$n$有平方因子，则$\mu(n)=0$。由定义显然。

命题：默比乌斯函数是积性函数。

证明：若$a,b$中任意一个有平方因子，则$ab$也有平方因子，因此$\mu(ab)=\mu(a)\mu(b)=0$。

若$a,b$均没有平方因子且互质，则$ab$也能写成一串不同质数的乘积，因此$\mu(n)=(-1)^{k_a+k_b}=(-1)^{k_a}(-1)^{k_b}=\mu(a)\mu(b)$。得证。

命题：$\sum_{d|n}\mu(n)=[n=1]$。

证明：$n=1$时显然成立。

当$n \neq 1$时，设$n$有$k$个不同的质因子$p_1,p_2,\cdots,p_k$，设该集合为$S=\{p_i|1 \leq i \leq k \} $则上式可写成
$$
\sum_{S \subseteq T}(-1)^{|S|}=\sum_{i=0}^{k}\binom ki(-1)^i=(1-1)^k=0
$$
注：由定义有$\mu$在有平方因子的$n$的因子处的取值为$0$。

命题：$n$以内的默比乌斯函数值可在$O(n \log n)$的时间复杂度内计算出来。

由上个命题可得$\mu(n)=[n=1]-\sum_{d|n,d \neq n}\mu(d)$

```cpp
int mu[N];
void get_mu(int n) {
    mu[1] = 1;
    for (int i = 1; i <= n; ++i)
        for (int j = 2 * i; j <= n; j += i)
            mu[j] -= mu[i];
}
```

### 狄利克雷卷积

定义（狄利克雷卷积）：设$f,g$为数论函数，则$f,g$的狄利克雷卷积$h = f * g$被定义为
$$
h(n)=\sum_{d|n}f(d)g\left( \frac n d \right)
$$
命题：若$f,g$都是积性函数，则$f,g$的狄利克雷卷积也是积性函数。

证明：
$$
h(n_1n_2)=\sum_{d|n_1n_2}f(d)g\left(\frac nd\right)=\sum_{d_1|n_1}\sum_{d_2|n_2}f(d_1d_2)g\left(\frac{n_1n_2}{d_1d_2}\right)\\
=\sum_{d_1|n_1}\sum_{d_2|n_2}f(d_1)f(d_2)g\left(\frac{n_1}{d_1}\right)g\left(\frac{n_2}{d_2}\right)=\sum_{d_1|n_1}f(d_1)g\left(\frac{n_1}{d_1}\right)\sum_{d_2|n_2}f(d_2)g\left(\frac{n_2}{d_2}\right)
=h(n_1)h(n_2)
$$
命题：狄利克雷卷积满足交换律。

证明：
$$
(f*g)(n)=\sum_{d|n}f(d)g\left( \frac n d \right)=\sum_{d|n}f\left( \frac n d \right)g(d)=(g*f)(n)
$$


命题：狄利克雷卷积满足结合律，即$(f*g)*h=f*(g*h)$。

证明：
$$
((f*g)*h)(n)=\sum_{d_1|n}\left(\sum_{e_1|d_1}f(e_1)g\left(\frac {d_1}{e_1}\right)\right)h\left(\frac n{d_1}\right)
$$
设$d_2=\frac{n}{e_1},e_2=\frac{n}{d_1}$，则有
$$
=\sum_{d_2|n}f\left(\frac{n}{d_2}\right)\left(\sum_{e_2|d_2}g\left(\frac{d_2}{e_2}\right)h(e_2)\right)=(f*(g*h))(n)
$$
命题：设$f_1,f_2,\cdots,f_k$为数论函数，则
$$
(f_1*f_2*\cdots*f_k)(n)=\sum_{d_1d_2\cdots d_k|n}f_1(d_1)f_2(d_2)\cdots f_k(d_k)
$$
证明：由定义一层层展开即可。

定义（单位函数）：$e(n)=[n=1]$。

命题：对于任意数论函数$f$有$f * e=1$。

定义（常数函数）：$1(n)=1$

定义（恒等函数）：$id (n)=n$，类似的有$id^k(n)=n^k$。

命题：$\mu * 1=e$。由定义显然。

这个命题能让展开形如$[\gcd(x,y,z,\cdots)=1]$的部分进行化简。

例1：求$m$以内与$n$互质的数的个数。
$$
f(m,n)=\sum_{i=1}^{m}[\gcd(i,n)=1]=\sum_{i=1}^me(\gcd(i,n))=\sum_{i=1}^m\sum_{e|\gcd(i,n)}\mu(e)\\=\sum_{i=1}^m\sum_{e|i \wedge e|n}\mu(e)=\sum_{e|n}\mu(e)\sum_{e|i\wedge i \leq m}1=\sum_{e|n}\mu(e)\left\lfloor\frac me \right\rfloor
$$
例2：求$m$以内与$n$的$\gcd$为$d$的数的个数
$$
g(m,n,d)=\sum_{i=1}^m[\gcd(i,n)=d]=\sum_{d|i \wedge i \leq m}[\gcd(i,n)=d]=\sum_{j=1}^{\lfloor m/d\rfloor}[\gcd(j,n/d)=1]=f(m/d,n/d);
$$
命题：$\varphi *1=id$。因为$(\varphi * 1)(n)=\sum_{d|n}\varphi(d)=n$。

命题：$\sum_{d|n}\mu(d)/d=\varphi(n)/n$

证明：
$$
\varphi(n)=(\varphi*1*\mu)(n)=(id*\mu)(n)=\sum_{d|n}\mu(d)\frac{n}{d}\\
$$
两边同除$n$即可。

定理（默比乌斯反演）：对于积性函数$f$，设$g=f*1$，则$f=g*\mu$。

将每项写出来即是
$$
g(n)=\sum_{d|n}f(d) \Leftrightarrow f(n)=\sum_{d|n}\mu(n/d)g(d)
$$
证明：$g*\mu=(f*1)*\mu=f*(1*\mu)=f$。

同时也有
$$
g(d)=\sum_{d|n}f(n) \Leftrightarrow f(d)=\sum_{d|n}\mu(n/d)g(n)
$$
证明：
$$
\sum_{d|n}\mu(n/d)g(n)=\sum_{d|n}\mu(n/d)\sum_{n|m}f(m)=\sum_{d|m}f(m)\sum_{d|n|m}\mu(n/d)\\=\sum_{d|m}f(m)\sum_{T|\frac md}\mu(t)=\sum_{d|m}f(m)[m/d=1]=f(d)
$$
例(CF1139D)：有一个空数列$\{ a \}$，每一轮向$\{ a \}$中加入一个范围在$[1,m]$内的随机整数，当$\{ a \}$中所有数的$\gcd$为$1$时停止，问停止时$\{ a \}$的长度$X$的期望值。
$$
E[X]=\sum_{i=1}^{\infty}iP[X=i]=\sum_{i=1}^{\infty}P[X\geq i]\\
$$
当$i=1$时$P[x\geq 1]=1$，否则$P[X \geq i+1]=P[X>i]$，即长度为$i$的值域在$[1,m]$内的随机数列的$\gcd$不为$1$的概率。

设$f(d)$为$d=\gcd \{a\}$的数列$\{a\}$数量，$g(d)$为$d|\gcd\{a\}$的数列$\{a\}$数量，则有
$$
g(d)=\lfloor m/d \rfloor^i =\sum_{d|n}f(n)\Rightarrow f(d)=\sum_{d|n}\mu(n/d)g(n)=\sum_{d|n}\mu(n/d)\lfloor m/n\rfloor^i
$$
于是
$$
P[X>i]=\frac {1}{m^i} \left(m^i-f(1)\right)=\frac {1}{m^i} \left(m^i-\sum_{n \geq 1}\mu(n)\lfloor m/n \rfloor ^i\right)\\=\frac {1}{m^i} \left(m^i-m^i-\sum_{n \geq 2}\mu(n)\lfloor m/n \rfloor ^i\right)=-\frac {1}{m^i}\sum_{n=2}^m\mu(n)\lfloor m/n \rfloor ^i
$$
因此
$$
\sum_{i=1}^{\infty}P[X\geq i]=1+\sum_{i=1}^{\infty}P[X>i]=1+\sum_{i=1}^{\infty}\left(-\frac {1}{m^i}\sum_{n=2}^m\mu(n)\lfloor m/n \rfloor ^i\right)\\=1-\sum_{n=2}^m\mu(n)\sum_{i=1}^{\infty}(\lfloor m/n \rfloor/m)^i=1-\sum_{n=2}^m\mu(n)\frac{\lfloor m/n \rfloor/m}{1-\lfloor m/n \rfloor/m}=1+\sum_{n=2}^m\mu(n)\frac{\lfloor m/n \rfloor}{\lfloor m/n \rfloor-m}
$$

## 筛法

### 区间筛



### 欧拉筛

前面提到的埃式筛的时间复杂度是$O(n \log \log n)$，没能做到线性复杂度的原因是有一些合数可能被筛去多次（如$36$被$2$和$3$分别筛了一次）。

欧拉筛做到了线性的时间复杂度，即$O(n)$内找到$n$以内的所有质数。

对于每个最小质因子为$p$的合数$i$，欧拉筛遍历小于等于$p$的所有质数$q$并将$qi$筛去。

因为对于最小质因子为$q$的合数$j$，$j/q$的最小质因子大于等于$q$，所以其必定会被$j/q$筛去。注意到这个分解有唯一性，所以其只会被$j/q$筛去。

每个合数只被筛去一次，因此欧拉筛时间复杂度是线性的。

欧拉筛的性质很适合用来处理一些积性函数的值。

考虑对于积性函数$f$，当筛去合数$qi$时如何计算$f(qi)$ 。

若$q \nmid i$，则$f(qi)=f(q)f(i)$。否则设$q$在$qi$中的幂次为$e$，则$f(qi)=f(q^e)f(qi/q^e)$。

因此需要快速获得积性函数在质数幂次处的取值，还需要预处理对于合数$i$，其最小质因子$p$在其中的幂次$e$和$p^e$。

```cpp
bool is_prime[N];
vector<int> primes;
int pe[N], pp[N];
int f[N];
int get_f(int p, int e, int q);	//	returns f[q]; q=pow(p,e);
void euler_sieve(int n) {
    fill_n(is_prime + 1, n, true);
   	is_prime[1] = 0;
    pe[1] = 0; pp[1] = 0;
    f[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
            pe[i] = 1;
            pp[i] = i;
            f[i] = get_f(i, 1, i);
        }
        for (int p : primes) {
            if (i * p > n) break;
            is_prime[i * p] = 0;
            if (i % p != 0) {
            	pe[i * p] = 1;
                pp[i * p] = p;
                f[i * p] = f[i] * f[p];
            }
            else {
                pe[i * p] = pe[i] + 1;
                pp[i * p] = pp[i] * p;
                f[i * p] = get_f(p, pe[i * p], pp[i * p]) * f[i / pp[i]];
                break;
            }
        }
    }
}

```

对于常见的积性函数，`get_f`的取值如下

$d(n):$`getf_(p,e,q)=(e+1)`

$\mu(n):$`get_f(p,e,q)=(e==1?-1:e==0)`

$\varphi(n):$`get_f(p,e,q)=(e==0?1:(q/p)*(p-1))`

上面的代码对于$\mu,\varphi$等简单的积性函数还能进一步简化。

### 杜教筛

有一些问题涉及到求解积性函数的前缀和，且线性的时间复杂度无法满足要求。

给定积性函数$f$，若存在积性函数$g,h$满足$f*g=h$且$g$和$h$的前缀和能够很快求出，则可用下面的式子
$$
S_h(n)=\sum_{i=1}^nh(i)
=\sum_{i=1}^n\sum_{d|i}g(d)f\left(\frac id\right)
=\sum_{d=1}^ng(d)\sum_{d|i,i \leq n}f\left(\frac id \right)\\
=\sum_{d=1}^ng(d)\sum_{j=1}^{\left \lfloor \frac nd \right \rfloor}f(j)
=\sum_{d=1}^ng(d)S_f\left(\left \lfloor \frac nd \right \rfloor\right)
=g(1)S_f(n)+\sum_{d=2}^ng(d)S_f\left(\left \lfloor \frac nd \right \rfloor\right)\\
S_f(n)=\left(S_h(n)-\sum_{d=2}^ng(d)S_f\left(\left \lfloor \frac nd \right \rfloor\right)\right)/g(1)=S_h(n)-\sum_{d=2}^ng(d)S_f\left(\left \lfloor \frac nd \right \rfloor\right)
$$
递归求出$S_f(n)$的值。利用整除分块并记忆化，则计算$f(n)$时需要进行$2\sqrt n$次求和，因此最终的复杂度是
$$
T(n)=\sum_{i=1}^\sqrt n{\sqrt i}+\sum_{i=1}^{\sqrt n}\sqrt {n/i}\leq \int_0^\sqrt n\left(\sqrt n+\sqrt {n/i}+C\right)=O(n^{3/4})
$$
注：其中$C$为某个小常数。

若利用欧拉筛提前筛出前$n^{2/3}$的值，则最终的时间复杂度为
$$
T(n)=O(n^{2/3})+\sum_{i=1}^{\sqrt[3]{n}}\sqrt{n/i}=O(n^{2/3})
$$

记忆化可以不用`unordered_map`，因为只需要存储$S_f$在$\lfloor n/x \rfloor$位置的取值，所以对于$x \geq \sqrt n$可以将$S_f(x)$放在记忆化数组的下标$\lfloor n/x\rfloor$处。

```cpp
namespace sieve {

const int N = 1000001;
int sf[N];				//	f的前缀和，用欧拉筛取得
int sg(int n);			//	计算g的前缀和
int sh(int n);			//	计算h的前缀和
void eulerian_sieve();	//	...
    
int m[N]; ll n;

int cal(ll x) {
    if (x < N) return sf[x];
    int& sum = m[n / x];
    if (sum != -1) return sum;
    sum = sh(x);
    for (ll l = 2, r; l <= x; l = r + 1) {
        r = x / (x / l);
        sum = sub(sum, mul(sub(sg(r % P), sg((l - 1) % P)), cal(x / l)));
    }
    return sum;
}

//	init之后可O(1)获得所有n/x位置的取值
void init(ll n_) {
    n = n_;
    fill_n(m, (int)sqrt(n) + 2, -1);
    cal(n);
}
    
}
```

例：设$f(n)=\varphi(n)n^2$，求$f(n)$的前缀和。$n \leq 10^9$。

解：设$g=id^2$，则
$$
(f*g)(n)=\sum_{d|n}d^2\varphi\left(\frac nd \right)\frac{n^2}{d^2}=n^2\sum_{d|n}\varphi\left(\frac nd \right)=n^3=h(n)
$$

### min25筛

给定积性函数$f$，若$f$在质数位置上的取值是一个多项式$P_f(x)$且对于任意质数$p$，$f(p^e)$可以快速求，则可在$O(n^{3/4}/\log n)$的时间复杂度内求出$\displaystyle \sum_{i =1}^nf(i)$。

算法共分两步。

第一步筛出$f$在$n$以内质数位置上的取值之和。

第二步将合数位置上的取值加回去。

#### 第一步：求$g_k(i,n)$

定义：

$Primes$为质数集合，$p_i$为第$i$个质数，$\pi(x)$为$x$以内的质数个数，$m_x$为$x$的最小质因子。特别的，$m_1=1$。

$S_f(n)$为$f$的前缀和，即$\displaystyle S_f(n)=\sum_{i =1}^nf(i)$

$S(i,n)$为埃式筛运行过程中筛掉$1$和最小质因子属于前$i$个质数的合数后剩下来的数集

$f(x)$在质数位置取值相同的多项式$\displaystyle P_f(x)=\sum_{k=0}^Ka_kx^k$

$g_k(i,n)$为$x^k$在$S(i,n)$处的取值之和，即$\displaystyle g_k(i,n)=\sum_{x \in S(i,n)}x^k$

下标从$1$开始的等幂求和$\displaystyle S_k(x)=\sum_{x=1}^nx^k$

则由定义有
$$
S(i,n)=[n]-\{x|x \notin Primes \wedge m_x \leq p_i\}=\{x|x\leq n \wedge (x \in Primes \vee m_x>p_i)\}
$$
接下来考虑递推$g(i,n)$。

初始条件为除了$1$之外所有$n$以内的正整数之和，即$g_k(0,n)=S_k(n)-1$。

因为大于$\sqrt n$的质数无法筛去$n$以内的任何合数，所以
$$
\sum_{x \in Primes \wedge x \leq n}x^k=g_k(\infty,n)=g_k(\pi(\sqrt x),n)
$$
因此递推到不大于$\sqrt n$的质数就可以终止了。

不难看出要求的即是
$$
\sum_{ x \leq n \wedge x \in Primes}f(x)=\sum_{ x \leq n \wedge x \in Primes}P_f(x)=\sum_{ x \leq n \wedge x \in Primes}\sum_{k=0}^Ka_kx^k\\=\sum_{k=0}^Ka_k\sum_{ x \leq n \wedge x \in Primes}x^k=\sum_{k=0}^Ka_k g_k(\infty,n)=\sum_{k=0}^Ka_k g_k(\pi(\sqrt x),n)
$$
从$g(i-1,n)$转移到$g(i,n)$过程中筛去的是最小质因子为$p_i$的合数（余下部分的质因子均大于等于$p_i$）。
$$
S(i,n)=S(i-1,n)-\{x|x \notin Primes \wedge m_x=p_i\}\\
g_k(i,n)=g_k(i-1,n)-\sum_{x\leq n \wedge m_x =p_i}x^k\\
=g_k(i-1,n)-p_i^k\sum_{x\leq \lfloor n / p_i \rfloor \wedge m_x \geq p_i }x^k
$$
注意到
$$
\{x|x\leq \lfloor n / p_i \rfloor \wedge m_x \geq p_i\}=S(i-1,\lfloor n / p_i \rfloor)-\{p_j|j \leq i-1\}=S(i-1,\lfloor n / p_i \rfloor)-S(i - 1,p_i)
$$
于是
$$
\sum_{x,{m_x \geq p_i \wedge x\leq \lfloor n / p_i \rfloor}}x^k=g_k(i,\lfloor n / p_i \rfloor)-g_k(i,p_i)\\
$$

代入即可得递推式
$$
g_k(i,n)=g_k(i-1,n)-p_i^{k}\left[g_k(i - 1,\lfloor n / p_i \rfloor)-g_k(i - 1,p_i)\right]
$$
因为有$\lfloor \lfloor a/b\rfloor/c\rfloor=\lfloor a/(bc)\rfloor$，所以只要求出所有$g_k(i,\lfloor n/x \rfloor)$即可。

注：因为$i \leq \pi(\sqrt n)$，所以$p_i \leq \sqrt n$，不需要特地求$g_k(i,p_i)$。

#### 第二步：求$s(i,n)$

这一步的目的是从第$\pi(\sqrt n)$个质数开始依次将$f(x)$在最小质因子为$p_i$的数的位置上的取值加回来。

定义$T(i,n)=\{x|x \leq n \wedge m_x \geq p_i\}=S(i-1,n)-\{p|p \in Primes \wedge p <p_{i-1}\}$

定义$\displaystyle s(i,n)= \sum_{x \in T(i,n)}f(x)=\sum_{x \leq n \wedge m_x\geq p_i}f(x)$

边界为$\displaystyle s(\pi(\sqrt n),n)=\sum_{k=0}^Ka_k\left(g_k(\pi(\sqrt n),n)-g_k(\pi(\sqrt n),\sqrt n\right)$

最终要求的即是$S_f(n)=1+s(0,n)$。

从$s(i,n)$转移到$s(i-1,n)$需要将最小质因子为$p_i$的数加上。因为$f$不一定像$x^k$一样是完全积性的，所以需要枚举$p_i$在这些合数中的幂次$e$。这一步需要快速求$f(p^e)$。
$$
s(i,n)=s(i+1,n)+\sum_{x\leq n\wedge m_x=p_i}f(x)\\
=s(i+1,n)+\sum_{e=1,p_i^e \leq n}\left(f(p_i^e)+\sum _{x\leq \lfloor n/p_i^e\rfloor \wedge m_x>p_i}f(xp_i^e)\right)\\
=s(i+1,n)+\sum_{e=1,p_i^e \leq n}f(p_i^e)\left(1+\sum _{x\leq \lfloor n/p_i^e\rfloor \wedge m_x>p_i}f(x)\right)
$$
注意到满足$x\leq \lfloor n/p_i^e\rfloor \wedge m_x>p_i$的数集
$$
\{x|x\leq \lfloor n / p_i^e \rfloor \wedge m_x > p_i\}=\{x|x\leq \lfloor n / p_i^e \rfloor \wedge m_x \geq p_{i+1}\}=T(i+1,\lfloor n / p_i^e \rfloor)
$$
因此
$$
s(i,n)=s(i+1,n)+\sum_{e=1,p_i^e \leq n}f(p_i^e)(1+s(i+1,\lfloor n/p_i^e\rfloor))
$$

```cpp
namespace sieve {

const int N = 1000005, K = 2, a[K] = { P - 1, 1 };

int s0(int n) { return n; }
int s1(int n) { return mul(mul(n, n + 1), i2); }
int s2(int n) { return mul(mul(n, n + 1), mul(2 * n + 1, i6)); }
int (*sk[3])(int n) = { s0, s1, s2 };

bool ip[N]; int ps[N], pc;
void eulerian_sieve(int n) {
    fill_n(ip + 1, n, 1); pc = 0; ip[1] = 0;
    for (int i = 2; i <= n; ++i) {
        if (ip[i]) ps[++pc] = i;
        for (int j = 1; j <= pc && i * ps[j] <= n; ++j) {
            ip[i * ps[j]] = 0;
            if (i % ps[j] == 0) break;
        }
    }
}

//  sq为sqrt(n)，r为小于等于sq的质数个数，即\pi(\sqrt n)
//  w[i]为第i大的n/x。w[1]=n, w[c]=1。
//  如果x>sq，则g_k(i,n)其在w中的位置为id1[x],否则为id2[n/x]
ll n, sq, w[N]; int c;
int id1[N], id2[N], g[K][N];

inline int& id(ll x) { return x <= sq ? id1[x] : id2[n / x]; }
inline int cal_f(int p, int e, ll q) { return p ^ e; }

void cal_g(ll n_) {
    n = n_; sq = sqrt(n_); c = 0;
    for (ll l = 1, r; l <= n; l = r + 1) {
        ll v = w[++c] = n / l; r = n / v; id(v) = c;
        for (int k = 0; k != K; ++k)
            g[k][c] = sub(sk[k](v % P), 1);
    }
    eulerian_sieve(2 * sq);
    while (ps[pc] > sq) pc--;
    for (int i = 1; i <= pc; ++i)
        for (int j = 1, p = ps[i]; 1ll * p * p <= w[j]; ++j)
            for (int k = 0, q = 1; k != K; ++k, q = mul(q, p))
                g[k][j] = sub(g[k][j], mul(q, sub(g[k][id(w[j] / p)], g[k][id(ps[i - 1])])));
}

int cal_s(int i, ll x) {
    int p = ps[i], res = 0;
    if (x <= 1 || p > x) return 0;
    if (1ll * p * p > x) {
        for (int k = 0; k != K; ++k)
            res = add(res, mul(a[k], sub(g[k][id(x)], g[k][id(ps[i - 1])])));
        if (p == 2) res = add(res, 2);
    }
    else {
        res = cal_s(i + 1, x);
        ll q = p;
        for (int e = 1; q <= x; e++, q *= p)
            res = add(res, mul(cal_f(p, e, q), add(1, cal_s(i + 1, x / q))));
    }
    return res;
}

int cal_sf(ll n) {
    cal_g(n);
    return add(cal_s(1, n), 1);
}
    
}
```

