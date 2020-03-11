# 数论

## 初等数论

### 整除与因子

对于正整数$a,b$，若存在正整数$k$使得$ak=b$，则称$$整除$b$，记作$a|b$。

对于正整数$a,b$，若$a|b$，则我称a$是$b$的因子。

性质1：若$a|b$，则$\frac ba |b$。证明显然。

性质2：对于$n$的任何一个因子$d$，要么$d \leq \sqrt n$，要么$n / d \leq \sqrt n$（当$\sqrt n$是整数时两者均成立）。

通过性质2可以在$O(\sqrt n)$内找到$n$的所有因子。即枚举$\sqrt n$内的所有正整数$d$，检查是否有$d|n$，若是则将$d$和$\frac nd$加入因子集（注意判断是否有$d = n/d$）。

性质3（传递性）：若$a | b, b | c$，则$a | c$。证明显然。

定义（质数）：若正整数$n \geq 2$只有$1$和$n$两个因子，则称$n$为质数（或素数），否则称为合数。

定义（质因子）：若正整数$n$能被质数$p$整除，则称$p$为$n$的质因子。

定理（唯一分解定理）：每个大于$1$的正整数$n$都能唯一的写成一串不递减的素数的乘积。

证明略。

定理（质数定理）：$n$以内的素数数量是$O(n / \ln n)$级别的。证明略。

算法（埃式筛）：

找到$n$以内的所有质数可以在$O(n \log \log n)$的时间复杂度内解决。证明略。

遍历从$2$到$n$的所有正整数，对于每个质数将其所有倍数筛去，这样就得到了原始的埃式筛。

因为所有质数里面只有$2$不是奇数，所以只需要筛掉大于等于$3$的奇合数。

因为大于$2$的偶数不可能是质数，所以对于质数$p$只需要筛掉其奇数倍的数，即$3p,5p,7p,\cdots$。

当找到质数$p$时，所有最小质因子小于$p$的合数均已被筛去，因此接下来第一个未被筛去的$p$的倍数应是$p^2$。

```cpp
vector<int> eratosthenes(int n) {
    vector<int> primes;
    vector<char> is_prime(n + 1, true);
   	is_prime[1] = 0;
    for (int i = 3; i <= n; i += 2) {
        if (!is_prime[i]) continue;
        primes.push_back(i);
        for (int j = i * i; j <= n; j += 2 * i)
            is_prime[j] = false;
    }
    return primes;
}
```

算法（朴素质因数分解）：

对于合数$n$来说，其大于$\sqrt n$的质因子至多只有一个，所以只需要检查小于等于$\sqrt n$的所有素数是否能够整除$n$。

$\sqrt n$以内的素数可以预处理，因此运行$T$次$O(n)$级别的质因数分解的时间复杂度a是$O\left(\sqrt {n} \log \log n+T \sqrt {n / \log n}\right)$。

```cpp
vector<int> prime_factorization(int n) {
    vector<int> factors;
    vector<int> primes = eratosthenes((int)sqrt(n));
    for (int p : primes) {
        if (p * p > n) break;
        while (n % p == 0) {
            factors.push_back(p);
            n /= p;
        }
    }
    if (n != 1)
        factors.push_back(n);
    return factors;
}
```

### 欧几里得算法

命题：对于正整数$a,b$，$\gcd (a,b)=\gcd(b, a \mod b)$。

证明：设$\gcd(a,b)=d_1,\gcd(b, a \mod b)=d_2$。

由定义得存在$k_1,k_2 \in \Z$使得$a=k_1d_1,b=k_2d_1$，则$a \mod b=a-\lfloor a/b \rfloor b=(k_1-\lfloor a/b \rfloor k_2)d_1$，因此$d_1|a \mod b$。所以$d_1 | d_2$。

由定义得存在$k_1,k_2 \in \Z$使得$b=k_1d_2,a \mod b=k_2d_2$，则$a=a \mod b+\lfloor a/b \rfloor b=(k_2+k_1\lfloor a/b \rfloor )d_2$，因此$d_2 | a$。所以$d_2 | d_1$。

于是$d_1=d_2$。

算法（欧几里得）：

因为每一轮迭代中$a$与$b$的和至少减少$1$，所以算法一定能够结束。

```cpp
int euclid(int a, int b) {
    if (b != 0) return euclid(b, a % b);
    else return a;
}
```

算法（扩展欧几里得）：

对欧几里得算法进行一些修改可以求得一组系数$u,v$使得$au+bv=\gcd (a,b)$。

设欧几里得算法运行时每一轮迭代参数分别是$a_1,a_2,\cdots a_n$，其中$a_n=0$，$a_{n-1}=\gcd(a,b)$，且对于$i \geq 1$有$a_{i+2}=a_i \mod a_{i-1}$。

设$a_i=q_i a_{i+1}+a_{i+2}$，其中$q_i=\lfloor a_i / a_{i+1} \rfloor$且$u_ia_i+v_ia_{i+1}=d$，则
$$
u_ia_i+v_ia_{i+1}=u_i(q_ia_{i+1}+a_{i+2})+v_ia_{i+1}=(u_iq_i+v_i)a_{i+1}+u_ia_{i+2}=u_{i+1}a_{i+1}+v_{i+1}a_{i+2}
$$
因此$u_i=v_{i+1},v_{i}=u_{i+1}-u_{i}q_{i}$。

迭代终点为$1 \cdot a_{n-1}+0\cdot a_n=\gcd(a,b)$。

```cpp
pair<int, int> ext_euclid(int a, int b) {
    if (b != 0) {
        pair<int, int> res = ext_euclid(b, a % b);
        return { res.second, res.first - (a / b) * res.second };
    }
    else return { 1, 0 };
}
```

### 一次不定方程

定义（二元一次不定方程）：形如$ax+by=c$的方程称为二元一次不定方程，其中$a,b,c$已知。

先考虑$c=0$的情况。（不要问为什么）

命题：设$x_0$为使得方程$ax_0+by=0$有整数解$y=y_0$的最小正整数，则对于方程$ax+by=0$的任意解$x,y$均存在$k$使得$x=kx_0,y=ky_0$。

证明：若$ax_1+by_1=0$且$x_0\nmid x_1$，将$x_1$除以$x_0$得$x_1=qx_0+r(0< r < x_0)$，则有$a(x_1-qx_0)+b(y_1-qy_0)=0$。令$x_2=x_1-qx_0,y_2=y_1-qy_0$，显然$x_2=r<x_0,y_2$是比$x_0$更小的正整数解，导出矛盾。因此不存在$x_1,y_1$使得$x_0 \nmid x_1$且$ax_1+by_1=0$，得证。

注意到$b|ax_0$，因此$ax_0$是$a$与$b$的最小公倍数。即$ax_0=\text{lcm} (a, b)=\frac{ab}{\gcd(a,b)}$

命题：不妨设$b>0,d = \gcd (a, b)$，则$x_0=\frac bd,y_0=-\frac ad$。

结论：二元一次不定方程$ax+by=0$的解集是$\{(x,y)|k \in \Z ,x=\frac {kb}{d},y=-\frac {ka}{d}\}$。其中$d=\gcd (a,b)$。

注意到若能够找到任意一组$x^*, y^*$使得$a x^*+by^*=c$，则方程$ax+by=c$的解集可以写成$\{(x,y)|k \in \Z ,x=\frac {kb}{d}+x^*,y=-\frac {ka}{d}+y^*\}$。其中$d=\gcd (a,b)$。

上一节中的扩展欧几里得算法可以用来解决寻找$x^*,y^*$的问题。

运行扩欧后得到$(u,v)$使得$au+bv=d$，若$d \nmid c$则方程无解（证明思路同上）。否则有$x^*=\frac {cu}{d},y^*=\frac{cv}{d}$。

总结一下，二元一次不定方程$ax+by=c$的解集是

$$
\{(x,y)|k \in \Z, x=\frac{kb+cu}{d},y=\frac{-ka+cv}{d}\}
$$

其中$u,v,d$为使用扩展欧几里得算法解$au+bv=d$所得。

## 模运算与一次同余方程

### 模运算

定义（同余）：对于模数$m$，可以定义一个$\Z$上的等价关系$\equiv_m$，$a \equiv_m b$当且仅当存在整数$k$使得$a+km=b$。一般记作$a \equiv b \mod m$。

等价关系$\equiv_m$将$\Z$划分成了$m$个剩余类$\bar 0,\bar 1,\cdots ,\overline{m-1}$。

等价类之间可以定义和整数一样的加法和乘法，这里称为模意义下的加法和乘法。

模意义下的加法和乘法满足交换律、结合律、分配律。

模意义下的除法问题可以描述为解一次同余方程$ax \equiv b \mod m$。

### 一次同余方程

定义（一次同余方程）：对于模数$m$，给定$a,b$求$x$使得$ax \equiv b \mod m$。

由同余的意义可得，一次同余方程的解集与不定方程$ax + my = b$的解集相同，即$x=\frac{km+bu}{d}$，也可写成$$x \equiv \frac {bu}{d} \mod \frac md$$

在模$m$的意义下为以下同余类：$$\overline{\frac md+\frac {bu}d},\overline{\frac {2m}d+\frac {bu}d}，\cdots$$

一次同余方程在当$b=1,\gcd(a,m)=1$时一定只有一个解$x$，定义$x$为$a$在模$m$意义下的乘法逆（逆元），记为$x=a^{-1}$。

命题：对于质数$p$，求$1$到$n$内所有数的逆元可用以下递推式求解：
$$
a^{-1}\equiv (p \mod a)^{-1}(p - \lfloor p/a\rfloor)
$$
设$p=qa+r(0 \leq r <a)$，则
$$
\left((p-q)r^{-1}\right)a=(p-q)r^{-1}(p-r)/q\equiv(p-q)(r^{-1}p-1)q^{-1}\mod p\\
=r^{-1}q^{-1}p^2-(r^{-1}+q^{-1})p+1\equiv1 \mod p
$$

### 一次同余方程组

定义（一次同余方程组）：形如

$$
\begin{cases}
x & \equiv & b_1 \mod m_1\\
x & \equiv & b_2 \mod m_2\\
& \vdots &\\
x & \equiv & b_n \mod m_n
\end{cases}
$$

称为一次同余方程组。

命题：一次同余方程组的解是一次同余式。证明略。

定理（中国剩余定理）：当$m_i$两两互质时，方程组有解

$$
x \equiv \sum_{i=1}^{n} b_iu_iv_i \mod M
$$

其中
$$
M=\prod_{i=1}^{n}m_i,u_i=\frac{M}{m_i},u_iv_i \equiv 1 \mod m_i
$$

证明：考虑第$i$项，因为$\gcd(u_i, m_i)=1$，所以$v_i$存在且唯一。因为当$i \neq j$时有$m_j|u_i$，所以$b_iu_iv_i \equiv 0 \mod m_j$，否则有$b_iu_iv_i \equiv b_i \mod m_i$。因此第$i$项对第$i$个方程的贡献是$b_i$，对其他方程都是$0$。

我们从另一个角度来考虑，即将两个一次同余方程合并。

考虑不定方程

$$
b_1+k_1m_1=b_2+k_2m_2
$$

移项

$$
k_1m_1-k_2m_2=b_2-b_1
$$

对$m_1u+(-m_2)v=d$运行扩展欧几里得算法解得$u,v,d$

则解为

$$
k_1=\frac{(b_2-b_1)u-km_2}{d},k_2=\frac{(b_2-b_1)v-km_1}{d}
$$

代入原式

$$
b_1+\frac{(b_2-b_1)um_1}{d}-k\frac{m_1m_2}{d}=b_2+\frac{(b_2-b_1)vm_2}{d}-\frac{km_1m_2}{d}
$$

解得

$$
x\equiv b_1+\frac{(b_2-b_1)um_1}{d} \mod \frac{m_1m_2}{d}
$$

将同余式依次合并至$x \equiv 0 \mod 1$即可解出模数任意的一次同余方程组。

## 积性函数与筛法

定义（数论函数）：定义域在正整数上，且值域中元素能互相做加法和乘法（即在某个交换环中）的函数称为数论函数。

定义（积性函数）：给定数论函数$f$，若对任意互质的$a,b$有$f(ab)=f(a)f(b)$，则称$f$为积性函数。

命题：$f(1)=1$。证明显然。

命题：给定任意一个积性函数在所有质数的幂次$p^e$上的取值，则可以唯一确定这个积性函数。

证明：对于任意正整数$n$，考虑$n$的质因子分解$n=\prod p_i^{q_i}$，则$f(n)=\prod f(p_i^{q_i})$。

### 除数函数

定义（除数函数）：$d(n)$为$n$的因子数量。

命题：设$n$的质因子分解为$n=\prod p_i^{q_i}$，则$d(n)=\prod(q_i+1)$。

证明：对于$n$的任意因子$x$，$x$的每个质因子的在$x$中幂次必然小于等于其在$n$中的幂次。因此对于在$n$中幂次为$q_i$的质因子，其在$n$的因子中有$q_i+1$种可能，且与其他质因子互相独立。

命题：除数函数是积性函数。

证明：由上式显然。

### 欧拉函数

定义（欧拉函数）：$\varphi(n)$为$1$到$n-1$中与$n$互质的数的数量。

命题：设$n$的质因子分解为$n=\prod p_i^{q_i}$，则$\varphi(n)=n\prod(1-\frac{1}{p_i})$。

证明：考虑容斥。设$n$的质因子分别是$p_1,p_2,\cdots,p_k$，则
$$
\varphi(n)=n-\sum_{1 \leq i \leq k}\frac n{p_i}+\sum_{1 \leq i < j \leq k}\frac n{p_ip_j}- \cdots=n\left(1-\frac{1}{p_1}\right)\left(1-\frac{1}{p_2}\right)\cdots\left(1-\frac{1}{p_k}\right)
$$
注：每一项的符号为$(-1)^{质因子数量}$

命题：欧拉函数是积性函数。

证明：由上式显然。

### 默比乌斯函数

定义（默比乌斯函数）：
$$
\sum_{d|n}\mu(d)=[n =1]
$$
命题：$\mu(1)=1$。证明显然。

命题：默比乌斯函数是积性函数。

证明：考虑第二类归纳法：
$$
\mu(n_1n_2)=\mu(n_1n_2)+\sum_{d_1|n_1}\mu(d_1)\sum_{d_2|n_2}\mu(d_2)=\mu(n_1n_2)+\sum_{d|n_1n_2,d \neq n_1n_2}\mu(d_1d_2)+\mu(n_1)\mu(n_2)\\
=\sum_{d|n_1n_2}\mu(d)+\mu(n_1)\mu(n_2)=\mu(n_1)\mu(n_2)
$$
命题：若$n$为质数，则$\mu(n)=-1$。证明显然。

命题：若$n$为$k$个质数的积，则$\mu(n)=(-1)^{k}$。证明显然。

命题：若$n=p^e(e \geq 2)$，则$\mu(n)=0$。

证明：考虑第二类归纳法：

当$e=2$时有$\mu(p^2)=-\mu(1)-\mu(p)=0$
$$
\mu(p^e)=-\sum_{i=0}^{e-1}\mu(p^i)=-\mu(1)-\mu(p)=0
$$
命题：若存在$p$使得$p^2|n$，则$\mu(n)=0$。证明显然。

由此得到默比乌斯函数的另一个等价定义：
$$
\mu(n)=\begin{cases}
(-1)^{k} & n=p_1p_2\cdots p_k\\
0&else
\end{cases}
$$

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

定义（单位函数）：$e(n)=[n==1]$。

命题：对于任意数论函数$f$有$f * e=1$。

定义（常数函数）：$1(n)=1$

定义（恒等函数）：$id (n)=n$，类似的有$id^k(n)=n^k$。

命题：$\mu * 1=e$。由定义显然。

命题：$\varphi *1=id$。证明：
$$
(\varphi * 1)(n)=\sum_{d|n}\varphi(d)=\sum_{d|n}\varphi\left(\frac nd\right)=\sum_{d|n}\sum_{i=1}^{\frac nd}[\gcd\left(i,\frac nd\right)=1]\\
=\sum_{d|n}\sum_{k=1}^{\frac nd}[\gcd\left(kd,n\right)=d]=\sum_{d|n}\sum_{i=1}^n[\gcd(i,n)=d]=n
$$
命题：
$$
\sum_{d|n}\frac {\mu(d)}{d}=\frac{\varphi(n)}{n}
$$
证明：
$$
\varphi(n)=(\varphi*1*\mu)(n)=(id*\mu)(n)=\sum_{d|n}\mu(d)\frac{n}{d}\\
$$
两边同除$n$即可。

定理（默比乌斯反演）：对于积性函数$f$，设$g=f*1$，则$f=g*\mu$。

证明：$g*\mu=(f*1)*\mu=f*(1*\mu)=f$。

## 原根与离散对数

## $\Z/n\Z$的结构

## 素性测试与因子分解

### 
