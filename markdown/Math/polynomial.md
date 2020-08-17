# 多项式

## 多项式基础

### 快速卷积原理

给定$f(x),g(x)$，计算$h(x)=f(x)g(x)$。

因为一个$n-1$次多项式的系数可以被其在$n$个不同位置上的取值唯一确定（考虑对其在$n$个位置上的取值列一个$n$元线性方程组，其系数矩阵即为范德蒙德矩阵，因范德蒙德行列式不为0当且仅当$x_i$两两不同，所以该方程组有唯一解），所以如果对于一个$n-1$次多项式我们能够快速求出$f(x)$在$n$个不同位置上的取值并根据其反推原多项式的系数，卷积即可快速完成。

如已知$f$与$g$之积是某个不超过$n-1$次的多项式，则可利用某种神奇算法快速计算出$f(x_0), \cdots, f(x_{n-1})$和$g(x_0),\cdots,g(x_{n-1})$，然后$O(n)$算出所有$h(x_i)=f(x_i)g(x_i)$，之后再用某种神奇算法从$h(x_0),\cdots,h(x_{n-1})$还原处$h$的各项系数。

### 快速傅里叶变换

快速傅里叶变换可在$O(n \log n)$内从$n-1$次多项式的$n$个系数计算出多项式在$n$个不同位置的取值。

设$n$次单位根$\omega_n$满足$\omega_n^{n}=1$且$\forall i \in \{ 1, 2, ..., n - 1\}$，$\omega_n^i \neq 1$，选取的位置即为$x_i=\omega_n^{i},i\in\{0,1,...,n-1\}$。

考虑利用单位根$\omega_n$分治对系数序列进行变换，即求出$f(\omega_n^k)，k\in \{0,1,...,n-1\}$，这里设$n=2^w$。

将序列按下标的奇偶分成两份，即

$$
f_0(x)=\sum_{i=0}^{\frac n2-1} a_{2i}x^{i}
$$

$$
f_1(x)=\sum_{i=0}^{\frac n2-1}a_{2i+1}x^i
$$

对每份分别进行FFT可得$f_0(\omega_{\frac n2}^{k}),f_1(\omega_{\frac n2}^{k}),k \in \{0,1,...,\frac n2-1\}$。

注：对$f_0$和$f_1$进行FFT时所用的单位根为$\omega_{\frac n2}=\omega_n^{2}$，即分治下去所得结果为$f_0,f_1$在$\omega_n^{2k},k\in\{0,1,...,\frac n2 -1\}$上的取值。

$$
f(\omega_n^k)=\sum_{i=0}^{n-1}a_i \omega_n^{ki}=\sum_{i=0}^{\frac n2-1}a_{2i} \omega_n^{2ki}+\sum_{i=0}^{\frac n2-1}a_{2i+1}\omega_n^{2ki+k}
$$

$$
=\sum_{i=0}^{\frac n2-1}a_{2i} \omega_n^{2ki}+\omega_n^k\sum_{i=0}^{\frac n2-1}a_{2i+1}\omega_n^{2ki}=f_0(\omega_n^{2k})+\omega_n^kf_1(\omega_n^{2k})
$$

$$
f(\omega_n^{k+\frac n2})=f(-\omega_n^k)=\sum_{i=0}^{n-1}a_i(-\omega_n^{k})^i=\sum_{i=0}^{\frac n2 - 1}a_{2i}(-\omega_n^k)^{2i}+\sum_{i=0}^{\frac n2 - 1} a_{2i+1}(-\omega_n^{k})^{2i+1}
$$

$$
=\sum_{i=0}^{\frac n2 - 1}a_{2i}\omega_n^{2ki}-\omega_n^k\sum_{i=0}^{\frac n2 - 1}a_{2i+1}\omega_n^{2ki}=f_0(\omega_n^{2k})-\omega_n^{k}f_1(\omega_n^{2k})
$$
注：上述过程的$k \in \{0,1,...,\frac n2 -1\}$

于是每层分治可以在$O(n)$复杂度内合并，总复杂度$O(n \log n)$。

### 快速傅里叶逆变换

快速傅里叶逆变换可在$O(n \log n)$内从多项式在$n$个不同位置的取值推出多项式的$n$个系数。

考虑对$f$利用单位根$\omega_n$进行FFT得到$g$，再对$g$利用单位根$\omega_n^{-1}$进行FFT得到$h$，并设$g,h$的系数分别为$\{b_j\},\{c_k\}$
$$
b_j=\sum_{i=0}^{n-1}a_i\omega_n^{ij}
$$

$$
c_k=\sum_{j=0}^{n-1}b_j\omega_n^{-j}=\sum_{j=0}^{n-1}\left(\sum_{i=0}^{n-1}a_i\omega_n^{ij}\right)\omega_n^{-jk}=\sum_{i=0}^{n-1}a_i\sum_{j=0}^{n-1}\omega_n^{ij}\omega_n^{-jk}
=\sum_{i=0}^{n-1}a_i\sum_{j=0}^{n-1}\omega_n^{(i-k)j}
$$
注意到当$i=k$时有
$$
\sum_{j=0}^{n-1}\omega_n^{(i-k)j}=\sum_{j=0}^{n-1}1^j=n
$$
当$i \neq k$时有
$$
\sum_{j=0}^{n-1}\omega_n^{(i-k)j}=\frac{1-\omega_n^{(i-k)n}}{1-\omega_n^{i-k}}=\frac{1-1}{1-\omega_n^{i-k}}=0
$$
所以$c_k=na_k$

### 快速数论变换

在一般的FFT中，取$\omega_n=\cos \frac {2 \pi}{n}+i \sin \frac {2 \pi}{n}$并在复数域上进行运算即可。在系数较大时会产生精度问题，但使用`long double`可大幅度提升精度与TLE风险。

有时我们面对的问题是模素数意义下的，这时我们可以考虑直接在模素数有限域下解决快速傅里叶变换的问题。对于素数$p=2^k r+1$，模$p$域乘法群为$p-1$次循环群$Z_{p-1}$，且有子群$Z_{2^k}=\{g^{ri}|i \in \{0,1,...,2^k-1\}\}$，其生成元的阶为$2^k$，于是我们可取$\omega_n=g^r$这里$g$为模$p$意义下的原根。因此对于素数$p=2^kr+1$，我们可在模$p$意义下进行长度为$2^k$或其因子的FFT。

```cpp
namespace NTT1 {

using ::mul;
using ::inv;

const int W = 18, S = 1 << W, g = 3;
int w[S + 1], rev[S << 1], *r[W + 1];
void init() {
    for (int s = 0; s <= W&&(r[s]=rev+(1<<s),1); ++s)
        for (int i = 0; i != (1 << s); ++i)
            r[s][i] = (r[s][i >> 1] >> 1) | ((i & 1) << (s - 1));
    w[0] = 1; w[1] = qpm(g, (P - 1) / S);
    for (int i = 2; i <= S; ++i) w[i] = mul(w[i - 1], w[1]);
}

int m, s, im;
int init(int n) {
    for (s = 0, m = 1; m < n; m <<= 1, ++s);
    im = inv(m); return m;
}

void ntt(int* p, int t) {
    for (int i = 0; i != m; ++i) if (i < r[s][i]) swap(p[i], p[r[s][i]]);
    for (int i = 1, z = 0; i != m; i <<= 1, z++)
        for (int j = 0; j != m; j += (i << 1))
            for (int k = 0, u, v; k != i; k++)
                u = p[j+k], v = mul(w[(t?(i<<1)-k:k)<<W-z-1], p[i+j+k]),
                p[j + k] = add(u, v), p[i + j + k] = sub(u, v);
    if (t) for (int i = 0; i != m; ++i) p[i] = mul(p[i], im);
}

int px[S], py[S];
vi pmul(const vi& p1, const vi& p2, int n = 0) {
    int n1 = p1.size(), n2 = p2.size(), n3 = n1 + n2 - 1;
    init(n3);
    copy_n(p1.begin(), n1, px); fill(px + n1, px + m, 0);
    copy_n(p2.begin(), n2, py); fill(py + n2, py + m, 0);
    ntt(px, 0); ntt(py, 0);
    for (int i = 0; i != m; ++i) px[i] = mul(px[i], py[i]);
    ntt(px, 1); vi p3(n3); copy_n(px, n3, p3.begin());
    if (n && n3 > n) p3.resize(n); return p3;
}

vi pinv(const vi& p1) {
    int n1 = p1.size(), n2 = (n1 + 1) >> 1;
    if (n1 == 1) return { inv(p1[0]) };
    else {
        vi p2 = pinv(vi(p1.begin(), p1.begin() + n2));
        init(n1 << 1);
        copy_n(p1.begin(), n1, px); fill(px + n1, px + m, 0);
        copy_n(p2.begin(), n2, py); fill(py + n2, py + m, 0);
        ntt(px, 0); ntt(py, 0);
        for (int i = 0; i != m; ++i)
            px[i] = mul(sub(2, mul(px[i], py[i])), py[i]);
        ntt(px, 1); vi p3(n1); copy_n(px, n1, p3.begin());
        return p3;
    }
}

}

using NTT1::init;
using NTT1::pmul;
using NTT1::pinv;
```

### 任意模数NTT

9次普通NTT，巨大常数，$O(n \log n)$

碰到了建议重新思考生成函数之外的做法或弃题（认真）。

```cpp
namespace NTT2 {

    typedef array<int, 3> mint;
    const int M[3] = { 167772161, 469762049, 998244353 }, g = 3;

    template<const int i> int add(int a, int b) { int r = a + b; return r < M[i] ? r : r - M[i]; }
    template<const int i> int sub(int a, int b) { int r = a - b; return r < 0 ? r + M[i] : r; }
    template<const int i> int mul(long long a, long long b) { return a * b % M[i]; }
    template<const int i> int inv(int a) { return a == 1 ? a : mul<i>(inv<i>(M[i] % a), M[i] - M[i] / a); }
    template<const int i> int qpm(int a, int b, int r = 1) {
        do if (b & 1) r = mul<i>(r, a); while (a = mul<i>(a, a), b >>= 1);
        return r;
    }

    inline mint mul(mint a, mint b) { return { mul<0>(a[0], b[0]), mul<1>(a[1], b[1]), mul<2>(a[2], b[2]) }; }

    inline int crt(const mint& x){
        static const int o0 = inv<1>(M[0]);
        static const int o1 = inv<2>(1ll * M[0] * M[1] % M[2]);
        long long a = 1ll * (x[1] - x[0] + M[1]) * o0 % M[1] * M[0] + x[0];
        return (1ll * (x[2] - a % M[2] + M[2]) * o1 % M[2] * M[0] % P * M[1] + a) % P;
    }

    const int W = 18, S = 1 << W;
    int rev[S << 1], *r[W + 1];
    mint w[S + 1];
    void init() {
        for (int s = 0; s <= W&&(r[s]=rev+(1<<s),1); ++s)
            for (int i = 0; i != (1 << s); ++i)
                r[s][i] = (r[s][i >> 1] >> 1) | ((i & 1) << (s - 1));
        w[0] = { 1, 1, 1 };
        w[1][0] = qpm<0>(g, (M[0] - 1) / S);
        w[1][1] = qpm<1>(g, (M[1] - 1) / S);
        w[1][2] = qpm<2>(g, (M[2] - 1) / S);
        for (int i = 2; i <= S; ++i) w[i] = mul(w[i - 1], w[1]);
    }

    int m, s; mint im;
    void init(int n) {
        for (s = 0, m = 1; m < n; m <<= 1, ++s);
        im = { inv<0>(m), inv<1>(m), inv<2>(m) };
    }

    void ntt(mint* p, int t) {
        for (int i = 0; i != m; ++i) if (i < r[s][i]) swap(p[i], p[r[s][i]]);
        for (int i = 1, z = 0; i != m; i <<= 1, z++)
            for (int j = 0; j != m; j += (i << 1))
                for (int k = 0; k != i; k++) {
                    mint u=p[j+k], v=mul(w[(t?(i<<1)-k:k)<<W-z-1], p[i+j+k]);
                    p[j + k] = { add<0>(u[0], v[0]), add<1>(u[1], v[1]), add<2>(u[2], v[2]) };
                    p[i + j + k] = { sub<0>(u[0], v[0]), sub<1>(u[1], v[1]), sub<2>(u[2], v[2]) };
                }
        if (t) for (int i = 0; i != m; ++i) p[i] = mul(p[i], im);
    }

    mint px[S], py[S];
    vi pmul(const vi& p1, const vi& p2, int n = 0) {
        int n1 = p1.size(), n2 = p2.size(), n3 = n1 + n2 - 1;
        init(n3);
        for (int i = 0; i != n1; ++i) px[i] = { p1[i] % M[0], p1[i] % M[1], p1[i] % M[2] };
        for (int i = n1; i != m; ++i) px[i] = { 0, 0, 0 };
        for (int i = 0; i != n2; ++i) py[i] = { p2[i] % M[0], p2[i] % M[1], p2[i] % M[2] };
        for (int i = n2; i != m; ++i) py[i] = { 0, 0, 0 };
        ntt(px, 0); ntt(py, 0);
        for (int i = 0; i != m; ++i) px[i] = mul(px[i], py[i]);
        ntt(px, 1);
        vi p3(n3);
        for (int i = 0; i != n3; ++i) p3[i] = crt(px[i]);
        if (n && n3 > n) p3.resize(n); return p3;
    }

    vi pinv(const vi& p1) {
        int n1 = p1.size(), n2 = (n1 + 1) >> 1;
        if (n1 == 1) return { ::inv(p1[0]) };
        vi p2 = pinv(vi(p1.begin(), p1.begin() + n2));
        return pmul(psub({2}, pmul(p1, p2)), p2, n1);
    }
}

using NTT2::init;
using NTT2::pmul;
using NTT2::pinv;
```

### 循环卷积



## 多项式操作

### 卷积

给定$f(x),g(x)$，求$h(x)=f(x)g(x)$。

对$f(x)$与$g(x)$分别进行FFT/NTT，计算逐点乘积后进行逆变换即可。

在$n$小的时候暴力计算可减小常数。

```cpp
int px[S], py[S];
vi mul(const vi& p1, const vi& p2, int n = 0) {
    int n1 = p1.size(), n2 = p2.size(), n3 = n1 + n2 - 1;
    init(n3);
    copy_n(p1.begin(), n1, px); fill(px + n1, px + m, 0);
    copy_n(p2.begin(), n2, py); fill(py + n2, py + m, 0);
    ntt(px, 0); ntt(py, 0);
    for (int i = 0; i != m; ++i) px[i] = mul(px[i], py[i]);
    ntt(px, 1); vi p3(n3); copy_n(px, n3, p3.begin());
    if (n && n3 > n) p3.resize(n); return p3;
}
```

### 分治卷积

给定多项式$f(x)$，求$g(x)$满足$g_i=\sum_{j=0}^{i-1}f_{i-j}g_j$，其中$g_0$已给定。

使用分治计算$g_{l ... r-1}$。

设$m = \lceil (l+r)/2\rceil,w=m-l$。

对左半边分治下去，即计算$g_{l...m-1}$。

考虑计算左边对右边的贡献（$i \geq m$）：
$$
\sum_{j=l}^{i-1}f_{i-j}g_j=\sum_{j=l}^{m-1}f_{i-j}g_j+\sum_{j=m}^{i-1}f_{i-j}g_j
$$
考虑用卷积计算第一项。令
$$
h_j=\begin{cases}
g_{l+j} & j < w\\
0 & j \geq w
\end{cases}
$$

$$
(h*g)_{w+i}=\sum_{j=0}^{w+i}h_{j}g_{w+i-j}=\sum_{j=0}^{w-1}h_{j}g_{w+i-j}
=\sum_{j=0}^{w-1}g_{l+j}f_{w+i-j}=\sum_{j=l}^{m-1}f_{m+i-j}g_j
$$
即第一项可以通过计算$\{g_l,g_{l+1},\cdots,g_{m-1}\}$与$f$的卷积得到。

时间复杂度$O(n \log ^2 n)$

```cpp
void conv(const vi& f, vi& g, int l, int r) {
    if (l + 1 == r) return;
    else {
        int m = (l + r) >> 1;
        conv(f, g, l, m);
        vi h = mul(
            vi(g.begin() + l, g.begin() + m),
            vi(f.begin(), f.begin() + r - l + 1)
        );
        for (int i = m; i != r; ++i) f[i] = add(f[i], h[i - l]);
        conv(f, g, m, r);
    }
}
```

### 有依赖的分治卷积

给定运算$\bigoplus$，函数$t$与递推初值$g_1$和递推关系

$$f_i=\bigoplus_{j=1}^{i} g_j,g_i=t\left(\sum_{j=1}^{i-1}f_{i-j}g_j\right)$$

使用分治计算$f_{l...r-1}$和$g_{l...r-1}$。且后者中不涉及区间$l...r-1$的贡献已经计算完成。

设$m = \lceil (l+r)/2\rceil,w=m-l$。

对左半边分治下去，即计算$g_{l...m-1}$。

考虑计算左边对右边的贡献$(i \geq m)$，即计算：
$$
\sum_{(l \leq j_1 \vee l \leq j_2 )\wedge j_1 < m \wedge j_2 < m \wedge j_1+j_2=i}f_{j_1}g_{j_2}
$$
当$l=1$时，该式等价于$u=\{f_1,f_2,...,f_{m-1}\}$与$v=\{g_1,g_2,...,g_{m-1}\}$的卷积的第$i$项。

当$l \neq 1$时，根据分治的性质，当$l \neq 1$时有$j_1+j_2 \geq 2l \geq r > i$，因此$l \leq j_1$ 与 $l \leq j_2$不能同时满足。

分类讨论得该式等于
$$
\sum_{j_1=l}^{m-1}f_{j_1}g_{i-j_1}+\sum_{j_2=l}^{m-1}f_{i-j_2}g_{j_2}
$$
第一项等价于$u=\{f_l,f_{l+1}, \cdots, f_{m-1}\}$与$v=\{g_1,g_2, \cdots, g_{r-l+1}\}$的卷积的第$i-l-1$项。

第二项等价于$u=\{f_1,f_{2}, \cdots, f_{r-l+1}\}$与$v=\{g_{l},g_{l+1}, \cdots, g_{m-1}\}$的卷积的第$i-l-1$项。

时间复杂度$O(n \log ^2 n)$

```cpp
void conv(vi& f, vi& g, int l, int r) {
    if (l + 1 == r) {
        //  g[l]=t(g[l])
        //  Calculate f[l]
    }
    else {
        int m = (l + r + 1) >> 1;
        conv(f, g, l, m);
        if (l == 1) {
            vi h = mul(
                vi(f.begin() + 1, f.begin() + m),
                vi(g.begin() + 1, g.begin() + m)
            );
            for (int i = m; i != r; ++i) g[i] = add(g[i], h[i - l - 1]);
        }
        else {
            assert(m - l + 1 <= l);
            vi h1 = mul(
                vi(f.begin() + l, f.begin() + m),
                vi(g.begin() + 1, g.begin() + r - l + 1)
            ),
               h2 = mul(
                vi(f.begin() + 1, f.begin() + r - l + 1),
                vi(g.begin() + l, g.begin() + m)
            );
            for (int i = m; i != r; ++i) g[i] = add(g[i], add(h1[i - l - 1], h2[i - l - 1]));
        }
        conv(f, g, m, r);
    }
}
```

### 缩放

给定$f(x)$，求$g(x)=f(kx)$

$$g(x)=f(kx)=\sum_{i=0}^{n-1}{a_ik^ix^i}$$

```cpp
vi scale(const vi& a, ll d) {
    vi b = a;
    for (int i = 0; i != b.size(); ++i)
        b[i] = mul(b[i], qpm(d, i));
    return b;
}
```

### 平移

前置：卷积

给定$f(x)$，求$g(x)=f(x+\delta)$。
$$
f(x+\delta)=\sum_{i=0}^{n-1}{a_i\sum_{j=0}^i{\binom ij x^j \delta^{i-j}}}
=\sum_{j=0}^{n-1} x^j \sum_{i=j}^{n-1}{\binom ij a_i \delta^{i-j}}
$$

$$
=\sum_{j=0}^{n-1} x^j \sum_{i=0}^{n-j-1} \binom {n-i-1}{j}a_{n-i-1}\delta^{n-i-j-1}
=\sum_{j=0}^{n-1} x^{n-j-1} \sum_{i=0}^{j} \binom {n-i-1}{n-j-1} a_{n-i-1}\delta^{j-i}
$$
令$c_i=a_{n-i-1}(n-i-1)!$，$d_i=\frac{\delta^i}{i!}$
$$
=\sum_{j=0}^n \frac{x^{n-j}}{(n-j)!}\sum_{i=0}^j \frac{(n-i)!}{(j-i)!}a_{n-i} \delta^{j-i}
=\sum_{j=0}^n\frac{x^{n-j}}{(n-j)!}\sum_{i=0}^j c_i d_{j-i}
$$
将序列c与序列d卷起来再乘以阶乘即可。

时间复杂度$O(n \log n)$，大常数。

```cpp
vi shift(const vi& a, ll d) {
    int n = a.size();
    vi b = a, c(n);
    reverse(b.begin(), b.end());
    for (int i = 0; i != n; ++i) {
        b[i] = mul(b[i], fac[n - i - 1]);
        if (!i) c[i] = 1;
        else c[i] = mul(c[i - 1], mul(d, invs[i]));
    }
    vi r = mul(b, c, n);
    reverse(r.begin(), r.end());
    return egf(r);
}
```

### 求导

给定$f(x)$，求$g(x)=f'(x)$。

$$g(x)=\sum_{i=0}^{n-2}(i+1)a_{i+1}x^i$$

时间复杂度$O(n)$。

```cpp
vi deriv(const vi& p1) {
    int n1 = p1.size();
    vi p2(n1 - 1);
    for (int i = 1; i != n1; ++i) p2[i - 1] = mul(i, p1[i]);
    return p2;
}
```

### 积分

给定$f(x)=\sum_{i=0}^{n-1}{a_ix^i}$，求$g(x)=\int_0^x{f(t)}dt$。

$$g(x)=\sum_{i=1}^{n}{\frac{a_{i-1}}{i}x^i}$$

时间复杂度$O(n)$。

```cpp
vi integ(const vi& p1) {
    int n1 = p1.size();
    vi p2(n1 + 1, 0);
    for (int i = 0; i != n1; ++i) p2[i + 1] = mul(p1[i], invs[i + 1]);
    return p2;
}
```

### 牛顿迭代

注：以下涉及到牛顿迭代过程中的$\frac n2$均自动向上取整。

给定多项式$t(g)$，求$g(x)$使得$t(g) \equiv 0 \mod x^n$。

设已求出$h(x)$使得$t(h) \equiv 0 \mod x^\frac n2$。

因为有$t(g) \equiv 0 \mod x^n$，所以有$t(g) \equiv 0 \mod x^\frac n2$，因此$g(x) \equiv h(x) \mod x^\frac n2$。

考虑$t(g)$在$h$处的泰勒展开

$$t(g) \equiv \sum_{k=0}^\infty \frac {t^{(k)}(h)}{k!}(g-h)^k \mod x^n$$

因为$g(x)-h(x) \equiv 0 \mod x^ \frac n2$，所以第二项往后的项全部为0，我们得到
$$
0 \equiv t(g) \equiv t(h)+t'(h)(g-h) \mod x^n
$$
解得
$$
g \equiv h- \frac{t(h)}{t'(h)} \mod x^n
$$
边界：当$n=1$时根据具体情况处理。

涉及到牛顿迭代的多项式算法基本都自带巨大常数，$n \log n$跑$n=10^5$要半秒左右。

### 求逆

前置：卷积

给定$f(x)$，求$g(x)$使得$f(x)g(x) \equiv 1 \mod x^{n}$。

原始形式：

$$g_0=f_0^{-1},\sum_{j=0}^{i}f_jg_{i-j}=0$$

即

$$g_i=\begin{cases}
f_0^{-1} & i=0\\
-f_0^{-1}\sum_{j=1}^{i}f_jg_{i-j} & i \geq 1
\end{cases}$$

考虑牛顿迭代：设$t(g)=\frac 1g-f$，则迭代方程为
$$
g(x) = h(x)-\frac{\frac{1}{h(x)}-f(x)}{-\frac{1}{h^2(x)}}=2h(x)-h^2(x)f(x) \mod x^n
$$
边界：当$n=1$时直接求常数的乘法逆即可。

时间复杂度$O(n \log n)$，大常数。

```cpp
vi inv(const vi& p1) {
    int n1 = p1.size(), n2 = (n1 + 1) >> 1;
    if (n1 == 1) return { inv(p1[0]) };
    else {
        vi p2 = inv(vi(p1.begin(), p1.begin() + n2));
        init(n1 << 1);
        copy_n(p1.begin(), n1, px); fill(px + n1, px + m, 0);
        copy_n(p2.begin(), n2, py); fill(py + n2, py + m, 0);
        ntt(px, 0); ntt(py, 0);
        for (int i = 0; i != m; ++i)
            px[i] = mul(sub(2, mul(px[i], py[i])), py[i]);
        ntt(px, 1); vi p3(n1); copy_n(px, n1, p3.begin());
        return p3;
    }
}
```

### 除法与取模

前置：求逆

给定$f(x)$与$g(x)$，求$q(x)$与$r(x)$使得$f(x)=q(x)g(x)+r(x)$。

其中$q(x)$为$n-m+1$次多项式，$r(x)$为$m-1$次多项式。

原始形式：直接模拟$O(n^2)$的多项式除法/取模

```cpp
pair<vi, vi> div(const vi& p1, const vi& p2) {
    int n1 = p1.size(), n2 = p2.size();
    if (n1 < n2) return { vi(), p1 };
    vi p3(n1 - n2 + 1, 0), p4 = p1;
    for (int i = n1 - 1; i >= n2 - 1; --i) {
        p3[i - n2 + 1] = mul(p4[i], inv(p2[n2 - 1]));
        for (int j = 0; j != n2; ++j)
            p4[i - n2 + 1 + j] = sub(p4[i - n2 + 1 + j], mul(p2[j], p3[i - n2 + 1]));
    }
    while (!p4.empty() && p4.back() == 0) p4.pop_back();
    return { p3, p4 };
}
```

考虑将$f(x)$的系数前后翻转

$$f_R(x)=x^{n-1}f(\frac 1x)=\sum_{i=0}^{n-1}a_ix^{n-i-1}$$

将$\frac 1x$作为参数带入上式并两边乘上$x^{n-1}$可得

$$x^{n-1}f(\frac 1x)=x^{n-m+1}q(\frac 1x) x^{m-1}g(\frac 1x)+x^{n-m+1} x^{m-1} r(\frac 1x)$$

$$f_R(x)=q_R(x)g_R(x)+x^{n-m+1}r_R(x)$$

$$f_R(x)=q_R(x)g_R(x) \mod x^{n-m+1}$$

$$q_R(x)=f_R(x)g_R^{-1}(x) \mod x^{n-m+1}$$

对$g_R$多项式求逆后卷上$f_R(x)$并抛掉多余系数即可得到$q_R(x)$，进一步可求出$r_R(x)$。

时间复杂度$O(n \log n)$，大常数。

```cpp
pair<vi, vi> div(const vi& p1, const vi& p2) {
    int n1 = p1.size(), n2 = p2.size(), n3 = n1 - n2 + 1;
    if (n3 <= 0) return { { 0 }, p1 };
    vi p1r = p1, p2r = p2;
    reverse(p1r.begin(), p1r.end());
    reverse(p2r.begin(), p2r.end());
    p1r.resize(n3, 0); p2r.resize(n3, 0);
    vi p3 = mul(p1r, inv(p2r), n3);
    reverse(p3.begin(), p3.end());
    vi p4 = sub(p1, mul(p2, p3));
    p4.resize(n2 - 1, 0);
    return { p3, p4 };
}
```

注：`first`是$q(x)$，`second`是$r(x)$。

### 开根号

前置：卷积，求逆

给定$f(x)$，求$g(x)$使得$g(x)^2 \equiv f(x) \mod x^{n}$。

考虑牛顿迭代：设$t(g)=g^2-f$，则迭代方程为
$$
g(x)=h(x)-\frac{h^2(x)-f(x)}{2h(x)}=\frac{h^2(x)+f(x)}{2h(x)} \mod x^n
$$
边界：当$n=1$时直接求常数的模意义下二次剩余即可。

时间复杂度$O(n \log n)$，大常数。

```cpp
int msqrt(int n) {
    if (!n) return 0;
    int q = P - 1, s = 0, z = 2;
    //while (~q & 1) q >>= 1, s++;
    q >>= (s = __builtin_ctzll(q));
    if (s == 1) return qpm(n, (P + 1) / 4);
    while(qpm(z, (P - 1) / 2) == 1) ++z;
    int c = qpm(z, q), t = qpm(n, q),
       r = qpm(n, (q + 1) / 2), m = s;
    while(t % P != 1) {
        ll i = 1; while(qpm(t, 1ll << i) != 1) ++i;
        ll b = qpm(c, 1ll << (m - i - 1));
        r = mul(r, b); c = mul(b, b);
        t = mul(t, c); m = i;
    }
    return min(r, P - r);
}

vi sqrt(const vi& p1) {
    int n1 = p1.size(), n2 = (n1 + 1) >> 1;
    if (n1 == 1) return { ::msqrt(p1[0]) };
    else {
        vi p2 = sqrt(vi(p1.begin(), p1.begin() + n2));
        vi p3 = mul(p2, 2); p3.resize(n1); p3 = inv(p3);
        return mul(add(mul(p2, p2, n1), p1), p3, n1);
    }
}
```

### 对数

前置：求逆，求导，积分

给定$f(x)$，求$g(x)=\ln f(x) \mod x^n$。

原始形式：

$$\ln f(x)=\sum_{j=1}^{+\infty}(-1)^{j+1}\frac{f(x)^j}{j}$$

$$g(x)=\ln f(x)=\int_0^{x}\frac{f'(t)}{f(t)}dt$$

时间复杂度$O(n \log n)$，大常数。

```cpp
vi log(const vi& p1) {
    return integ(mul(deriv(p1), inv(p1), p1.size() - 1));
}
```

### 指数

前置：卷积，对数

给定$f(x)$，求$g(x)=\exp f(x) \mod x^n$。其中$f(0)=a_0=0$。

原始形式：
$$
\exp f(x)=\sum_{i=0}^{\infty}\frac{f(x)^i}{i!}\\
b_i=\sum_{}
$$
考虑牛顿迭代：设$t(g)=\ln g-f$，则迭代方程为
$$
g(x)=h(x)-\frac{\ln h-x}{\frac{1}{h}}=h(x)(1- \ln h+f) \mod x^n
$$
边界：当$n=1$时返回常数1。

时间复杂度$O(n \log n)$，大常数。

```cpp
vi exp(const vi& p1) {
    if (p1.size() == 1) return { 1 };
    else {
        vi p2 = exp({p1.begin(),p1.begin()+(p1.size()+1>>1)});
        p2.resize(p1.size(), 0);
        return mul(p2, add(sub({ 1 }, log(p2)), p1), p1.size());
    }
}
```

### 快速幂

前置：指数

给定$f(x)$，求$g(x)=f(x)^k$。

注意到$g(x)=f(x)^k=\exp \ln f(x)^k=\exp k \ln f(x)$，因此快速幂可在$O(n \log n)$内完成。

若$a_0 \neq 1$，则在取对数前需进行如下特判：

若$a_i, i \in \{0, 1, ..., j - 1\}$均为$0$，则可除以一个$x^j$后将$g$像后移$nj$位。

若进行上一步后仍有$a_0 \neq 1$，则将整个序列除以$a_0$，后将整个序列乘以$a_0^k$。

若$nk$较小则在NTT后直接对点值表示快速幂再逆NTT即可(CF1096G)。

时间复杂度$O(n \log n)$，巨大常数。

```cpp
vi pow(const vi& p1, int k) {
    int n1 = p1.size(), n2 = n1;
    while (n2 && !p1[n1 - n2]) n2--;
    int n3 = max(n1 - 1ll * (n1 - n2) * k, 0ll);
    if (!n2 || !n3) return vi(n1, 0);
    vi p2(p1.begin() + n1 - n2, p1.begin() + n1 - n2 + n3);
    ll c = p2[0]; p2 = mul(exp(mul(log(mul(p2, ::inv(c))), k)), qpm(c, k));
    p2.resize(n1, 0); rotate(p2.begin(), p2.begin() + n3, p2.end());
    return p2;
}
```

### 欧拉变换

前置：指数

给定$f(x)$，求$\mathscr E(f)(x)=\prod_{i=1}^{+\infty}\left(1-x^i\right)^{-a_i}$

用$\exp$展开后展开$- \ln (1-x^i)$后可得
$$
\mathscr E(f)(x)=\exp\left(\sum_{j=1}^{+\infty}\frac {f(x^j)}{j}\right)
$$
$\exp$内的多项式可在$O(n \log n)$内求出。

注：可扩展到形如($\prod_{i=1}^{+\infty}\left(1-b_ix^i\right)^{a_i}$)的多项式。

时间复杂度$O(n \log n)$。

```cpp
vi eul(const vi& p1) {
    int n = p1.size();
    vi p2 = p1;
    for (int j = 2; j != n; ++j)
        for (int i = 1; i * j < n; ++i)
            p2[i * j] = add(p2[i * j], mul(p1[i], invs[j]));
    return exp(p2);
}
```

### 多点求值

前置：除法

给定$f(x)$，和$x_i, i \in \{1, 2, ..., m\}$，求$f(x_i), i \in \{1,2,...,m\}$

构造多项式$g(x)=\prod_{i=1}^m(x-x_i)$，注意到$g(x_i)=0$，考虑多项式除法

$$f(x_i)=q(x_i)g(x_i)+r(x_i)=r(x_i)$$

对$g$取模后只需对$r$求值即可。

$g_{l,r}=\prod_{i=l}^{r}(x-x_i)$可分治求出，将中间结果保存至线段树上后再进行分治，即将对$f$求其在$x_i,i \in \{l,...,r\}$上的值分治为求$f \mod g_{l,mid-1}$在$x_i,i \in \{l,...,mid-1\}$上的值和求$f \mod g_{mid,r}$在$x_i,i \in \{mid,...,r\}$上的值。因NTT常数巨大，所以在分治到一定程度时直接转秦九韶暴力求值可降低常数。

时间复杂度$O(n \log^2 n)$，巨大常数。

```cpp
vi est[N];

int eval(const vi& p, int x) {
    int r = 0;
    for (int i = (int)p.size() - 1; i >= 0; --i)
        r = add(p[i], mul(r, x));
    return r;
}

void eval0(const vi& x, int p, int lb, int rb) {
    if (lb + 1 == rb) est[p] = { sub(0, x[lb]), 1 };
    else {
        int mid = (lb + rb) >> 1;
        eval0(x, p << 1, lb, mid);
        eval0(x, p << 1 | 1, mid, rb);
        est[p] = mul(est[p << 1], est[p << 1 | 1]);
    }
}

void eval1(const vi& r, const vi& x,
           int p, int lb, int rb) {
    vi w = div(r, est[p]).second;
    if (lb + 100 >= rb)
        for (int i = lb; i != rb; ++i)
            est[0][i] = eval(w, x[i]);
    else {
        int mid = (lb + rb) >> 1;
        eval1(w, x, p << 1, lb, mid);
        eval1(w, x, p << 1 | 1, mid, rb);
    }
}

vi eval(const vi& p, const vi& x) {
    eval0(x, 1, 0, x.size());
    est[0].resize(x.size());
    eval1(p, x, 1, 0, x.size());
    return est[0];
}
```

### 快速连续插值

前置：卷积

给定$f(0),f(1),...,f(n)$，求$f(m),f(m+1),...,f(m+n)$模意义下的值。其中$m>n$。

由插值公式可得
$$
f(m+x)=\sum_{i=0}^{n}f(i)\prod_{j \neq i}\frac{m+x-j}{i-j}=\sum_{i=0}^{n}f(i) \frac{(m+x)!/(m+x-n-1)!}{(m+x-i)(-1)^{n-i}i!(n-i)!}
$$
令
$$
u_i=\frac{f(i)}{(-1)^{n-i}i!(n-i)!}
$$
当$i>n$时$u_i=0$
$$
v_i=\frac{1}{m-n+i}
$$
可得
$$
(u * v)_{x}=\sum_{i=0}^{x}\frac{f(i)}{(-1)^{n-i}i!(n-i)!(m-n+x-i)}
$$

$$
(u * v)_{n+x}=\sum_{i=0}^{n}\frac{f(i)}{(-1)^{n-i}i!(n-i)!(m+x-i)}
$$
即
$$
f(m+x)=(u*v)_{n+x}\prod_{i=m+x-n}^{m+x}i
$$

于是可用一次(任意模数)NTT求出$f(m),f(m+1),...,f(m+n)$。

总复杂度$O(n \log n)$。

```cpp
vi lintp(const vi& y, int m) {
    int n = y.size() - 1; vi z1(y), z2(2 * n + 1, 0);
    for (int i = 0; i <= 2 * n; ++i) z2[i] = inv(m - n + i);
    for (int i = 0; i <= n; ++i) {
        int c = mul(ifac[n - i], ifac[i]);
        if ((n - i) & 1) c = sub(0, c);
        z1[i] = mul(z1[i], c);
    }
    z1 = mul(z1, z2);
    z1.erase(z1.begin(), z1.begin() + n);
    z1.resize(n + 1);
    int f = 1;
    for (int i = m - n - 1; i <= m - 1; ++i) if (i) f = mul(f, i);
    for (int i = 0; i <= n; ++i) {
        f = mul(f, m + i);
        if (m - n - 1 + i) f = mul(f, inv(m - n - 1 + i));
        z1[i] = mul(z1[i], f);
    }
    return z1;
}
```
