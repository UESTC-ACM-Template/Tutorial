
# 高精度计算

## 无符号整数

```cpp
#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

namespace NTT {

const int P = 998244353, g = 3;

inline int add(int a, int b) { int r = a + b; return r < P ? r : r - P; }
inline int sub(int a, int b) { int r = a - b; return r < 0 ? r + P : r; }
inline int mul(ll a, ll b) { ll r = a * b; return r % P; }
inline int inv(int a) { return a == 1 ? a : mul(inv(P % a), P - P / a); }
inline int qpm(int a, int b) {
    int r = 1;
    do if (b & 1) r = mul(r, a);
    while (a = mul(a, a), b >>= 1);
    return r;
}

const int W = 23, S = 1 << W;
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

int p1[S], p2[S];
vector<int> conv(const vector<int>& b1, const vector<int>& b2) {
    int n1 = b1.size(), n2 = b2.size(), n3 = n1 + n2; init(n3);
    copy_n(b1.begin(), n1, p1); fill(p1 + n1, p1 + m, 0);
    copy_n(b2.begin(), n2, p2); fill(p2 + n2, p2 + m, 0);
    ntt(p1, 0); ntt(p2, 0);
    for (int i = 0; i != m; ++i) p1[i] = mul(p1[i], p2[i]);
    ntt(p1, 1); return vector<int>(p1, p1 + n3 - 1);
}
}

const size_t rad = 10;
const char* rch = "0123456789ABCDEF";
struct biu : vector<int> {

    biu(ll x = 0) {  while (x) push_back(x % rad), x /= rad; }

    void trim() { while (!empty() && !back()) pop_back(); }

    friend int cmp(const biu& b1, const biu& b2) {
        int n1 = b1.size(), n2 = b2.size();
        if (n1 != n2) return n1 < n2 ? -1 : 1;
        for (int i = n1 - 1; i >= 0; --i)
            if (b1[i] != b2[i]) return b1[i] < b2[i] ? -1 : 1;
        return 0;
    }

    friend biu operator+(const biu& b1, const biu& b2) {
        int n1 = b1.size(), n2 = b2.size();
        biu b; b.assign(max(n1, n2) + 1, 0);
        for (int i = 0; i != max(n1, n2); ++i) {
            b[i] += (i < n1 ? b1[i] : 0) + (i < n2 ? b2[i] : 0);
            b[i + 1] += b[i] / rad, b[i] %= rad;
        }
        return b.trim(), b;
    }

    friend biu operator-(const biu& b1, const biu& b2) {
        int n2 = b2.size(), j; biu b = b1;
        for (int i = n2 - 1; i >= 0; --i) {
            if (b[i] < b2[i]) {
                for (j = i + 1; !b[j]; ++j) b[j] = rad - 1;
                b[j]--, b[i] += rad;
            }
            b[i] -= b2[i];
        }
        return b.trim(), b;
    }

    friend ll to_llong(const biu& b) {
        ll r = 0, w = 1;
        for (int i : b) r += i * w, w *= 10;
        return r;
    }

    // NTT
    friend biu operator*(const biu& b1, const biu& b2) {
        if (b1.empty() || b2.empty()) return {};
        vector<int> v = NTT::conv(b1, b2);
        int n3 = v.size(); biu b; b.assign(n3 + 1, 0);
        for (int i=0;i!=n3;++i)b[i]+=v[i],b[i + 1]+=b[i]/rad,b[i]%=rad;
        return b.trim(), b;
    }

    friend biu operator/(const biu& b1, const biu& b2) {
        int r = cmp(b1, b2);
        if (r == -1) return 0;
        else if (r == 0) return 1;
        int n1 = b1.size(), n2 = b2.size(), s;
        biu y; y.assign(1, 0); y.back()=100/(b2.back()*10+(n2>1?b2.end()[-2]:0));
        for (s = 1; s / 4 <= n1; s <<= 1) {
            biu w(b2); int t = s + 1;
            if (n2 < t) w.insert(w.begin(), t - n2, 0);
            else w.erase(w.begin(), w.begin() + n2 - t);
            w = w * y; w.erase(w.begin(), w.begin() + t - 1);
            biu z; z.assign(s + 1, 0); z.back() = 2;
            y = y * (z - w);
        }
        biu b = y * b1; b.erase(b.begin(), b.begin() + s - 1 + n2);
        biu w = b * b2, w2 = w + b2;
        while (cmp(w2, b1) <= 0) b = b + 1, w = move(w2), w2 = w + b2;
        return b;
    }
};

string to_string(const biu& b) {
    int n = b.size(); if (!n) return "0";
    string s; for (int i = n - 1; i >= 0; --i) s += rch[b[i]];
    return s;
}

biu from_string_biu(const string& s) {
    biu b; for (char ch : s) b.push_back(strchr(rch, ch) - rch);
    reverse(b.begin(), b.end()); return b;
}
```

## 有符号整数

```cpp
struct bis {
    biu b; int s;
    bis operator-() const { return { b, -s }; }
    friend bis operator+(bis s1, bis s2) {
        if (!s1.s || !s2.s) return { s1.b + s2.b, s1.s + s2.s };
        if (s1.s == s2.s) return { s1.b + s2.b, s1.s };
        else {
            switch(cmp(s1.b, s2.b)){
                case 1: return { s1.b - s2.b, s1.s };
                case 0: return { biu{}, 0 };
                case -1: return { s2.b - s1.b, s2.s };
            }
            assert(false);
            return {};
        }
    }
    friend bis operator-(bis s1, bis s2) { return s1 + (-s2); }
    friend bis operator*(bis s1, bis s2) { return { s1.b * s2.b, s1.s * s2.s }; }
    friend bis operator/(bis s1, bis s2) {
        biu r = s1.b / s2.b;
        return { r, r.empty() ? 0 : s1.s * s2.s };
    }
};

string to_string(const bis& b) {
    string s = to_string(b.b);
    if (b.s == -1) s = '-' + s;
    return s;
}

bis from_string_bis(const string& s) {
    bool neg = (s[0] == '-');
    string s_ = s.substr(neg);
    return { from_string_biu(s_), neg ? -1 : s_[0] == '0' ? 0 : 1 };
}
```

## 实数

```cpp
const int maxd = 200;
struct bf {
    bis b;
    int e;
    void trim(int tgt = maxd) {
        int t = (int)b.b.size() - tgt;
        if (t < 0) {
            b.b.resize(tgt, 0);
            rotate(b.b.begin(), b.b.end() + t, b.b.end());
            e += t;
        }
        else {
            rotate(b.b.begin(), b.b.end() - tgt, b.b.end());
            b.b.resize(tgt);
            e += t;
        }
    }
    static void align(bf& f1, bf& f2) {
        bool neg = f1.e > f2.e;
        if (neg) swap(f1, f2);
        int d = f2.e - f1.e; f2.e -= d;
        f2.b.b.resize(f2.b.b.size() + d, 0);
        rotate(f2.b.b.begin(), f2.b.b.end() - d, f2.b.b.end());
        if (neg) swap(f1, f2);
    }
    friend bf operator+(const bf& b1, const bf& b2) {
        bf f1 = b1, f2 = b2; align(f1, f2);
        bf res = { f1.b + f2.b, f1.e };
        res.trim(); return res;
    }
    friend bf operator-(const bf& b1, const bf& b2) {
        bf f1 = b1, f2 = b2; align(f1, f2);
        bf res = { f1.b - f2.b, f1.e };
        res.trim(); return res;
    }
    friend bf operator*(const bf& b1, const bf& b2) {
        bf res = { b1.b * b2.b, b1.e + b2.e };
        res.trim(); return res;
    }
    friend bf operator/(const bf& b1, const bf& b2) {
        bf f1 = b1; f1.trim(2 * maxd);
        bf res = { f1.b / b2.b, f1.e - b2.e };
        res.trim(); return res;
    }
};

string to_string(const bf& b) {
    string s = to_string(b.b); int e = b.e;
    if (s[0] == '0') return s;
    while (s.back() == '0') s.pop_back(), e++;
    bool neg = (s[0] == '-');
    if (neg) s = s.substr(1);
    int sz = s.size();
    if (e < 0) {
        if (sz + e < 0) s = string(-e-sz, '0') + s;
        s.insert(s.size() + e, ".");
        if (s[0] == '.') s = '0' + s;
    }
    else if (e > 0)
        s = s + string(e, '0');
    if (neg) s = '-' + s;
    return s;
}

bf from_string_bf(const string& s) {
    size_t pos = s.find('.', 0);
    if (pos == string::npos) {
        int e = 0;
        while (e + 1 <= s.size() && s[s.size() - e - 1] == '0') e++;
        return { from_string_bis(s.substr(0, s.size() - e)), e };
    }
    else {
        int e = -(s.size() - pos - 1);
        string t = s.substr(0, pos) + s.substr(pos + 1);
        bool neg = (t[0] == '-');
        int l = t.find_first_of('0', neg);
        int r = t.find_first_not_of('0', neg);
        if (t[neg] == '0') t.erase(l, r - l);
        return { from_string_bis(t), e };
    }
}
```