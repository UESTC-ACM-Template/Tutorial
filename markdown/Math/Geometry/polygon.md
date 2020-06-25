给定$n$条长度分别为$a_1,\cdots,a_n$的线段，求由这些线段围成的凸多边形的最大面积。

定理(Cramer)：当面积最大时，多边形的顶点都在圆上，且这样的圆半径唯一。

证明：咕咕咕

设
$$
\displaystyle a_p=\max\{a_i\}\\
t_i(r)=\arccos\left(1-\frac{a_i^2}{2r^2}\right)\\
f^+(r)=\sum_{i=1}^nt_i(r)-2\pi\\
f^-(r)=\sum_{i=1}^nt_i(r)-2t_p(r)
$$
则圆心在多边形内外的解分别对应$f^+(r)=0$和$f^-(r)=0$。

注：由$-1 \leq 1-\frac{a_i^2}{2r^2} \leq 1$解得$r \geq a_p/2$。

因为$t_i(r)$单调减，且$\displaystyle \lim_{r \to+\infty}f^+(r)=-2\pi$，所以若$f^+(a_p/2)<0$则$f^+(r)=0$无解。



此时有
$$
f^-(a_p/2)=f^+(a_p/2)+2\pi-2t_p(a_p/2)<0\\
\lim_{r \to \infty}f^-(r)=2 \pi +\lim_{r \to \infty}f^+(r)=0
$$
因为
$$
\lim _{x \to \infty}x \arccos(1-\frac{a^2}{2x^2})=\lim_{t \to 0}\frac{\arccos(1-t^2a^2/2)}{t}=\lim_{t \to 0}\frac{\frac{a^2 t}{\sqrt{1 - (1 - a^2 t^2/2)^2}}}{1}\\=\sqrt{\lim_{t \to 0}\frac{a^4t^2}{1-(1-a^2t^2/2)^2}}=\sqrt{\lim_{t \to 0}\frac{a^4t^2}{a^2t^2-a^4t^4/4}}=\sqrt{\lim_{t \to 0}\frac{a^4}{a^2-a^4t^2/4}}=\sqrt{a^2}=a
$$
所以有
$$
\lim_{r\to\infty}rf^-(r)=\sum_{i=1}^n a_i-2a_p>0
$$
所以$f^-(r)$在无穷远处符号为正。

所以$f^-(r)$零点

