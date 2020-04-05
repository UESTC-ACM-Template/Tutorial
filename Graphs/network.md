# 网络优化问题

若无特别声明，本章中的图均为大小有限的有向图。顶点集大小$|V|=n$,边集大小$|E|=m$。

在网络优化领域内，图中的顶点(vertex)被称为节点(node)，边(edge)被称为弧(arc)。因此在某些专著中网络的点集和边集分别用$N$和$A$表示。

因为作者非常懒不愿意提供例子，这里提供一本带有大量例子的网络流专著《Network-Flows-Theory-Algorithms-and-Applications》。

## 最短路

定义（$s-t$路）：给定一张图$G=(V,E)$和$s,t$两点，一条$s-t$路是一个点和边交错形成的序列$s,(s,v_1),v_1,(v_1,v_2),\cdots,(v_{k-1},v_k),v_k,(v_k,t),t$。（废话，但还是要走形式）

若给定权函数$w$，则称该路径上所有边权之和为该路径的权值/长度。

定义（简单路径）：若路径经过的点不重复，则称其为简单路径。

命题：简单路径中至多只有$|V|-1$条边。

证明：显然。

定义（单源最短路径）：给定一张图$G=(V,E)$和权函数$w:E \to \R$和源点$s$，对于其他每个点$t\in V-\{s\}$找出最短的$s-t$路。

定义（圈）：若一条$s-t$路的起点和终点相同，则称其为圈。

命题：当图中存在长度为负数的圈且存在一条$s-t$路能够经过该圈上的点，则不存在最短的$s-t$路。

证明：略。

命题：存在一条最短$s-t$简单路径。

证明：若某条最短$s-t$路经过了重复的点，则其在第一次和第二次碰到该点之间走过了一个圈。因为最短$s-t$路存在所以图中不存在负权圈，所以将该圈删去不会增加权值。因此要么这条路径不是最短路导出矛盾，要么可以一直迭代下去将所有圈全部删掉得到一条简单路径。

以下如无特别声明则默认存在最短$s-t$路。

最短路问题有一个良好的性质。

命题：对于任意一条最短$s-t$路其中任意一条边$(u,v)$，该最短路中$s-u$部分和$v-t$部分分别是最短$s-u$路和最短$v-t$路。

利用这个性质可以设计出非常朴素的最短路算法。

算法（Bellman-Ford）：

定义DP状态`dp[i][u]`表示由$i$条边组成的最短$s-u$路长度。

边界为`dp[0][s]=0,dp[0][v]=inf`，其中`v!=s`。

对于$i$从$0$到$n-1$，对s于$u\in V$，枚举$u$的出边$(u,v)$并将其接在最短$s-v$路后面，将其转移至`dp[i+1][v]`。

最终的`dp[n-1][t]`即是最短$s-t$路长度。在转移过程中维护一个`pre`数组即可从`t`开始回溯至`s`来将具体方案构造出来（若发生转移则将`pre[v]`设成`u`）。

因为每一轮只有发生转移的`v`才会对下一轮产生新影响，所以可以利用一个队列来维护每一轮发生转移的节点。

由此获得的算法被称为SPFA。这两个算法的复杂度均为$O(nm)$。

正确性证明：考虑前面的命题，不妨设最短$s-t$路中的最后一条边为$(v,t)$，而在前一轮迭代中求出了最短$s-v$路，由归纳法可得最短$s-t$路一定会被求出。且因为简单路径中边数不超过$n-1$，所以在$n-1$轮枚举后要么求出了一条最短$s-t$简单路，要么图中存在负圈。

如未保证图中不存在负圈，则可在$n-1$轮迭代的基础上继续迭代，若仍然发生了状态更新，则说明存在负圈。在SPFA中体现为节点入队次数超过$n-1$次。

对于无负权边的图，存在另一个精妙的贪心算法。

算法（Dijkstra）：

维护一个点集$S$和到其他所有点的已知最短距离$d$。

其中$d(u)$为除$u$外不经过任何$S$以外的点的最短$s-u$路长度。

最初$S=\{s\}$。

每一轮迭代中选出属于$V-S$的$d(u)$最小的$u$将其放入$S$，并用其所有出边$(u,v)$更新$d(v):=\min\{d(v),d(u)+w(u,v)\}$。

正确性证明：

即证每一轮迭代中将$u$拿出时已经求出最短$s-u$路。

若将$u$从$V-S$中取出时$d(u)$不是最短$s-u$路长度，设最短$s-u$路最后一条边为$(v,u)$。若$v$之前已被从$V-S$中取出，则必进行过$(v,u)$这条边对$u$的更新，因此$v$此时仍属于$V-S$。所以此时$d(v) \geq d(u)$。而因为任意一条$v-u$路径长度非负，所以$v$更新得到的$d(u)$大于等于此时的$d(u)$。得证。

实现：

该过程本质是维护一个集合，支持取出权值最小元素和修改权值。

可以考虑pbds，也可以用一个`std::priority_queue`，不修改权值而是直接将新权值丢进去。

因为新权值一定比老权值小，所以碰到老权值时直接丢掉即可。



## 最大流

定义（$s-t$流）：给定一张图$G=(V,E)$和容量函数$c:E\rightarrow \R$，图$G$上的$s-t$流被定义为$f:E\rightarrow \R$，满足
$$
\forall (u,v) \in E,0 \leq f(u,v) \leq c(u,v)\\
\forall v \in V, \sum_{(u,v) \in E}f(u,v)-\sum_{(v,u) \in E}f(v,u)=
\begin{cases}
-x&v=s\\
x&v=t\\
0&else
\end{cases}\\
$$
其中$s$被称为源点，$t$被称为汇点，$x$被称为$f$的大小。

第一行描述的是每条边的流量限制，第二行描述的是每个节点的流量平衡。

第二行的第一个和式即从该点出去的流量，第二个和式即进入该点的流量。

定义（最大流）：对于一个带容量的图$G$，其最大流即大小最大的流。（废话）

寻找最大流的一个朴素思路是每次选定一条路径，将这条路径上的流量尽可能增大，然后再寻找下一条，直到找不到为止。

定义（残量网络）：给定一个图$G=(V,E)$和容量函数$c:E\rightarrow \R$，对于其上的一个$s-t$流$f$，定义残量网络为将容量函数更换成$c':E\rightarrow \R$的流量网路且满足
$$
c'(u,v)=\begin{cases}
c(u,v)-f(u,v) & (u,v) \in E\\
f(v,u) & (v,u) \in E
\end{cases}
$$

第一种情况很好理解。对于第二种情况，考虑如下情景：

存在一条流$s\xrightarrow{p_1}u\rightarrow v\xrightarrow {p_2}t$和两条路径$s\xrightarrow{p_3}v$与$u\xrightarrow{p_4}t$。

不妨设$c(s,u)=c(s,v)=c(u,t)=c(v,t)=c(u,v)=\delta$。

原网络的最大流应是$s\rightarrow u\rightarrow t$有$c$的流量，$s\rightarrow v \rightarrow t$也有$c$的流量，大小为$c+c=2c$。

但这时前面的朴素算法走进了死胡同无法寻找更大的流，原因是无法将$(u,v)$上的流撤销。

残量网络通过引入反向边来实现撤销机制。即若一条边$(u,v)$已经有$f(u,v)$的流，而通过$s \rightarrow v$和$u \rightarrow t$都可以增加$\delta \leq f(u,v)$的流，那么将$f(u,v)$减少$\delta$并将那两条路径增加$\delta$即可在保持流量平衡的前提下找到更大的流。

定义（增广路）：给定容量网络$G$和流$f$，其上的$s-t$路被称为增广路当且仅当其在残量网络中的任意一条边剩余容量均大于$0$。

命题：若当前流已经是最大流，则残量网络中不存在增广路。

证明：若残量网络中存在增广路则可将增广路上流量增加来获得更大的流，矛盾。

后面可以证明残量网络中不存在增广路是当前流是最大流的充分条件。

定义（割）：给定一个图$G=(V,E)$和权函数$w:E\rightarrow \R$，定义$s-t$割为是顶点集$V$的一个2-划分。即$V=S \cup \bar S,S \cap \bar S = \emptyset,s \in S,t \in \bar S$。该割集的权被定义为$\displaystyle \sum_{u \in S \wedge v \in \bar S}w(u,v)$，割边集即为$\{(u,v)|u \in S \wedge v \in \bar S\}$。

命题：给定一个$s-t$割，则任意一条$s-t$路上有至少一条边属于该割的割边集。

定义（最小割）：略。（我不想走形式了）

命题：给定一个图，边权即为容量，其上的任意一个$s-t$流的大小$\leq$任意一个$s-t$割的大小。

证明：对于大小为$x$的流，将$S$中所有点的流量平衡式子加起来即得从点集$S$流出的净流量为$x$，同理流入$\bar S$的净流量也为$x$。

同时这等于该流在割边集上的流量之和。流量和$\leq $容量和，得证。

定理（最大流最小割）：对于任意一张容量网络和给定的源汇点$s,t$，最大$s-t$流的流量大小等于最小$s-t$割割边集的边权和。

考虑一个最大流，因为其残量网络中不存在增广路，所以任意一条$s-t$路的边集中都有至少一条满流边。

考虑由每一条$s-t$路上的第一条满流边组成的边集，将其之前的所有点组成$S$，之后的所有点组成$\bar S$。这是一个割集。因为最大流在该割边集上均满流（容量等于流量），所以此时流量等于割边权值和。即最大流的流量等于该割边集的权值和。由前一个命题中的$\leq$可知，这个割是最小$s-t$割。

---

前置：线性规划及其对偶

定理（最大流的线性规划表述）：

对于带权图$G$，其上的最大$s-t$流等价于以下线性规划问题的最优解。
$$
\max f\\
s.t.\\
\forall v \in V, \sum_{(u,v) \in E}x_{uv}-\sum_{(v,u) \in E}x_{vu}=
\begin{cases}
-f&v=s\\
f&v=t\\
0&else
\end{cases}\\
x_{uv}+y_{uv}=c_{uv}\\
x_{uv},y_{uv},f\geq 0\\
$$
考虑建立其对偶问题。

对于点$u$的流量平衡，其对偶变量为$w_u$。

对于边$(u,v)$的流量限制，其对偶变量为$z_{uv}$。

对每条边$(u,v)$上的流量$x_{uv}$，其在点$u$的流量平衡中的系数为$-1$，在点$v$的流量平衡中系数为$1，$在边$(u,v)$的流量限制中中系数为$1$。

对于松弛变量$y_{uv}$，其在边$(u,v)$的流量限制中系数为$1$。

对于总流量$f$，其在点$s$的流量平衡中系数为$1$，在点$t$的流量平衡中系数为$-1$，其代价为$1$

定理（最小割的线性规划表述）：

对于带权图$G$，其上的最小$s-t$割等价于以下线性规划问题的最优解。
$$
\min c_{uv}z_{uv}\\
s.t.\\
\forall (u,v)\in E,w_v-w_u+z_{uv} \geq 0\\
\forall (u,v)\in E,z_{uv} \geq 0\\
w_s-w_t \geq 1
$$
因为$c_{uv}>0$，所以$z_{uv}$应该要尽可能小。

但有限制$z_{uv} \geq \max(0,w_u-w_v)$，所以$w_u-w_v$应该尽可能小。

注意到$w_s-w_t\geq1$，于是钦定$w_t=0$，$w_s=1$，这样可以使每条边的$w_u-w_v\leq 1$。

不难发现$w_u$的值决定其属于$S$还是$\bar S$，而$z_{uv}$由$w_u-w_v$，即由两个端点是否属于同一个集合所决定，即表示$(u,v)$是否属于割边集。

例：

有$n$个物品，对于第$i$个物品，将其标为$A$类可获得$a_i$的价值，将其标为$B$类可获得$b_i$的价值。

有$m$个关系，对于第$i$对关系，将$u_i$和$v_i$分在同一类可获得$c_i$的价值。

现在必须将这$n$个物品分为两类，求最大价值。

用$w_i$表示将第$i$个物品分到哪一类，$z_{uv}$表示$u,v,c$这组关系有没有被破坏（即$u$，$v$未被分到同一类）。

目标函数为$\min z_{uv}c_{uv}$。

于是有

$z_{uv}\geq |w_u-w_v|=\max\{w_u-w_v,w_v-w_u\}$。

拆成两个约束即
$$
z_{uv}+w_u-w_v\geq0  \\
z_{uv}+w_v-w_u\geq 0
$$
同时有$z_{uv}\geq 0$。

新建两个物品$s,t$。钦定$s$为$A$类，$t$为$B$类。

增加了$4n$个约束，形如
$$
z_{su}+w_s-w_u\geq0\\
z_{su}+w_u-w_s\geq0\\
z_{ut}+w_u-w_t\geq0\\
z_{ut}+w_t-w_u\geq0\\
$$
注意到$w_s=1,w_t=0$。

因为$z_{su}\geq 0,w_u \leq 1$，所以$z_{su}+1\geq w_u$属于废话。

因为$z_{ut}+w_u\geq0$，所以$z_{ut}+w_u-w_t\geq0$也属于废话。

对每个容量约束加入松弛量后

对偶回去即是一个$s$为源点，$t$为汇点，从$s$到$i$有容量为$a_i$的单向边，从$i$到$t$有容量为$b_i$的单向边，从$u_i$到$v_i$有容量为$c_i$的双向边或两条容量均为$c_i$的双向边。

最终答案即为所有$abc$之和减去被破坏掉的关系的价值（即最小割）。

例2：

有$n$个物品，对于第$i$个物品，将其标为$A$类可获得$a_i$的价值，将其标为$B$类可获得$b_i$的价值。

有$m$个关系，对于第$i$个关系，将$S_i$中物品同时分为$A$类可获得$c_i$的价值，同时分在$B$类可获得$d_i$的价值。




## 费用流