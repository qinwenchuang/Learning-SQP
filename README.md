# Learning-SQP
==还没研究明白GitHub对LaTeX的支持，ReadMe中公式的显示很乱，可以参考ReadMe.png吧。==
#### 0. 依赖开源库
> - qpOASES
> - Eigen3

#### 1. 编译和运行
> -  mkdir build
> -  cd build
> -  make
> - ./Result
> - 可以发现和知乎内容一致（2024-09-21更新部分）https://zhuanlan.zhihu.com/p/464676135

#### 2. 初步的简要说明
##### 2.1 LBFGS的计算公式
非线性问题的Hessian矩阵一般很难求解，并且耗时严重，一般会采用Hessian矩阵的近似，考虑到求解耗时和收敛性，优先选择L-BFGS。L-BFGS的计算使用到前后帧的$x$和目标函数的偏导，为了方便描述进行如下简写，
$$\begin{align} s_k = x_{k+1}-x_k \end{align}$$
$$\begin{align} y_k = \bigtriangledown_xf(x^{k+1})-\bigtriangledown_xf(x^{k}) \end{align}$$
$$\begin{align} r_k = \theta_ky_k + (1-\theta_k)B_k s_k \end{align}$$
$$\begin{align} \theta_k = \left\{\begin{matrix} 1 & if s_k^T y_k \ge 0.2s_k^TB_ks_k\\ (0.8s_k^TB_k s_k)/(s_k^TB_k s_k-s_k^Ty_k) & if s_k y_k < 0.2s_k^TB_k s_k \end{matrix}\right. \end{align}$$
$B_k$的更新公式为，$B_k$的更新公式为，
$$\begin{align} B_{k+1} = B_k - \frac{B_ks_k s_k^TB_k}{s_k^TB_k s_k} + \frac{r_kr_k^T}{s_k^Tr_k} \end{align}$$

#### 3. TODO
> - 验证积分器的偏导
> - 补全轨迹规划的示例
> - 等