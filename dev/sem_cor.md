# Standard error of the mean of correlated time series

$$m = \sum_i x_i$$

We look for $\operatorname{std}(m) = \sqrt{Var(m)}$. Assumme all $x$ being normally distributed with the same spread: $x_i \sim N(\mu_i, \sigma)$, i.e. $\sigma^2 = Var(x)$. Then (unbiased estimator only if mean is not estimated) $s^2 \approx Var(m)$ sums all the 
entries in the Correlation matrix:

$$
\mathbf{\Sigma} = \sigma^2 \begin{bmatrix}
1 &  \rho_1 & \rho_2 & \ldots & \rho_{n-1}\\
\rho_1 &  1 & \rho_1 & \ldots & \rho_{n-2}\\
\rho_2 &  \rho_1 & 1 & \ldots & \rho_{n-3}\\
\ldots & \ldots& \ldots& \ldots& \ldots\\
 \rho_{n-1} & \rho_{n-2} & \rho_{n-3} & \ldots & 1
\end{bmatrix}
$$

where $\rho_k$ is the autocorrelation of lag $k$:

$$
\begin{aligned}
s^2 &= \sigma^2 \left( n + 2 \sum_{k=1}^{n-1} (n-k) \rho_k\right)
\end{aligned}
$$

The effective number of observations $n_{eff}$ is the number so that the variance of the mean scales as the variance of the uncorrelated mean.

$$
\begin{aligned}
Var(m) &= {\sigma^2 \over n_{eff}} = {s^2 \over n^2}
\\ 
&= \sigma^2 \left( {1\over n} + {2 \over n^2} \sum_{k=1}^{n-1} (n-k) \rho_k\right)
\\
n_{eff} &= \frac{n}{1+{2 \over n} \sum_{k=1}^{n-1} (n-k) \rho_k}
\\
s^2 &= n Var(x) {n \over n_{eff}} = n Var_{uncor}(x) {n-1 \over n_{eff}-1}
\\
Var(m) &= {\sigma^2 \over n_{eff}} = {n-1 \over n (n_{eff} -1)} Var_{uncorr}(x)
\end{aligned}
$$

$n_{eff}$ is smaller than $n$. 

The same formula for $n_{eff}$ is given in Zieba 2011, but they provide an unbiased formula for $Var(x)$ that deviates slightly from unbiased estimate of uncorrelated variance:

$$
\begin{aligned}
Var_{uncor}(x) &= {1\over (n-1)} \sum \left( x_i - \bar{x} \right)^2
\\
\sigma^2 = Var(x) &= \frac{n_{eff}}{n (n_{eff}-1)} \sum \left( x_i - \bar{x} \right)^2 
= {(n-1) n_{eff} \over n (n_{eff}-1)} Var_{uncor}(x)
\end{aligned}
$$

Futher, Zieba argues, if autocorrelationfunction is estimated, only. Then use only the components before the first negative correlation.

### Treat Missings values in the series

The simple approach of just skipping rows of missing values does not work,
because the distance between records changes leading to different entries in
the covariance matrix. In the following example a dot in covariane matrix $\Sigma$ represents a missing value for a series [1,2,.,.,5].
$$
\mathbf{\Sigma} = \sigma^2 \begin{bmatrix}
1 &  \rho_1 & . & . & \rho_4\\
\rho_1 &  1 & . & . & \rho_3\\
. & .& .& .& .\\
. & .& .& .& .\\
 \rho_4 & \rho_3 & . & . & 1
\end{bmatrix}
$$

Hence, the covariance matrix for the vector with removed missing is
$$
\mathbf{\Sigma} = \sigma^2 \begin{bmatrix}
1 &  \rho_1  & \rho_4\\
\rho_1 &  1  & \rho_3\\
 \rho_4 & \rho_3  & 1
\end{bmatrix}
$$

The count of each correlation term has to subtract the number of missing
terms for a given lag.
$$
\begin{aligned}
s^2 &= \sigma^2 \left( n_F + 2 \sum_{k=1}^{n-1} (n-k-m_k) \rho_k\right)
\text{,}
\end{aligned}
$$
where $m_k$ is the number of missing terms for lag $k$, i.e. each record where
either the original series or the lagged series has a missing. Note, that there are both, $n$, the total number of records, and $n_F$ the number of non-missing, i.e. finite, records.

Repeating the derivation of $n_{eff}$ by the same scaling of $s^2$.

$$
\begin{aligned}
Var(m) &= {\sigma^2 \over n_{eff}} = {s^2 \over n_F^2}
\\ 
&= \sigma^2 \left( {1 \over n_F} + {2 \over n_F^2} \sum_{k=1}^{n-1} (n-k-m_k) \rho_k\right)
\\
n_{eff} &= \frac{n_F}{1+{2 \over n_F} \sum_{k=1}^{n-1} (n-k-m_k) \rho_k}
\end{aligned}
$$
