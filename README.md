# IsoFuns


# Isotonic regression and its generalizations

IsoFuns.jl provides basic functions for solving isotonic regression problems:

+ (standard) isotonic regression problem
+ generalized isotonic regression problem
+ nearly isotonic regression problem
+ generalized nearly isotonic regression problem

## Installation

Run Julia, enter ] to bring up Julia's package manager, and add the IsoFuns.jl package:

```
julia> ]
(v1.9) pkg> add https://github.com/yutomiyatake/IsoFuns.jl
```

## Simple examples

```
using IsoFuns

n = 10
x = (1:n)/n + randn(n)
y = iso(x)
```


## References

+ R. E. Barlow, D. J. Bartholomew, J. M. Bremner and H. D. Brunk: Statistical Inference Under Order Restrictions (1972)
+ P. Groeneboom and G. Jongbloed: Nonparametric Estimation Under Shape Constraints (2014)
+ T. Matsuda and Y. Miyatake: Generalized nearly isotonic regression (2022)
+ T. Robertson, F. T. Wright and R. L. Dykstra: Order Restricted Statistical Inference (1988)
+ R. J. Tibshirani, H. Hoefling and R. Tibshirani: Nearly-isotonic regression (2011)
+ C. van Eeden: Restricted Parameter Space Estimation Problems (2006)

## Isotonic regression

Suppose that we have ``n`` independent normal observations ``X_i \sim \mathrm{N}(\mu_i,1)`` for ``i=1,\dots,n``, where
```math
\mu_1\leq \mu_2 \leq \dots \leq \mu_n
```
is a monotone sequence of means of a normal distribution.
The maximum likelihood estimate (MLE) of ``\mu`` is the solution to the constrained optimization problem:
```math
\hat{\mu} = \argmax_{\mu_1\leq \dots \leq \mu_n} \sum_{i=1}^n\log p(x_i\mid \mu_i) = \argmin_{\mu_1\leq \dots \leq \mu_n} \sum_{i=1}^n \frac{1}{2} (x_i-\mu_i)^2.
```
This formulation is 
the **isotonic regression** of ``x_1,\dots,x_n`` with uniform weights (``w_i=1`` in the formulation below).

The weighted formulation reads
```math
\hat{\mu} = \argmin_{\mu_1\leq \dots \leq \mu_n} \sum_{i=1}^n \frac{1}{2} w_i (x_i-\mu_i)^2,
```
where ``w_i > 0`` ``(i=1,\dots,n)``.

The function ```iso``` or ```iso!``` solves these problems using the algorithm called PAVA (Pool-Adjacent-Violators algorithm).

```@docs
iso
```

```@docs
iso!
```