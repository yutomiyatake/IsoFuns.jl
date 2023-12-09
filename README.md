# IsoFuns


# Isotonic regression and its generalizations

[IsoFuns.jl](https://yutomiyatake.github.io/IsoFuns.jl/build/index.html) provides basic functions for solving isotonic regression problems:

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

