# IsoFuns

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yutomiyatake.github.io/IsoFuns.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yutomiyatake.github.io/IsoFuns.jl/dev/)
[![Build Status](https://github.com/yutomiyatake/IsoFuns.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yutomiyatake/IsoFuns.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yutomiyatake/IsoFuns.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yutomiyatake/IsoFuns.jl)

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
(v1.9) pkg> add IsoFuns
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