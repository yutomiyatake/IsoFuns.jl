var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Standard-isotonic-regression","page":"Examples","title":"Standard isotonic regression","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> using IsoFuns\n\njulia> n = 5\n5\n\njulia> x = (1:n)/n + randn(n)\n5-element Vector{Float64}:\n 0.7382449115188521\n 0.0017480727195733348\n 1.8535356206530813\n 0.440015631618591\n 1.1163182317933293\n\njulia> y = iso(x)\n5-element Vector{Float64}:\n 0.3699964921192127\n 0.3699964921192127\n 1.1366231613550006\n 1.1366231613550006\n 1.1366231613550006","category":"page"},{"location":"examples/#Generalized-nearly-isotonic-regression-(binomial)","page":"Examples","title":"Generalized nearly isotonic regression (binomial)","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Plots\nusing Distributions","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"N = 10\nn = 100\np = zeros(n)\np[1:50] = range(0.2,0.8,length=50)\np[51:100] = range(0.2,0.8,length=50)\ntrial = N*ones(n)\nsuccess = zeros(n)\nfor i=1:n\n    success[i] = rand(Binomial(N,p[i]))\nend\nx = success./trial","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"y = neariso_Binomial(success,1.0,trial)[1];","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"plot(x,st=:scatter,label=\"data\")\nplot!(y,label=\"λ=1.0\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"aic_λ,aic_value=IsoFuns.neariso_AIC_Binomial(success,trial)\nplot(aic_λ,aic_value,xlabel=\"λ\",ylabel=\"AIC\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"λ_opt = aic_λ[argmin(aic_value)]\nz = neariso_Binomial(success,λ_opt,trial)[1]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"plot(x,st=:scatter,label=\"data\")\nplot!(y,label=\"λ=1.0\")\nplot!(z,label=\"λ=opt\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"func/iso/#Isotonic-regression","page":"Isotonic regression","title":"Isotonic regression","text":"","category":"section"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Suppose that we have n independent normal observations X_i sim mathrmN(mu_i1) for i=1dotsn, where","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"mu_1leq mu_2 leq dots leq mu_n","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"is a monotone sequence of means of a normal distribution. The maximum likelihood estimate (MLE) of mu is the solution to the constrained optimization problem:","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hatmu = argmax_mu_1leq dots leq mu_n sum_i=1^nlog p(x_imid mu_i) = argmin_mu_1leq dots leq mu_n sum_i=1^n frac12 (x_i-mu_i)^2","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"This formulation is  the isotonic regression of x_1dotsx_n with uniform weights (w_i=1 in the formulation below).","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"The weighted formulation reads","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hatmu = argmin_mu_1leq dots leq mu_n sum_i=1^n frac12 w_i (x_i-mu_i)^2","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"where w_i  0 (i=1dotsn).","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"The function iso or iso! solves these problems using the algorithm called PAVA (Pool-Adjacent-Violators algorithm).","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso","category":"page"},{"location":"func/iso/#IsoFuns.iso","page":"Isotonic regression","title":"IsoFuns.iso","text":"iso(x::Vector) -> y\n\nSame as iso! with w=Nothing, but allocates an output vector y.\n\n\n\n\n\niso(x::Vector, w::Vector) -> y\n\nSame as iso!, but allocates an output vector y.\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso!","category":"page"},{"location":"func/iso/#IsoFuns.iso!","page":"Isotonic regression","title":"IsoFuns.iso!","text":"iso!(x::Vector, w::Union{Vector, Nothing}=nothing) -> x\n\nPerform isotonic regression \n\nArguments\n\nx: input vector \nw: weights, default to ones if not provided\n\nOutputs\n\nx: output vector, which is monotone\n\nAlgorithm\n\nPAVA (Pool-Adjacent-Violators algorithm)\n\n\n\n\n\n","category":"function"},{"location":"func/iso/#Generalized-isotonic-regression","page":"Isotonic regression","title":"Generalized isotonic regression","text":"","category":"section"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"The above formulation can be generalized to accommodate distributions other than the normal distribution with a uniform variance. Such a formulation is often called the generalized isotonic regression.","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"When X_isim p_i(x_imid theta_i) for i=1dotsn, where p_i(x_imidtheta_i ) is a one parameter exponential family defined by","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"p_i(x_imid theta_i) = h_i(x_i) exp (theta_i x_i -w_i psi(theta_i))","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"the problem ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hattheta = argmax_theta_1leq dots leq theta_n sum_i=1^nlog p_i(x_imid theta_i)","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"can also be solved by the function iso or iso!, as the formulation is equivalent to","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hateta = argmin_eta_1leq dots leq eta_n sum_i=1^n fracw_i2 Big(fracx_iw_i-eta_iBig)^2 quad textwith quad eta = mathrmE_theta x = psi ^prime (theta)","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Below, we list four examples of one parameter exponential families. For each case, the relation between the natural parameter theta_i and the expectation parameter eta_i = mathrmE_theta x_i = psi ^prime (theta_i) is specified. ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"","category":"page"},{"location":"func/iso/#Normal","page":"Isotonic regression","title":"Normal","text":"","category":"section"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"The probability distribution function (PDF) of the Normal Distribution with mean mu_i and standard deviation sigma_igeq 0 is given by","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"p_i(x_imid mu_i) = frac1sqrt2pisigma_i^2 exp bigg( -frac(x_i-mu_i)^22sigma_i^2 bigg)","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Since","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"beginaligned\nlog p_i(x_imid mu_i) = mu_i fracx_isigma_i^2 - fracmu_i^22sigma_i^2 - fracx_i^22sigma_i^2 - frac12log(2pisigma_i^2) \n= theta_i tildex_i - sigma_i^-2fractheta_i^22 - \nfracsigma_i^2 tildex_i^22 - frac12log(2pisigma_i^2) quad bigg(theta_i = fracmu_isigma_i^2quad  tildex_i=fracx_isigma_i^2bigg)\nendaligned","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"we can infer the following:","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"natural parameter: displaystyle theta = fracmusigma^2\ndisplaystyle psi(theta) = fracsigma^22theta^2\nweight: displaystyle w_i = sigma_i^-2\nexpectation parameter: eta_i = psi ^prime (theta_i) = sigma^2 theta = mu_i","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Suppose we have n observations X_i sim mathrmN(mu_isigma_i^2). The functions  iso_Normal and iso_Normal! solve the corresponding problem. Their input and output are as follows.","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Input ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"x = (x_1dotsx_n)\nmathrmvariance = (sigma_1^2dotssigma_n^2)","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Output","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hateta = (hateta_1dotshateta_n) = (hatmu_1dotshatmu_n), which is monotone","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Normal","category":"page"},{"location":"func/iso/#IsoFuns.iso_Normal","page":"Isotonic regression","title":"IsoFuns.iso_Normal","text":"iso_Normal(x::Vector, variance::Vector) -> y\n\nSame as iso_Normal!, but allocates an output vector y.\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Normal!","category":"page"},{"location":"func/iso/#IsoFuns.iso_Normal!","page":"Isotonic regression","title":"IsoFuns.iso_Normal!","text":"iso_Normal!(x::Vector, variance::Vector) -> x\n\nPerform isotonic regression (Normal)\n\nArguments\n\nx: input vector \nvariance a vector consisting of variances, default to ones if not provided\n\nOutputs\n\nx: output vector, which is monotone\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"","category":"page"},{"location":"func/iso/#Binomial","page":"Isotonic regression","title":"Binomial","text":"","category":"section"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"The probability mass function (PMF) of the Binomial Distribution with the number of trials N_i and the probability of success r_i in an individual trial is given by","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"p_i(x_imid r_i) = beginpmatrix N_i  x_i endpmatrix r_i^x_i (1-r_i)^N_i-x_i","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Since","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"beginaligned\nlog p_i(x_imid r_i) = big( log r_i - log(1-r_i) big) x_i + N_i log (1-r_i) + log beginpmatrix N_i  x_i endpmatrix \n= theta_i x_i - N_i log(1+mathrme^theta_i) + log beginpmatrix N_i  x_i endpmatrix\nendaligned","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"we can infer the following:","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"natural parameter: displaystyle theta_i = log r_i - log(1-r_i) quad bigg( r_i = fracmathrme^theta_i1+mathrme^theta_i bigg)\ndisplaystyle psi(theta_i) =  log(1+mathrme^theta_i)\nweight: displaystyle w_i = N_i\nexpectation parameter: displaystyle eta_i  = psi ^prime (theta) =  fracmathrme^theta_i1+mathrme^theta_i = r_i","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Suppose we have n observations X_i sim mathrmBi(N_ir_i). The functions  iso_Binomial and iso_Binomial! solve the corresponding problem. Their input and output are as follows.","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Input ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"x = (x_1dotsx_n)\nw = (N_1dotsN_2)","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Output","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hateta = (hateta_1dotshateta_n) = (hatr_1dotshatr_n), which is monotone","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Binomial","category":"page"},{"location":"func/iso/#IsoFuns.iso_Binomial","page":"Isotonic regression","title":"IsoFuns.iso_Binomial","text":"iso_Binomial(success::Vector, trial::Vector) -> x\n\nSame as iso_Binomial!, but allocates an output vector x.\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Binomial!","category":"page"},{"location":"func/iso/#IsoFuns.iso_Binomial!","page":"Isotonic regression","title":"IsoFuns.iso_Binomial!","text":"iso_Binomial!(success::Vector, trial::Vector) -> success\n\nPerform isotonic regression (Binomial)\n\nArguments\n\nsuccess: input vector consisting of the number of success\ntrials: a vector consisting of the number of trials\n\nOutputs\n\nsuccess: output vector, which is monotone\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"","category":"page"},{"location":"func/iso/#Poisson","page":"Isotonic regression","title":"Poisson","text":"","category":"section"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"The PMF of the Poisson Distribution with the average rate of occurrence lambda is gicen by","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"p(x_imid lambda_i) = fraclambda_i^x_i mathrme^-lambda_ix_i","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Since","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"beginaligned\nlog p(x_imid lambda_i) = x_i log lambda_i-lambda_i-log x_i \n= theta_i x_i-mathrme^theta_i-log x_i\nendaligned","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"we can infer the following:","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"natural parameter: displaystyle theta_i = loglambda_i\ndisplaystyle psi(theta_i) = mathrme^theta_i\nweight: displaystyle w_i = 1\nexpectation parameter: displaystyle eta_i  = psi ^prime (theta_i) = mathrme^theta_i = lambda_i","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Suppose we have n observations X_i sim mathrmPo(lambda_i). The functions  iso_Poisson and iso_Poisson! solve the corresponding problem. Their input and output are as follows.","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Input ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"x = (x_1dotsx_n)","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Output","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hateta = (hateta_1dotshateta_n) = (hatlambda_1dotshatlambda_n), which is monotone","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Poisson","category":"page"},{"location":"func/iso/#IsoFuns.iso_Poisson","page":"Isotonic regression","title":"IsoFuns.iso_Poisson","text":"iso_Poisson(x::Vector) -> y\n\nSame as iso_Poisson!, but allocates an output vector y.\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Poisson!","category":"page"},{"location":"func/iso/#IsoFuns.iso_Poisson!","page":"Isotonic regression","title":"IsoFuns.iso_Poisson!","text":"iso_Poisson!(x::Vector) -> x\n\nPerform isotonic regression (Poisson)\n\nArguments\n\nx: input vector consisting of the number of events\n\nOutputs\n\nx: output vector, which is monotone\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"","category":"page"},{"location":"func/iso/#Chisq","page":"Isotonic regression","title":"Chisq","text":"","category":"section"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"The PDF of the Chi Squared Distribution (typically written chi^2) with d_i degrees of freedom and scale parameter s_i is given by","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"p_i(x_imid s_i) = frac1Gamma(d_i2) (2s_i)^d_i2 x^d_i2-1 mathrme^-fracx_i2s_i ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Since ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"beginaligned\nlog p_i(x_imid s_i) = -frac12s_i x_i - fracd_i2 log (2s_i) + bigg( fracd_i2-1bigg)x_i - log Gamma (d_i2)  \n= theta_i x_i - fracd_i2 (-log (-theta_i)) + bigg( fracd_i2-1bigg)x_i - log Gamma (d_i2) \nendaligned","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"we can infer the following:","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"natural parameter: displaystyle theta_i = -frac12s_i\ndisplaystyle psi(theta_i) = - log (theta_i)\nweight: displaystyle w_i = fracd_i2\nexpectation parameter: displaystyle eta_i  = psi ^prime (theta_i) = -frac1theta_i = 2s_i","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Suppose we have n observations X_i sim s_i chi^2(d_i). The functions  iso_Chisq and iso_Chisq! solve the corresponding problem. Their input and output are as follows.","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Input ","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"x = (x_1dotsx_n)\nw = (d_1dotsd_2)","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"Output","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"hateta = (hateta_1dotshateta_n) = (2hats_1dots2hats_n), which is monotone","category":"page"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Chisq","category":"page"},{"location":"func/iso/#IsoFuns.iso_Chisq","page":"Isotonic regression","title":"IsoFuns.iso_Chisq","text":"iso_Chisq(x::Vector, d::Vector) -> y\n\nSame as iso_Chisq!, but allocates an output vector y.\n\n\n\n\n\n","category":"function"},{"location":"func/iso/","page":"Isotonic regression","title":"Isotonic regression","text":"iso_Chisq!","category":"page"},{"location":"func/iso/#IsoFuns.iso_Chisq!","page":"Isotonic regression","title":"IsoFuns.iso_Chisq!","text":"iso_Chisq!(x::Vector, d::Vector) -> x\n\nPerform isotonic regression (Chisq)\n\nArguments\n\nx: input vector \nd: a vector consisting of the degrees of freedom\n\nOutputs\n\nx: output vector, which is monotone\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = IsoFuns","category":"page"},{"location":"#Isotonic-regression-and-its-generalizations","page":"Home","title":"Isotonic regression and its generalizations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"IsoFuns.jl provides basic functions for solving isotonic regression problems:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(standard) isotonic regression problem\ngeneralized isotonic regression problem\nnearly isotonic regression problem\ngeneralized nearly isotonic regression problem","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Run Julia, enter ] to bring up Julia's package manager, and add the IsoFuns.jl package:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]\n(v1.9) pkg> add https://github.com/yutomiyatake/IsoFuns.jl","category":"page"},{"location":"#Simple-examples","page":"Home","title":"Simple examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using IsoFuns\n\nn = 10\nx = (1:n)/n + randn(n)\ny = iso(x)","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"R. E. Barlow, D. J. Bartholomew, J. M. Bremner and H. D. Brunk: Statistical Inference Under Order Restrictions (1972)\nP. Groeneboom and G. Jongbloed: Nonparametric Estimation Under Shape Constraints (2014)\nT. Matsuda and Y. Miyatake: Generalized nearly isotonic regression (2022)\nT. Robertson, F. T. Wright and R. L. Dykstra: Order Restricted Statistical Inference (1988)\nR. J. Tibshirani, H. Hoefling and R. Tibshirani: Nearly-isotonic regression (2011)\nC. van Eeden: Restricted Parameter Space Estimation Problems (2006)","category":"page"},{"location":"func/neariso/#Nearly-isotonic-regression","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"","category":"section"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"Suppose that we have n independent normal observations X_i sim mathrmN(mu_i1) for i=1dotsn, where mu_1mu_2dotsmu_n is a piecewize monotone sequence of normal means. The nearly isotonic regression problem is formulated as","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"beginaligned\nhatmu =argmin_mu_1 dots mu_nleft(sum_i=1^n -log p(x_imid mu_i) +lambdasum_i=1^n-1 (mu_i-mu_i+1)_+right)\n= argmin_mu_1dots  mu_nleft(sum_i=1^nfrac12 (x_i-mu_i)^2+lambdasum_i=1^n-1 (mu_i-mu_i+1)_+right)\nendaligned","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"where (a)_+=mathrmmax(a0) and lambda0 is a regularization parameter.","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"The weighted formulation reads","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"hatmu \n=argmin_mu_1 dots mu_n left( sum_i=1^n frac12 w_i (x_i-mu_i)^2+lambdasum_i=1^n-1 (mu_i-mu_i+1)_+right)","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"where w_i  0 (i=1dotsn).","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"The function neariso solves these problems. Note that the result hatmu = (hatmu_1dots hatmu_n) consists of several clusters.  This function outputs hatmu and the number of clusters. ","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso","category":"page"},{"location":"func/neariso/#IsoFuns.neariso","page":"Nearly isotonic regression","title":"IsoFuns.neariso","text":"neariso(x::Vector, λ, w::Vector) -> y, K\n\nPerform weighted nearly isotonic regression \n\nArguments\n\nx: input vector\nλ: parameter\nw: weights, default to ones if not provided\n\nKeywords\n\nvar: variances of the Gaussian distributions\n\nOutputs\n\ny: output vector, which is piecewise monotone\nK: the number of clusters\n\nAlgorithm\n\nmodified PAVA \n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"The number of clusters decreases if one chooses a bigger lambda. The function neariso_path outputs the lambda values at which several clusters merge, and the corresponding number of clusters.","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_path","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_path","page":"Nearly isotonic regression","title":"IsoFuns.neariso_path","text":"neariso_path(x::Vector, w::Vector) -> knot, K\n\nλ-path of weighted nearly isotonic regression\n\nArguments\n\nx: input vector\nw: weights, default to ones if not provided\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nK : the number of clusters at each checkpoint\n\nAlgorithm\n\nmodified PAVA \n\n\n\n\n\n","category":"function"},{"location":"func/neariso/#Generalized-nearly-isotonic-regression","page":"Nearly isotonic regression","title":"Generalized nearly isotonic regression","text":"","category":"section"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"The above formulation can be generalized to one parameter exponential families, as is done for the isotonic regression. Such a formulation is called the generalized nearly isotonic regression.","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"When X_isim p_i(x_imid theta_i) for i=1dotsn, where p_i(x_imidtheta_i ) is a one parameter exponential family defined by","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"p_i(x_imid theta_i) = h_i(x_i) exp (theta_i x_i -w_ipsi(theta_i))","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"the problem ","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"hattheta = argmin_theta_1 dots  theta_n left(sum_i=1^n -log p_i(x_imid theta_i)+lambdasum_i=1^n-1 (theta_i-theta_i+1)_+right)","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"can also be solved by the function neariso, as the formulation is equivalent to","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"hateta = argmin_eta_1dots  eta_nleft(sum_i=1^nfrac12 w_i (x_i-eta_i)^2+lambdasum_i=1^n-1 (eta_i-eta_i+1)_+right) quad textwith quad eta = mathrmE_theta x = psi ^prime (theta)","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"Further, apart from neariso_DistibutionName and neariso_path_DistributionName, the function neariso_AIC_value_DistributionName provides the AIC value at lambda. Additionally, neariso_AIC_DistributionName outouts the collection of the AIC values at each checkpoint of the relaxation parameter. Here, the Akaike Information Criterion (AIC) is defined by","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"mathrmAIC(lambda) = -2 sum_i=1^n log p_i(x_imid (hattheta_lambda)_i) + 2K_lambda","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"where hattheta_lambda = (hattheta_1dotshattheta_n) is the solution for lambda, and K_lambda signifies the corresponding number of clusters.","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_Normal","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_Normal","page":"Nearly isotonic regression","title":"IsoFuns.neariso_Normal","text":"neariso_Normal(x::Vector,λ,variance::Vector) -> y, K\n\nPerform weighted nearly isotonic regression  (Normal)\n\nArguments\n\nx: input vector\nλ: parameter\nvariance: weights, default to ones if not provided\n\nOutputs\n\ny: output vector, which is piecewise monotone\nK: the number of clusters\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_path_Normal","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_path_Normal","page":"Nearly isotonic regression","title":"IsoFuns.neariso_path_Normal","text":"neariso_path_Normal(x::Vector, variance::Vector) -> knot, K\n\nλ-path of weighted nearly isotonic regression (Normal)\n\nArguments\n\nx: input vector\nvariance a vector consisting of variances, default to ones if not provided\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nK : the number of clusters at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_Normal","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_Normal","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_Normal","text":"neariso_AIC_Normal(x::Vector, variance::Vector) -> knot, AIC\n\nλ-path and AIC of weighted nearly isotonic regression (Normal)\n\nArguments\n\nx: input vector\nw: weights, default to ones if not provided\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nAIC : AIC at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_value_Normal","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_value_Normal","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_value_Normal","text":"neariso_AIC_value_Normal(x::Vector, λ, variance::Vector) -> AIC\n\nAIC value of weighted nearly isotonic regression (Normal)\n\nArguments\n\nx: input vector\nw: weights, default to ones if not provided\n\nOutputs\n\nAIC : AIC at λ\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_Binomial","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_Binomial","page":"Nearly isotonic regression","title":"IsoFuns.neariso_Binomial","text":"nieariso_Binomial(success::Vector, λ, trial::Vector) -> y, K\n\nPerform weighted nearly isotonic regression  (Binomial)\n\nArguments\n\nsuccess: input vector consisting of the number of success\nλ: parameter\ntrials: a vector consisting of the number of trials\n\nOutputs\n\ny: output vector, which is piecewise monotone\nK: the number of clusters\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_path_Binomial","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_path_Binomial","page":"Nearly isotonic regression","title":"IsoFuns.neariso_path_Binomial","text":"neariso_path_Binomial(success::Vector, trial::Vector) -> knot, K\n\nλ-path of weighted nearly isotonic regression (Binomial)\n\nArguments\n\nsuccess: input vector consisting of the number of success\ntrials: a vector consisting of the number of trials\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nK : the number of clusters at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_Binomial","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_Binomial","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_Binomial","text":"neariso_AIC_Binomial(success::Vector, trial::Vector) -> knot, AIC\n\nλ-path and AIC of weighted nearly isotonic regression (Binomial)\n\nArguments\n\nsuccess: input vector consisting of the number of success\ntrials: a vector consisting of the number of trials\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nAIC : AIC at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_value_Binomial","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_value_Binomial","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_value_Binomial","text":"neariso_AIC_value_Binomial(success::Vector, λ, trial::Vector) -> AIC\n\nAIC value of weighted nearly isotonic regression (Binomial)\n\nArguments\n\nsuccess: input vector consisting of the number of success\ntrials: a vector consisting of the number of trials\n\nOutputs\n\nAIC : AIC at λ\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_Poisson","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_Poisson","page":"Nearly isotonic regression","title":"IsoFuns.neariso_Poisson","text":"nieariso_Poisson(x::Vector, λ) -> y, K\n\nPerform weighted nearly isotonic regression  (Poisson)\n\nArguments\n\nx: input vector consisting of the number of events\nλ: parameter\n\nOutputs\n\ny: output vector, which is piecewise monotone\nK: the number of clusters\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_path_Poisson","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_path_Poisson","page":"Nearly isotonic regression","title":"IsoFuns.neariso_path_Poisson","text":"neariso_path_Poisson(x::Vector) -> knot, K\n\nλ-path of weighted nearly isotonic regression (Poisson)\n\nArguments\n\nx: input vector consisting of the number of events\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nK : the number of clusters at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_Poisson","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_Poisson","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_Poisson","text":"neariso_AIC_Poisson(x::Vector) -> knot, AIC\n\nλ-path and AIC of weighted nearly isotonic regression (Poisson)\n\nArguments\n\nx: input vector consisting of the number of events\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nAIC : AIC at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_value_Poisson","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_value_Poisson","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_value_Poisson","text":"neariso_AIC_value_Poisson(x::Vector, λ) -> AIC\n\nAIC value of weighted nearly isotonic regression (poisson)\n\nArguments\n\nx: input vector consisting of the number of events\n\nOutputs\n\nAIC : AIC at λ\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"","category":"page"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_Chisq","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_Chisq","page":"Nearly isotonic regression","title":"IsoFuns.neariso_Chisq","text":"neariso_Chisq(x::Vector, λ, d::Vector) -> x, K\n\nPerform weighted nearly isotonic regression  (Chisq)\n\nArguments\n\nx: input vector \nλ: parameter\nd: a vector consisting of the degrees of freedom\n\nOutputs\n\ny: output vector, which is piecewise monotone\nK: the number of clusters\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_path_Chisq","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_path_Chisq","page":"Nearly isotonic regression","title":"IsoFuns.neariso_path_Chisq","text":"neariso_path_Chisq(x::Vector, d::Vector) -> knot, K\n\nλ-path of weighted nearly isotonic regression (Chisq)\n\nArguments\n\nx: input vector \ndefault: a vector consisting of the degrees of freedom\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nK : the number of clusters at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_Chisq","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_Chisq","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_Chisq","text":"neariso_AIC_Chisq(x::Vector, d::Vector) -> knot, AIC\n\nλ-path and AIC of weighted nearly isotonic regression (Chisq)\n\nArguments\n\nx: input vector \nd: a vector consisting of the degrees of freedom\n\nOutputs\n\nknot: checkpoints of the relaxation parameter\nAIC : AIC at each checkpoint\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"neariso_AIC_value_Chisq","category":"page"},{"location":"func/neariso/#IsoFuns.neariso_AIC_value_Chisq","page":"Nearly isotonic regression","title":"IsoFuns.neariso_AIC_value_Chisq","text":"neariso_AIC_value_Chisq(x::Vector, λ, d::Vector) -> AIC\n\nAIC value of weighted nearly isotonic regression (Chisq)\n\nArguments\n\nx: input vector \nd: a vector consisting of the degrees of freedom\n\nOutputs\n\nAIC : AIC at λ\n\n\n\n\n\n","category":"function"},{"location":"func/neariso/","page":"Nearly isotonic regression","title":"Nearly isotonic regression","text":"","category":"page"}]
}
