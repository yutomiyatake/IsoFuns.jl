export neariso, neariso_path, neariso_Normal, neariso_path_Normal, neariso_AIC_Normal, neariso_AIC_value_Normal, neariso_Binomial, neariso_path_Binomial, neariso_AIC_Binomial, neariso_AIC_value_Binomial, neariso_Poisson, neariso_path_Poisson, neariso_AIC_Poisson, neariso_AIC_value_Poisson, neariso_Chisq, neariso_path_Chisq, neariso_AIC_Chisq, neariso_AIC_value_Chisq

mutable struct abmlc
    a::Vector
    β::Vector
    m::Vector
    c::Vector
    λ_tmp::Float64
    t_min::Float64
    idx::Vector{Int64}
    K::Int64
end

function init(v::abmlc,x::Vector,w::Vector,n)
    (v.β[1], v.a[1], v.c[1]) = (x[1], w[1], 1)
    v.K = 1
    @inbounds for l = 2:n
        if x[l] == x[l - 1]
            v.a[v.K] = v.a[v.K] + w[l]
            v.c[v.K] = v.c[v.K] + 1
        else
            v.K = v.K + 1
            (v.β[v.K], v.a[v.K], v.c[v.K]) = (x[l], w[l], 1)
        end
    end
    (v.β, v.a, v.c) = (v.β[1:v.K], v.a[1:v.K], v.c[1:v.K])
end

function compute_t_min(v::abmlc)
    s = sign.(max.(v.β[1:v.K - 1] - v.β[2:v.K], zeros(v.K - 1)))
    v.m = ([0;s] - [s;0]) ./ v.a

    if v.K == 1
        v.t_min = Inf
    else
        # compute t_min
        t = (v.β[2:v.K] - v.β[1:v.K - 1]) ./ (v.m[1:v.K - 1] - v.m[2:v.K]) + v.λ_tmp * ones(v.K - 1);
        v.t_min = minimum(t[t .> v.λ_tmp])

        min_first_index = findfirst(isequal(v.t_min),t)
        min_first_value = v.β[min_first_index] + v.m[min_first_index]*(v.t_min-v.λ_tmp) 
        for i = 1:length(t)
            if t[i]-v.t_min < 1e-14 && abs(v.β[i] + v.m[i]*(t[i]-v.λ_tmp) - min_first_value) <= eps(min_first_value)
                t[i] = v.t_min
            end
        end

        v.idx = findall(x -> x == v.t_min, t)
    end

end

function update_a(v::abmlc)
    for i = 1:length(v.idx)
        v.a[v.idx[end-i+1]] = v.a[v.idx[end-i+1]] + v.a[v.idx[end-i+1]+1]
    end
    v.a = v.a[setdiff(1:end, v.idx + ones(Int64, length(v.idx)))]
end

function update_c(v::abmlc)
    for i = 1:length(v.idx)
        v.c[v.idx[end-i+1]] = v.c[v.idx[end-i+1]] + v.c[v.idx[end-i+1]+1]
    end
    v.c = v.c[setdiff(1:end, v.idx + ones(Int64, length(v.idx)))]
end

function update_β(v::abmlc)
    v.β = v.β + (v.t_min - v.λ_tmp) * v.m
    v.β = v.β[setdiff(1:end, v.idx + ones(Int64, length(v.idx)))]
end

function update_λ_tmp(v::abmlc)
    v.λ_tmp = v.t_min
end

function update_K(v::abmlc)
    v.K = v.K - length(v.idx)
end



"""
    neariso(x::Vector, λ, w::Vector) -> y, K

Perform weighted nearly isotonic regression 

# Arguments
- `x`: input vector
- `λ`: parameter
- `w`: weights, default to ones if not provided

# Keywords
- `var`: variances of the Gaussian distributions

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

# Algorithm
- modified PAVA 
"""
function neariso(x::Vector{<:AbstractFloat},λ, w::Vector{<:Real}=ones(length(x)); var=[])

    n = length(x)
    n == length(w) || throw(DimensionMismatch("Lengths of input vector and weights mismatch"))
    
    if var != []
        n == length(var) || throw(DimensionMismatch("Lengths of input vector and var mismatch"))
        minimum(var) > 0 || throw(DomainError(var, "Every element of var must be positive"))
        x .= x ./ var
        w .= w .* var
    end

    minimum(w) >= 0 || throw(DomainError(w, "Every element of w must be positive"))


    # if λ == 0
    #     return x, param.K
    # end
        
    param = abmlc(zeros(n),zeros(n),zeros(n),zeros(Int64,n),0.0,0.0,[],0)
    init(param,x,w,n)

    if λ == 0
        return x, param.K
    end

        while param.λ_tmp < λ
            compute_t_min(param)

            if λ <= param.t_min
                param.β = param.β + (λ - param.λ_tmp) * param.m
                output = []
                for l = 1:param.K
                    append!(output, param.β[l] * ones(param.c[l]))
                end
                return output, param.K
            end

            update_β(param)
            update_a(param)
            update_c(param)
            update_λ_tmp(param)
            update_K(param)

            if param.K == 1
                return param.β[1] * ones(param.c[1]), param.K
            end
        end

end

"""
    neariso_Normal(x::Vector,λ,variance::Vector) -> y, K

Perform weighted nearly isotonic regression  (Normal)

# Arguments
- `x`: input vector
- `λ`: parameter
- `variance`: weights, default to ones if not provided

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
neariso_Normal(x::Vector,λ,variance::Vector) = neariso(x,λ,1 ./variance)
neariso_Normal(x::Vector,λ) = neariso(x,λ)


"""
    nieariso_Binomial(success::Vector, λ, trial::Vector) -> y, K

Perform weighted nearly isotonic regression  (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `λ`: parameter
- `trials`: a vector consisting of the number of trials

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
neariso_Binomial(success::Vector, λ, trial::Vector) = neariso(success./trial, λ, trial)

"""
    nieariso_Poisson(x::Vector, λ) -> y, K

Perform weighted nearly isotonic regression  (Poisson)

# Arguments
- `x`: input vector consisting of the number of events
- `λ`: parameter

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
neariso_Poisson(x::Vector, λ) = neariso(Vector{Float64}(x), λ)

"""
    neariso_Chisq(x::Vector, λ, d::Vector) -> x, K

Perform weighted nearly isotonic regression  (Chisq)

# Arguments
- `x`: input vector 
- `λ`: parameter
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
neariso_Chisq(x::Vector, λ, d::Vector) = neariso(x./d, λ, d/2)


"""
    neariso_path(x::Vector, w::Vector) -> knot, K

λ-path of weighted nearly isotonic regression

# Arguments
- `x`: input vector
- `w`: weights, default to ones if not provided

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint

# Algorithm
modified PAVA 
"""
function neariso_path(x::Vector{<:AbstractFloat},w::Vector{<:Real}=ones(length(x)))

    n = length(x)
    n == length(w) || throw(DimensionMismatch("Lengths of input vector and weights mismatch"))

    param = abmlc(zeros(n),zeros(n),zeros(n),zeros(Int64,n),0.0,0.0,[],0)
    knot = []
    K = []
    init(param,x,w,n)

    while true
        append!(knot, param.λ_tmp)
        append!(K, param.K)
        compute_t_min(param)

        if param.t_min == Inf
            return knot, K
        end

        update_β(param)
        update_a(param)
        update_c(param)
        update_λ_tmp(param)
        update_K(param)

    end
end

"""
    neariso_path_Normal(x::Vector, variance::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Normal)

# Arguments
- `x`: input vector
- `variance` a vector consisting of variances, default to ones if not provided

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
neariso_path_Normal(x::Vector,variance::Vector) = neariso_path(x,1 ./variance)
neariso_path_Normal(x::Vector) = neariso_path(x)

"""
    neariso_path_Binomial(success::Vector, trial::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
neariso_path_Binomial(success::Vector, trial::Vector) = neariso_path(success./trial, trial)


"""
    neariso_path_Poisson(x::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
neariso_path_Poisson(x::Vector) = neariso_path(Vector{Float64}(x))


"""
    neariso_path_Chisq(x::Vector, d::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `default`: a vector consisting of the degrees of freedom

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
neariso_path_Chisq(x::Vector, d::Vector) = neariso_path(x./d, d/2)



"""
    neariso_AIC_Normal(x::Vector, variance::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Normal)

# Arguments
- `x`: input vector
- `w`: weights, default to ones if not provided

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Normal(x::Vector, variance::Vector)
    knot, K = neariso_path_Normal(x,variance)
    AIC = Vector{Float64}(2*K)
    n = length(x)

    for i=1:length(knot)
        η, _ = neariso_Normal(x,knot[i],variance)
        for j=1:n
            AIC[i] += -2 * (-(x[j] - η[j])^2 / (2 * variance[j]) - 1/2 * log(2 * π * variance[j]))
        end
    end

    return knot, AIC
end
neariso_AIC_Normal(x::Vector) = neariso_AIC_Normal(x::Vector,ones(length(x)))

"""
    neariso_AIC_value_Normal(x::Vector, λ, variance::Vector) -> AIC

AIC value of weighted nearly isotonic regression (Normal)

# Arguments
- `x`: input vector
- `w`: weights, default to ones if not provided

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Normal(x::Vector, λ, variance::Vector)
    η, K = neariso_Normal(x,λ,variance)
    AIC = 2*K
    n = length(x)

    for j=1:n
        AIC += -2 * (-(x[j] - η[j])^2 / (2 * variance[j]) - 1/2 * log(2 * π * variance[j]))
    end

    return AIC
end
neariso_AIC_value_Normal(x::Vector, λ) = neariso_AIC_value_Normal(x::Vector, λ, ones(length(x)))


"""
    neariso_AIC_Binomial(success::Vector, trial::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Binomial(success::Vector, trial::Vector)
    knot, K = neariso_path_Binomial(success,trial)
    AIC = Vector{Float64}(2*K)
    n = length(success)

    for i=1:length(knot)
        η, _ = neariso_Binomial(success,knot[i],trial)
        for j=1:n
            if η[j]!=0. && η[j]!=1.
                AIC[i] += -2 * (success[j]*log(η[j]) + (trial[j]-success[j])*log(1-η[j]) + log_binomial(trial[j],success[j]))
            end
        end
    end

    return knot, AIC
end

function log_binomial(n, k)
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
end


"""
    neariso_AIC_value_Binomial(success::Vector, λ, trial::Vector) -> AIC

AIC value of weighted nearly isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Binomial(success::Vector, λ, trial::Vector)
    η, K = neariso_Binomial(success,λ,trial)
    AIC = 2*K
    n = length(success)

    for j=1:n
        if η[j]!=0. && η[j]!=1.
            AIC += -2 * (success[j]*log(η[j]) + (trial[j]-success[j])*log(1-η[j]) + log_binomial(trial[j],success[j]))
        end
    end

    return AIC
end

"""
    neariso_AIC_Poisson(x::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Poisson(x::Vector)
    knot, K = neariso_path_Poisson(x)
    AIC = Vector{Float64}(2*K)
    n = length(x)

    for i=1:length(knot)
        η, _ = neariso_Poisson(x,knot[i])
        for j=1:n
            AIC[i] += -2 * (x[j]*log(η[j]) - η[j] - lgamma(x[j]+1))
        end
    end

    return knot, AIC
end

"""
    neariso_AIC_value_Poisson(x::Vector, λ) -> AIC

AIC value of weighted nearly isotonic regression (poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Poisson(x::Vector, λ)
    η, K = neariso_Poisson(x,λ)
    AIC = 2*K
    n = length(x)

    for j=1:n
        AIC += -2 * (x[j]*log(η[j]) - η[j] - lgamma(x[j]+1))
    end

    return AIC
end

"""
    neariso_AIC_Chisq(x::Vector, d::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Chisq(x::Vector, d::Vector)
    knot, K = neariso_path_Chisq(x,d)
    AIC = Vector{Float64}(2*K)
    n = length(x)

    for i=1:length(knot)
        η, _ = neariso_Chisq(x,knot[i],d)
        for j=1:n
            AIC[i] += -2 * ((d[j]/2-1)*x[j] - x[j]/η[j] - lgamma(d[j]/2) - d[j]/2*log(η[j]))
        end
    end

    return knot, AIC
end

"""
    neariso_AIC_value_Chisq(x::Vector, λ, d::Vector) -> AIC

AIC value of weighted nearly isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Chisq(x::Vector, λ, d::Vector)
    η, K = neariso_Chisq(x,λ,d)
    AIC = 2*K
    n = length(x)

    for j=1:n
        AIC += -2 * ((d[j]/2-1)*x[j] - x[j]/η[j] - lgamma(d[j]/2) - d[j]/2*log(η[j]))
    end

    return AIC
end

