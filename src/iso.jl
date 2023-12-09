export iso!, iso, iso_Normal!, iso_Normal, iso_Binomial!, iso_Binomial, iso_Poisson!, iso_Poisson, iso_Chisq!, iso_Chisq

"""
    iso!(x::Vector, w::Union{Vector, Nothing}=nothing) -> x

Perform isotonic regression 

# Arguments
- `x`: input vector 
- `w`: weights, default to ones if not provided

# Outputs
- `x`: output vector, which is monotone

# Algorithm
- PAVA (Pool-Adjacent-Violators algorithm)
"""
function iso!(x::Vector{<:AbstractFloat}, w::Union{Nothing,Vector{<:Real}}=nothing)
    
    n = length(x)

    if typeof(x) == Vector{Int64}
        x = map(Float64, x)
    end

    # Basic checks
    if n == 1
        return x
    end

    if w !== nothing
        length(w) == n || throw(DimensionMismatch("Lengths of input vector and weights mismatch"))
        minimum(w) > 0 || throw(DomainError(w, "Every element of w must be positive"))
    end

  
    @inbounds begin
        n -= 1
        while true
            i = 1
            pooled = 0
            while i <= n
                k = i
                while k <= n && x[k] >= x[k+1]
                    k += 1
                end

                if x[i] != x[k]
                    numerator = 0.0
                    denominator = 0.0
                    for j in i : k
                        wj = w === nothing ? 1.0 : w[j]
                        numerator += x[j] * wj
                        denominator += wj
                    end

                    for j in i : k
                        x[j] = numerator / denominator
                    end
                    pooled = 1
                end
                i = k + 1
            end
            if pooled == 0
                break
            end
        end
    end
    return x
end

"""
    iso(x::Vector) -> y

Same as ```iso!``` with ```w=Nothing```, but allocates an output vector ```y```.
"""
iso(x::Vector) = iso!(copy(x))

"""
    iso(x::Vector, w::Vector) -> y

Same as ```iso!```, but allocates an output vector ```y```.
"""
iso(x::Vector, w::Vector) = iso!(copy(x), w)


"""
    iso_Normal!(x::Vector, variance::Vector) -> x

Perform isotonic regression (Normal)

# Arguments
- `x`: input vector 
- `variance` a vector consisting of variances, default to ones if not provided

# Outputs
- `x`: output vector, which is monotone

"""
iso_Normal!(x::Vector, variance::Vector) = iso!(x, 1 ./variance)
iso_Normal!(x::Vector) = iso!(x)


"""
    iso_Normal(x::Vector, variance::Vector) -> y

Same as ```iso_Normal!```, but allocates an output vector ```y```.
"""
iso_Normal(x::Vector, variance::Vector) = iso(x, 1 ./variance)
iso_Normal(x::Vector) = iso(x)

"""
    iso_Binomial!(success::Vector, trial::Vector) -> success

Perform isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `success`: output vector, which is monotone

"""
iso_Binomial!(success::Vector, trial::Vector) = iso!(success./trial, trial)

"""
    iso_Binomial(success::Vector, trial::Vector) -> x

Same as ```iso_Binomial!```, but allocates an output vector ```x```.
"""
iso_Binomial(success::Vector, trial::Vector) = iso(success./trial, trial)

"""
    iso_Poisson!(x::Vector) -> x

Perform isotonic regression (Poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `x`: output vector, which is monotone

"""
iso_Poisson!(x::Vector) = iso!(Vector{Float64}(x))

"""
    iso_Poisson(x::Vector) -> y

Same as ```iso_Poisson!```, but allocates an output vector ```y```.
"""
iso_Poisson(x::Vector) = iso(Vector{Float64}(x))

"""
    iso_Chisq!(x::Vector, d::Vector) -> x

Perform isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `x`: output vector, which is monotone

"""
iso_Chisq!(x::Vector, d::Vector) = iso!(x./d, d/2)

"""
    iso_Chisq(x::Vector, d::Vector) -> y

Same as ```iso_Chisq!```, but allocates an output vector ```y```.
"""
iso_Chisq(x::Vector, d::Vector) = iso(x./d, d/2)
