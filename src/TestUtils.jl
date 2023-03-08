"""
    Hurst_RS(Z::AbstractVector) -> Tuple{Float64,Float64,Float64}

Estimate the Hurst exponent via rescaled range analysis.

# Arguments
- `Z::AbstractVector`: the real-valued time series vector

# Keywords
- `adjust::Bool = true`: correct the estimate via the Anis-Lloyd adjustment

# Returns
- `Tuple{Float64,Float64,Float64}`: the Hurst exponent in the first position, the
    log(n) values in the second position, and the log(R/S) values in the third position.

# Throws
- `ArgumentError`: throws an error if the length of `Z` is less than 500.

# Notes
For several different lengths n, divide the time series `Z`` of length N into d subseries of length n, where n
is an integral divisor of N. For each subseries m = 1, 2, ..., d, follow the steps presented in CITE.
per confint paper -> log2, no n < 50, at least 512 observations

# References
https://watermark.silverchair.com/63-1-111.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAvUwggLxBgkqhkiG9w0BBwagggLiMIIC3gIBADCCAtcGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMZn9ZfD4HT8p48oH7AgEQgIICqFFc432q0RZTSB95JYv8I7i-nrSWUQdoHy2CfmhMYjbyoUXNCkoCV1Hu_MbJC28IYbaz0Li3n2PcmUeZjBkqoqec_OZLQMypFp5PCR2xxqRqOUW3VQP2lrvPH-pJe9znb-RBrgY2wNF0zkvxadG1f7GAntRh2JBchmRWoCLKANBerIZnzKlaufIpRMOAq6uajOo36gzoNvlBGI2WghtJhsZuQqCFwHYRcTtmefCGIJRr4xOXs7D1ap9F5OwL8NhsRuZ8pPq_0YQuQctiSjU2DRmyD0nJI1BybRb-678jD06sgI0XH7Hu_CWoUAVSqYHjuX1gvdYac1LHZ3nIq4AsNmVSPJzr3cZ9-C4rwJ-kmbZCjO2TXavD28xpLcwFuOzhjrwoQGLfhLFvYiyjJWLyUMnH4sFDS2QKR39IVX37qpN_6B3xipuowiAeY13dyzWw_zE3TjIGBumESD5UHOu3Q05yMi7DHlGZq-R9G2SxYJB4PnOOqgWvR7uMgCGw1Xjg0jWs9UrpTcDuQMTe0QoKjWUwjNTI5sodgI8_YBogBCCXankFtuqYSsiElfl5CI7CcWj33VX6H7sMv2zSE8j-XyuOZKGHA5VnvF4N2mjBupHbId7rh6olQ8ANub2jqKsv8colB2Zv6ygtVq0INA760Q3oduR2fl1TAc5Y-VxJ17GtWUjRPLkbgsBhITu6KpqRdQZIZnQIl7dUGh9TdLtxcuMkGEkLqSk3qsCQ02Vj7zhb5zwdk4L3dZ7UMU8Ov1YlUJExbkuRnDKJfTW_i9FOkDnYTbapjLbsWsbW4hJA_T1es65nEOfQVs2IaXu2Sbg6QLMiqHgPnVyauzMVFMhPtm-j2Coj0dc4zLFtRDct9aUxX5NNhrThVjwJurtXh9MhcKtWyAFMdq5u
https://arxiv.org/pdf/1805.08931.pdf
https://arxiv.org/pdf/cond-mat/0103510.pdf

"""
function Hurst_RS(Z::AbstractVector{T}; adjust::Bool=true) where T<:Real

    # check series length -
    N = length(Z)
    if N < 500 throw(ArgumentError("Size Error: series length N = $(N) < 500")) else nothing end
    
    # divide the data Z into various lengths n (for n > 50) and take log -
    log_n=Array(range(log2(51),stop=log2(N),step=0.25))
    
    # collect each length n to evaluate -
    n_list = []
    for x in log_n
        push!(n_list,round(Int64,exp2(x),RoundDown))
    end

    # ensure the full data length N is included in evaluation set -
    if !(length(Z) in n_list)
        push!(n_list, N)
        push!(log_n, log2(N))
    end

    # calculate expected value for R/S over various n -
    RS_n = []
    for n in n_list
        rs_m = []
        for start in (range(0,stop=N,step=n))
            if (start+n) > N
                break
            end
            # set the rolling subseries m within m = 1,2,...,d
            Z_m = Z[start+1:start+n]
            
            # calculate the mean of the subseries
            E_m = Statistics.mean(Z_m)
            
            # normalize the data by subtracting the sample mean
            deviations = Z_m .- E_m
            
            # sum the deviations from the mean
            Y_m = cumsum(deviations)
            
            # compute the range R
            R_m = maximum(Y_m) - minimum(Y_m)
            
            # compute the standard deviation S
            S_m = Statistics.std(Z_m)
            
            # calculate the rescaled range R/S
            if R_m==0 || S_m==0
                RS_part = 0
            else
                RS_part = R_m / S_m
            end
            
            # add RS of block
            if RS_part != 0
                push!(rs_m, RS_part)
            end
        end

        # store the mean RS value for all d subseries of length n -
        if length(rs_m)>0
            push!(RS_n, Statistics.mean(rs_m))
        end
    end
    
    # use the relation (R/S)_n ~ cn^H to compute an estimate of H -
    if adjust == true
        # correct the biased estimate via the Anis-Lloyd adjustment
        E_RR = []
        for n in eachindex(n_list)
            # for n <= 340
            if n_list[n] <= 340
                AL_adjust = ((n_list[n] - 0.5)/n_list[n]) * ((gamma(((n_list[n] - 1)/2))) / 
                    (sqrt(pi)*gamma(n_list[n]/2))) * sum(sqrt((n_list[n] - i) / i) for i in 1:(n_list[n]-1))
                push!(E_RR, AL_adjust)
            # for n > 340
            else
                AL_adjust = ((n_list[n] - 0.5)/n_list[n]) * (1 / 
                    (sqrt(n_list[n]*(pi/2)))) * sum(sqrt((n_list[n] - i) / i) for i in 1:(n_list[n]-1))
                push!(E_RR, AL_adjust)
            end
        end

        # compute adjusted rescaled range statistic -
        RS_AL = RS_n .- E_RR .+ sqrt.(0.5 .* pi .* n_list)
        
        # run simple linear regression -> slope is estimate of the hurst exponent
        A=Array{Float64}([log_n ones(length(RS_AL))])
        RSlog=[]
        for r in RS_AL
            push!(RSlog,log2(r))
        end
        B=Array{Float64}(RSlog)
        H,c=A\B
        c=exp2(c)
        
        # collect results -
        return H, log_n, RSlog
    else
        # run simple linear regression -> slope is estimate of the hurst exponent
        A=Array{Float64}([log_n ones(length(RS_n))])
        RSlog=[]
        for r in RS_n
            push!(RSlog,log2(r))
        end
        B=Array{Float64}(RSlog)
        H,c=A\B
        c=exp2(c)

        # collect results -
        return H, log_n, RSlog
    end
end


"""
    Hurst_confint(Z::AbstractVector; level::Symbol=:ninetyfive) -> Tuple{Float64,Float64}

Returns the confidence interval for Hurst exponents estimated via rescaled range analysis.

Returns a Tuple with the lower and upper bound of the Hurst exponent (according to the
confidence interval).

# Arguments
- `Z::AbstractVector`: the real-valued time series vector

# Keywords
- `level::Symbol=:ninetyfive`: the confidence level (options -> `:ninety`,
    `:ninetyfive` (default), `:ninetynine`)

# Returns
- `Tuple{Float64,Float64,Float64}`: the lower bound of the Hurst exponent in the first
    position and the upper bound of the Hurst exponent in the second position (according to the
    confidence interval).

# Throws
- `ArgumentError`: throws an error if the length of `Z` is less than 500.

# References
https://arxiv.org/pdf/cond-mat/0103510.pdf
"""
function Hurst_confint(Z::AbstractVector{T}; level::Symbol=:ninetyfive) where T<:Real

    # check series length -
    N = length(Z)
    if N < 500 throw(ArgumentError("Size Error: series length N = $(N) < 500")) else nothing end
    
    # set adjusted sample size -
    M = log2(N)

    if level == :ninety
        # 90% confidence interval
        H_lowerconfint = 0.5 - exp(-7.35 * log(log(M)) + 4.06)
        H_upperconfint = exp(-7.07 * log(log(M)) + 3.75) + 0.5
    elseif level == :ninetyfive
        # 95% confidence interval
        H_lowerconfint = 0.5 - exp(-7.33 * log(log(M)) + 4.21)
        H_upperconfint = exp(-7.20 * log(log(M)) + 4.04) + 0.5
    elseif level == :ninetynine
        # 99% confidence interval
        H_lowerconfint = 0.5 - exp(-7.19 * log(log(M)) + 4.34)
        H_upperconfint = exp(-7.51 * log(log(M)) + 4.58) + 0.5
    else
        throw(ArgumentError("Invalid Argument: confidence interval $(level) not defined."))
    end
    
    # return -
    return H_lowerconfint, H_upperconfint
end

"""
    Hurst_confint(N::Int; level::Symbol=:ninetyfive) -> Tuple{Float64,Float64}

Returns the confidence interval for Hurst exponents estimated via rescaled range analysis.

# Arguments
- `N::Int`: the length of the time series

# Keywords
- `level::Symbol=:ninetyfive`: the confidence level (options -> `:ninety`,
    `:ninetyfive` (default), `:ninetynine`)

# Returns
- `Tuple{Float64,Float64,Float64}`: the lower bound of the Hurst exponent in the first
    position and the upper bound of the Hurst exponent in the second position (according to the
    confidence interval).

# References
https://arxiv.org/pdf/cond-mat/0103510.pdf
"""
function Hurst_confint(N::Int; level::Symbol=:ninetyfive)
    
    # set adjusted sample size -
    M = log2(N)

    if level == :ninety
        # 90% confidence interval
        H_lowerconfint = 0.5 - exp(-7.35 * log(log(M)) + 4.06)
        H_upperconfint = exp(-7.07 * log(log(M)) + 3.75) + 0.5
    elseif level == :ninetyfive
        # 95% confidence interval
        H_lowerconfint = 0.5 - exp(-7.33 * log(log(M)) + 4.21)
        H_upperconfint = exp(-7.20 * log(log(M)) + 4.04) + 0.5
    elseif level == :ninetynine
        # 99% confidence interval
        H_lowerconfint = 0.5 - exp(-7.19 * log(log(M)) + 4.34)
        H_upperconfint = exp(-7.51 * log(log(M)) + 4.58) + 0.5
    else
        throw(ArgumentError("Invalid Argument: confidence interval $(level) not defined."))
    end
    
    # return -
    return H_lowerconfint, H_upperconfint
end

# define some nonlinear functions -
abs_returns(returns::AbstractVector) = abs.(returns)
sqr_returns(returns::AbstractVector) = (returns).^2
sqr_abs_returns(returns::AbstractVector) = (abs_returns(returns)).^2
cos_returns(returns::AbstractVector) = cos.(returns)
log1p_sqr_returns(returns::AbstractVector) = log1p.(sqr_returns(returns))

"""
    compute_return_array(data::DataFrame; cc=true) -> DataFrame

Compute the return array of the price series DataFrame

# Arguments
- `data::DataFrame`: the price series DataFrame

# Keywords
- `cc::Bool=true`: denotes the use of close-to-close returns series (set `cc=false` if
    open-to-close returns series is desired)

# Returns
- `DataFrame`: the first column `:timestamp` denotes the timestamp and the second column
`:R` denotes the return instance corresponding to the timestamp
"""
function compute_return_array(data::DataFrame; cc=true)::DataFrame
    
    # compute log-returns -
    if cc == true
        p = data[:,:close]
        r = diff(log.(p[:,1])) # close-to-close returns
    else
        r = [log(data[i,"close"]) -
                            log(data[i,"open"]) for i in 1:nrow(data)] # open-to-close returns
    end 
    
    # build table -
    if cc == true
        t = data[1:(end - 1), :timestamp]
        return_table = DataFrame(timestamp = t, R = r)
    else
        t = data[:, :timestamp]
        return_table = DataFrame(timestamp = t, R = r)
    end

    # return -
    return return_table;
end