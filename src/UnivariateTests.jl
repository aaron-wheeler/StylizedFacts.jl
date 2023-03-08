"""
    heavy_tails_test(r::AbstractVector) -> Bool

Compute the kurtosis and Jarque-Bera statistic to test for the presence of heavy tails in
the real-valued vector of returns `r`.

# Arguments
- `r::AbstractVector{T}`: the real-valued vector of returns

# Returns
- `Bool`: true if the sample exhibits heavy tails and false otherwise.

# References
- 
"""
function heavy_tails_test(r::AbstractVector{T}) where T<:Real

    # test for excess kurtosis -
    κ = StatsBase.kurtosis(r)
    excess_kurtosis = if κ > 0 true else false end
    
    # test for rejection of Jarque-Bera test -
    JB = JarqueBeraTest(r)
    reject_normality = if pvalue(JB) > 0.05 false else true end
    
    # collect results -
    exhibit_heavytail = if excess_kurtosis && reject_normality true else false end
    
    # return -
    return exhibit_heavytail
end

"""
    absence_lin_autocor_test(r::AbstractVector) -> Bool

Compute the autocorrelation function and estimate the Hurst exponent to test for the
absence of linear autocorrelation in the real-valued vector of returns `r`.

# Arguments
- `r::AbstractVector{T}`: the real-valued vector of returns

# Returns
- `Bool`: true if the sample exhibits the absence of linear autocorrelation and false
    otherwise.

# Notes
The test...

# References
- 
"""
function absence_lin_autocor_test(r::AbstractVector{T}) where T<:Real
    
    # test for rapid decay to zero of autocorrelation function -
    acor_returns = StatsBase.autocor(r)
    conf = 1.96 / sqrt(length(r))
    rapid_decay = if abs(mean(acor_returns[2:end])) < conf true else false end
    
    # estimate the Hurst exponent using R/S analysis
    H,_,_ = Hurst_RS(r)
    H_LB, H_UB = Hurst_confint(r)
    approx_Brownian = if H < H_UB && H > H_LB true else false end
    
    # TODO: estimate the Hurst exponent using Detrended Fluctuation Analysis (DFA)
    # ...
    
    # collect results -
    exhibit_absencelinautocorr = if rapid_decay && approx_Brownian true else false end
    
    # return -
    return exhibit_absencelinautocorr
end

"""
    volatility_clustering_test(r::AbstractVector) -> Bool

Compute the autocorrelation function and estimate the Hurst exponent to test for the
presence of volatility clustering in the real-valued vector of returns `r`.

# Arguments
- `r::AbstractVector{T}`: the real-valued vector of returns

# Returns
- `Bool`: true if the sample exhibits volatility clustering and false otherwise.

# Notes
The test...

# References
- 
"""
function volatility_clustering_test(r::AbstractVector{T}) where T<:Real

    # calculate price volatility -
    p_volatility = abs.(r / std(r))
    
    # test for slow decay of positive autocorrelation -
    acor_volatility = StatsBase.autocor(p_volatility)
    u_conf = 1.96 / sqrt(length(r))
    slow_decay = if mean(acor_volatility[2:end]) > 0 &&  any(acor_volatility[2:end] .> u_conf) true else false end
    
    # estimate the Hurst exponent using R/S analysis and check for persistence -
    H,_,_ = Hurst_RS(p_volatility)
    _, H_UB = Hurst_confint(r)
    persistent = if H > (H_UB) true else false end
    
    # TODO: estimate the Hurst exponent using Detrended Fluctuation Analysis (DFA)
    # ...
    
    # collect results -
    exhibit_volatilityclustering = if slow_decay && persistent true else false end
    
    # return -
    return exhibit_volatilityclustering
end

"""
    nonlinear_dependence_test(r::Vector) -> Bool

Compute nonlinear transforms to test for the presence of nonlinear dependence in the
real-valued vector of returns `r`.

# Arguments
- `r::AbstractVector{T}`: the real-valued vector of returns

# Returns
- `Bool`: true if the sample exhibits nonlinear dependence and false otherwise.

# Notes
The test...

# References
- 
"""
function nonlinear_dependence_test(r::AbstractVector{T}) where T<:Real

    # compute nonlinear functions -
    r_abs = abs_returns(r)
    r_sqr = sqr_returns(r)
    r_sqr_abs = sqr_abs_returns(r)
    r_cos = cos_returns(r)
    r_log1p_sqr = log1p_sqr_returns(r)
    
    # test for slow decay of positive autocorrelation -
    u_conf = 1.96 / sqrt(length(r))    
    test1(r) = mean((StatsBase.autocor(r))[2:end]) > 0 && any((StatsBase.autocor(r))[2:end] .> u_conf) ? true : false
#     slow_decay = test1(r_abs) && test1(r_sqr) && test1(r_sqr_abs) && test1(r_cos) && test1(r_log1p_sqr)
    slow_decay = test1(r_abs) || test1(r_sqr) || test1(r_sqr_abs) || test1(r_cos) || test1(r_log1p_sqr)
    
    # estimate the Hurst exponent using R/S analysis and check for persistence -
    _, H_UB = Hurst_confint(r)
    test2(r) = Hurst_RS(r)[1] > (H_UB) ? true : false
#     persistent = test2(r_abs) && test2(r_sqr) && test2(r_sqr_abs) && test2(r_cos) && test2(r_log1p_sqr)
    persistent = test2(r_abs) || test2(r_sqr) || test2(r_sqr_abs) || test2(r_cos) || test2(r_log1p_sqr)
    
    # collect results -
    exhibit_nonlineardependence = if slow_decay && persistent true else false end
    
    # return -
    return exhibit_nonlineardependence
end

"""
    longmemory_volume_test(u::Vector) -> Bool

Compute the autocorrelation function and estimate the Hurst exponent to test for the
presence of long memory in the real-valued vector of volume `u`.

# Arguments
- `u::AbstractVector{T}`: the real-valued vector of volume

# Returns
- `Bool`: true if the sample exhibits long memory and false otherwise.

# Notes
The test...

# References
- 
"""
function longmemory_volume_test(u::AbstractVector{T}) where T<:Real
    
    # test for slow decay of positive autocorrelation -
    acor_volume = StatsBase.autocor(u)
    u_conf = 1.96 / sqrt(length(u))
    slow_decay = if mean(acor_volume[2:end]) > 0 &&  any(acor_volume[2:end] .> u_conf) true else false end
        
    # estimate the Hurst exponent using R/S analysis and check for persistence -
    H = Hurst_RS(u)[1]
    _, H_UB = Hurst_confint(u)
    persistent = if H > (H_UB) true else false end
    
    # collect results -
    exhibit_longmemoryvolume = if slow_decay && persistent true else false end
    
    # return -
    return exhibit_longmemoryvolume
end

"""
    crosscor_volume_volatility_test(u::Vector, r::Vector)

Compute the cross-correlation function to test for the presence of volatility clustering
in the real-valued vector of volume `u` and real-valued vector of returns `r`.

# Arguments
- `u::AbstractVector{T}`: the real-valued vector of volume
- `r::AbstractVector{T}`: the real-valued vector of returns

# Returns
- `Bool`: true if the sample exhibits volatility clustering and false otherwise.

# Notes
The test...

# References
- 
"""
function crosscor_volume_volatility_test(u::AbstractVector{T}, r::AbstractVector{T}) where T<:Real
    
    # calculate price volatility -
    p_volatility = abs.(r / std(r))
    
    # test for cross-correlation between trading volume and volatility -
    crosscor_vol = StatsBase.crosscor(u, p_volatility)
    u_conf = 1.96 / sqrt(length(u))
    cor_present = if mean(crosscor_vol) > 0 && any(crosscor_vol .> u_conf) true else false end
    
    # collect results -
    exhibit_crosscorvolumevolatility = if cor_present true else false end
    
    # return -
    return exhibit_crosscorvolumevolatility
end

"""
    cond_heavy_tails_test(r::Vector) -> Bool

Compute the kurtosis and Jarque-Bera statistic of the volatility-corrected residual
time series to test for the presence of conditional heavy tails in the real-valued
vector of returns `r`.

# Arguments
- `r::AbstractVector{T}`: the real-valued vector of returns

# Returns
- `Bool`: true if the sample exhibits conditional heavy tails and false otherwise.

# Notes
The test...

# References
- 
"""
function cond_heavy_tails_test(r::AbstractVector{T}) where T<:Real

    # correct returns series for volatility clustering via GARCH-type model -
    garch1_1 = ARCHModels.fit(GARCH{1,1}, r; meanspec = NoIntercept)
    residuals_condvol = residuals(garch1_1) # standardized by conditional volatility
    
    # test for excess kurtosis -
    κ = StatsBase.kurtosis(residuals_condvol)
    excess_kurtosis = if κ > 0 true else false end
    
    # test for rejection of Jarque-Bera test -
    JB = JarqueBeraTest(residuals_condvol)
    reject_normality = if pvalue(JB) > 0.05 false else true end
    
    # collect results -
    exhibit_condheavytail = if excess_kurtosis && reject_normality true else false end
    
    # return -
    return exhibit_condheavytail
end

"""
    leverage_effect_test(r::Vector) -> Bool

Compute the cross-correlation of returns with subsequent squared returns to test for the
presence of the leverage effect in the real-valued vector of returns `y`.

# Arguments
- `r::AbstractVector{T}`: the real-valued vector of returns

# Returns
- `Bool`: true if the sample exhibits the leverage effect and false otherwise.

# Notes
The test...

# References
- 
"""
function leverage_effect_test(r::AbstractVector{T}) where T<:Real
    
    # estimate volatility -
    volatility_est = sqr_abs_returns(r)
    
    # test for cross-correlation of returns with subsequent squared returns -
    leverage = StatsBase.crosscor(volatility_est, r)
    l_conf = -1.96 / sqrt(length(r))
    effect_present = if mean(leverage) < 0 && any(leverage .< l_conf) true else false end
    
    # collect results -
    exhibit_leverageeffect = if effect_present true else false end
    
    # return -
    return exhibit_leverageeffect
end

"""
    univariate_test_summary(ticker_symbol_array::Vector{String},
                            price_data_dictionary::Dict{String,DataFrame},
                            return_data_dictionary::Dict{String,DataFrame};
                            cc::Bool=true,
                            print_msg::Bool=true)

Compute the result of each univariate stylized facts test for each ticker in
`ticker_symbol_array`.

# Arguments
- `ticker_symbol_array::Vector{String}`: list of ticker symbols to be tested
- `price_data_dictionary::Dict{String,DataFrame}`: dictionary with a ticker symbol
    key and corresponding value of a DataFrame with OHLC price and volume data
- `return_data_dictionary::Dict{String,DataFrame}`: dictionary with a ticker symbol
    key and corresponding value of a DataFrame with returns series data (see: 
    compute_return_array() method)

# Keywords
- `cc::Bool=true`: denotes the use of close-to-close returns series (set `cc=false` if
    open-to-close returns series)
- `print_msg::Bool=true`: print out a summary of the stylized fact tests results

# Returns
- `Tuple{DataFrame,DataFrame}`: the first DataFrame stores the result `true` if the ticker
    exhibits the stylized fact for that respective test column and false otherwise; the
    second DataFrame is a summary of the percentage of tickers that were able to reproduce
    the stylized fact of each respective test column.

# References
- 
"""
function univariate_test_summary(ticker_symbol_array::Vector{String},
                                price_data_dictionary::Dict{String,DataFrame},
                                return_data_dictionary::Dict{String,DataFrame};
                                cc::Bool=true, print_msg::Bool=false)

    # build the DataFrames -
    styfacts_df = DataFrame(ticker = ticker_symbol_array)
    SF_aggregate_results = DataFrame(    
        SF1_HT_pct = [0.0], # % of tickers that exhibits stylized fact 1 -> heavy tails
        SF2_ALA_pct = [0.0], # % of tickers that exhibits stylized fact 2 -> absence of linear autocorrelation
        SF3_VC_pct = [0.0], # % of tickers that exhibits stylized fact 3 -> volatility clustering
        SF4_NLD_pct = [0.0], # % of tickers that exhibits stylized fact 4 -> nonlinear dependence
        SF5_LMV_pct = [0.0], # % of tickers that exhibits stylized fact 5 -> long memory of volume
        SF6_GLA_pct = [0.0], # % of tickers that exhibits stylized fact 6 -> gain/loss assymmetry
        SF7_CVV_pct = [0.0], # % of tickers that exhibits stylized fact 7 -> cross-correlation of volume and volatility
        SF8_AG_pct = [0.0], # % of tickers that exhibits stylized fact 8 -> aggregational gaussianity
        SF9_CHT_pct = [0.0], # % of tickers that exhibits stylized fact 9 -> conditional heavy tails
        SF10_LE_pct = [0.0] # % of tickers that exhibits stylized fact 10 -> leverage effect
    )
    
    # instantiate variables -
    N = length(ticker_symbol_array)
    SF1_HT = fill(false, N) # stylized fact 1 -> heavy tails
    SF2_ALA = fill(false, N) # stylized fact 2 -> absence of linear autocorrelation
    SF3_VC = fill(false, N) # stylized fact 3 -> volatility clustering
    SF4_NLD = fill(false, N) # stylized fact 4 -> nonlinear dependence
    SF5_LMV = fill(false, N) # stylized fact 5 -> long memory of volume
#     SF6_GLA = fill(false, N) # stylized fact 6 -> gain/loss assymmetry
    SF7_CVV = fill(false, N) # stylized fact 7 -> cross-correlation of volume and volatility
#     SF8_AG = fill(false, N) # stylized fact 8 -> aggregational gaussianity
    SF9_CHT = fill(false, N) # stylized fact 9 -> conditional heavy tails
    SF10_LE = fill(false, N) # stylized fact 10 -> leverage effect

    # conduct tests for each ticker -
    i = 0;
    for ticker_symbol ∈ ticker_symbol_array

        # instantiate time series -
        p_open_series = price_data_dictionary[ticker_symbol][!, :open]
        p_close_series = price_data_dictionary[ticker_symbol][!, :close]
        volume_series = price_data_dictionary[ticker_symbol][!, :volume]
        returns_series = return_data_dictionary[ticker_symbol][!, :R]

        # compute result -
        reproduce_heavytails = heavy_tails_test(returns_series)
        reproduce_absencelinautocorr = absence_lin_autocor_test(returns_series)
        reproduce_volatilityclustering = volatility_clustering_test(returns_series)
        reproduce_nonlineardependence = nonlinear_dependence_test(returns_series)
        reproduce_longmemoryvolume = longmemory_volume_test(volume_series)
#         reproduce_gainlossasym = GainLossAsymTest(returns_series)
        if cc == true
            reproduce_crosscorvolumevolatility = crosscor_volume_volatility_test(volume_series[1:end-1], returns_series)
        else
            reproduce_crosscorvolumevolatility = crosscor_volume_volatility_test(volume_series, returns_series)
        end
    #     reproduce_agggauss = AggGaussTest(p_open_series, p_close_series, :minute)
        reproduce_condheavytails = cond_heavy_tails_test(returns_series)
        reproduce_leverageeffect = leverage_effect_test(returns_series)

        # store result -
        i += 1
        SF1_HT[i] = reproduce_heavytails
        SF2_ALA[i] = reproduce_absencelinautocorr
        SF3_VC[i] = reproduce_volatilityclustering
        SF4_NLD[i] = reproduce_nonlineardependence
        SF5_LMV[i] = reproduce_longmemoryvolume
#         SF6_GLA[i] = reproduce_gainlossasym
        SF7_CVV[i] = reproduce_crosscorvolumevolatility
    #     SF8_AG[i] = reproduce_agggauss
        SF9_CHT[i] = reproduce_condheavytails
        SF10_LE[i] = reproduce_leverageeffect
    end

    # store results -
    styfacts_df[!, :SF1_HT] = SF1_HT
    SF_aggregate_results.SF1_HT_pct[1] = round((count(styfacts_df.SF1_HT) / N), digits=3)
    styfacts_df[!, :SF2_ALA] = SF2_ALA
    SF_aggregate_results.SF2_ALA_pct[1] = round((count(styfacts_df.SF2_ALA) / N), digits=3)
    styfacts_df[!, :SF3_VC] = SF3_VC
    SF_aggregate_results.SF3_VC_pct[1] = round((count(styfacts_df.SF3_VC) / N), digits=3)
    styfacts_df[!, :SF4_NLD] = SF4_NLD
    SF_aggregate_results.SF4_NLD_pct[1] = round((count(styfacts_df.SF4_NLD) / N), digits=3)
    styfacts_df[!, :SF5_LMV] = SF5_LMV
    SF_aggregate_results.SF5_LMV_pct[1] = round((count(styfacts_df.SF5_LMV) / N), digits=3)
#     styfacts_df[!, :SF6_GLA] = SF6_GLA
#     SF_aggregate_results.SF6_GLA_pct[1] = round((count(styfacts_df.SF6_GLA) / N), digits=3)
    styfacts_df[!, :SF7_CVV] = SF7_CVV
    SF_aggregate_results.SF7_CVV_pct[1] = round((count(styfacts_df.SF7_CVV) / N), digits=3)
#     styfacts_df[!, :SF8_AG] = SF8_AG
#     SF_aggregate_results.SF8_AG_pct[1] = round((count(styfacts_df.SF8_AG) / N), digits=3)
    styfacts_df[!, :SF9_CHT] = SF9_CHT
    SF_aggregate_results.SF9_CHT_pct[1] = round((count(styfacts_df.SF9_CHT) / N), digits=3)
    styfacts_df[!, :SF10_LE] = SF10_LE
    SF_aggregate_results.SF10_LE_pct[1] = round((count(styfacts_df.SF10_LE) / N), digits=3)

    # print output -
    if print_msg == true
        println("1. Number of assets that reproduce heavy tails = $(count(styfacts_df.SF1_HT)) / $(N)")
        if count(styfacts_df.SF1_HT) < N
            println("    Assets that fail to produce heavy tails => ", styfacts_df[styfacts_df.SF1_HT .== false, :ticker])
        end
        println("\n2. Number of assets that reproduce absence of linear autocorrelation = $(count(styfacts_df.SF2_ALA)) / $(N)")
        if count(styfacts_df.SF2_ALA) < N
            println("    Assets that fail to produce absence of linear autocorrelation => ", styfacts_df[styfacts_df.SF2_ALA .== false, :ticker])
        end
        println("\n3. Number of assets that reproduce volatility clustering = $(count(styfacts_df.SF3_VC)) / $(N)")
        if count(styfacts_df.SF3_VC) < N
            println("    Assets that fail to produce volatility clustering => ", styfacts_df[styfacts_df.SF3_VC .== false, :ticker])
        end
        println("\n4. Number of assets that reproduce nonlinear dependence = $(count(styfacts_df.SF4_NLD)) / $(N)")
        if count(styfacts_df.SF4_NLD) < N
            println("    Assets that fail to produce nonlinear dependence => ", styfacts_df[styfacts_df.SF4_NLD .== false, :ticker])
        end
        println("\n5. Number of assets that reproduce long memory of volume = $(count(styfacts_df.SF5_LMV)) / $(N)")
        if count(styfacts_df.SF5_LMV) < N
            println("    Assets that fail to produce long memory of volume => ", styfacts_df[styfacts_df.SF5_LMV .== false, :ticker])
        end
#         println("\n6. Number of assets that reproduce gain/loss asymmetry = $(count(styfacts_df.SF6_GLA)) / $(N)")
#         if count(styfacts_df.SF6_GLA) < N
#             println("    Assets that fail to produce gain/loss asymmetry => ", styfacts_df[styfacts_df.SF6_GLA .== false, :ticker])
#         end
        println("\n7. Number of assets that reproduce cross-correlation of volume and volatility = $(count(styfacts_df.SF7_CVV)) / $(N)")
        if count(styfacts_df.SF7_CVV) < N
            println("    Assets that fail to produce cross-correlation of volume and volatility => ", styfacts_df[styfacts_df.SF7_CVV .== false, :ticker])
        end
#         println("\n8. Number of assets that reproduce aggregational gaussianity = $(count(styfacts_df.SF8_AG)) / $(N)")
#         if count(styfacts_df.SF8_AG) < N
#             println("    Assets that fail to produce aggregational gaussianity => ", styfacts_df[styfacts_df.SF8_AG .== false, :ticker])
#         end
        println("\n9. Number of assets that reproduce conditional heavy tails = $(count(styfacts_df.SF9_CHT)) / $(N)")
        if count(styfacts_df.SF9_CHT) < N
            println("    Assets that fail to produce conditional heavy tails => ", styfacts_df[styfacts_df.SF9_CHT .== false, :ticker])
        end
        println("\n10. Number of assets that reproduce the leverage effect = $(count(styfacts_df.SF10_LE)) / $(N)")
        if count(styfacts_df.SF10_LE) < N
            println("    Assets that fail to produce the leverage effect => ", styfacts_df[styfacts_df.SF10_LE .== false, :ticker])
        end
    end

    # return -
    return styfacts_df, SF_aggregate_results
end