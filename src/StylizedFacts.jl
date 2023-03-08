module StylizedFacts

# include -
include("Include.jl")

# methods that we export -
export heavy_tails_test
export absence_lin_autocor_test
export volatility_clustering_test
export nonlinear_dependence_test
export longmemory_volume_test
export crosscor_volume_volatility_test
export cond_heavy_tails_test
export leverage_effect_test
export univariate_test_summary
export Hurst_RS
export Hurst_confint
export compute_return_array

end # module