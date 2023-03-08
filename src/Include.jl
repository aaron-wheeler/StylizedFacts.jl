# setup internal paths -
_PATH_TO_SRC = dirname(pathof(@__MODULE__))

# load external packages that we depend upon -
using DataFrames
using Distributions
using StatsBase
using Statistics
using HypothesisTests
using ARCHModels
using SpecialFunctions

# load my codes -
include(joinpath(_PATH_TO_SRC, "UnivariateTests.jl"))
include(joinpath(_PATH_TO_SRC, "TestUtils.jl"))