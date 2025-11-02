# using Logging # Show all messages
# global_logger(ConsoleLogger(stderr, Logging.Debug))
# # Show only info and above
# global_logger(ConsoleLogger(stderr, Logging.Info))
# # Show only warnings and errors
# global_logger(ConsoleLogger(stderr, Logging.Warn))

# package
using eegfun
using GLMakie
# using CairoMakie
# using BenchmarkTools
# while preprocessing routine
eegfun.preprocess("pipeline.toml")
