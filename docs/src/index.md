# Dante.jl

Dante.jl is a toy finite volume ideal MHD model on structured mesh. It is the descendant of [the original version in MATLAB](https://github.com/henry2004y/FVMHD-Dante).

## Benchmark

Riemann test 1 (Sods problem) without plotting (`test/PARAM_test.toml`):
```julia
4.859 ms (1182 allocations: 691.80 KiB)
```
