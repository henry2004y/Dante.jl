# Dante

<p align="center">
  <a href="https://github.com/henry2004y/Dante.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/henry2004y/Dante.jl/CI">
  </a>
  <a href="https://codecov.io/gh/henry2004y/Dante.jl">
    <img src="https://img.shields.io/codecov/c/github/henry2004y/Vlasiator.jl">
  </a>
  <a href="https://henry2004y.github.io/Dante.jl/dev">
    <img src="https://img.shields.io/badge/docs-latest-blue">
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue">
  </a>
</p>

Finite volume ideal MHD model with structured mesh.

## Benchmark

Riemann test 1 without plotting:
```julia
4.859 ms (1182 allocations: 691.80 KiB)
```

## Author

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

It originates from the [MATLAB version](https://github.com/henry2004y/FVMHD-Dante), and is then rewritten into Julia with improved performance and capabilities.