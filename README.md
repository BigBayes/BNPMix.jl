# MCMC for Normalized Random Measure Mixture Models

Unix | CodeCov | License
---- | ------- | -------
[![Travis](https://travis-ci.org/emilemathieu/NRMMM.jl.svg?branch=master)](https://travis-ci.org/emilemathieu/NRMMM.jl) | [![CodeCov](http://codecov.io/gh/emilemathieu/NRMMM.jl/coverage.svg?branch=master)](https://codecov.io/gh/emilemathieu/NRMMM.jl?branch=master) | [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## What's this about

This package is a Julia port of the Java code released with the article *MCMC for Normalized Random Measure Mixture Models*.
Do not hesitate to create pull requests for enhancements or to open issues.

## Installation and requirements

Requirements:

* Julia in `[0.6.x]`
* 64-bit architecture

In the Julia REPL:

```julia
Pkg.clone("https://github.com/emilemathieu/NRMMM.jl")
using NRMMM
```

## Algorithms implemented

* Marginalized Samplers:
  * Neal’s Algorithm 8 generalized
  * The Reuse algorithm
* Conditional Slice Sampler

**Notes**
- Conjugate version not implemented
- Only Normal-Gamma-Independent emission implemented
- Multidimensional data observation not implemented

## Example

You can run the Reuse conditional sampler on the galaxy dataset modelled as a mixtures of Gaussian with a normalized generalized Gamma prior:

```bash
cd examples
julia galaxy.jl
```

## Reference

* Favaro, Teh, *MCMC for Normalized Random Measure Mixture Models*, [link](https://www.stats.ox.ac.uk/~teh/research/npbayes/FavTeh2013a.pdf), 2013.
