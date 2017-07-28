# MCMC for Normalized Random Measure Mixture Models

Unix | CodeCov | License
---- | ------- | -------
[![Travis](https://travis-ci.org/BigBayes/BNPMix.jl.svg?branch=master)](https://travis-ci.org/BigBayes/BNPMix.jl) | [![CodeCov](http://codecov.io/gh/BigBayes/BNPMix.jl/coverage.svg?branch=master)](https://codecov.io/gh/BigBayes/BNPMix.jl?branch=master) | [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## What's this about

This package is a Julia port of the Java code released with the article *MCMC for Normalized Random Measure Mixture Models*.
Do not hesitate to create pull requests for enhancements or to open issues. You can write to me at: emile.mathieu -at- stats -dot- ox -dot- ac -dot- uk

## Installation and requirements

Requirements:

* Julia in `[0.6.x]`
* 64-bit architecture

In the Julia REPL:

```julia
Pkg.clone("https://github.com/BigBayes/BNPMix.jl")
using BNPMix
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

You can run the Reuse conditional sampler on the galaxy dataset modeled as a mixtures of Gaussian with a normalized generalized Gamma prior:

```bash
cd examples
julia galaxy.jl
```

## Tests

To localy run the tests, run in Julia:

```julia
Pkg.test("BNPMix")
```

## Reference

*  [Favaro, Teh, *MCMC for Normalized Random Measure Mixture Models*, 2013](https://www.stats.ox.ac.uk/~teh/research/npbayes/FavTeh2013a.pdf)
