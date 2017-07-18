# MCMC for Normalized Random Measure Mixture Models

Unix | CodeCov | License
---- | ------- | -------
[![Travis](https://travis-ci.org/emilemathieu/NRMMM.jl.svg?branch=master)](https://travis-ci.org/emilemathieu/NRMMM.jl) |  | [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Installation and requirements

Requirements:

* Julia in `[0.6.x]`
* 64-bit architecture

In the Julia REPL:

```julia
Pkg.clone("https://github.com/emilemathieu/NRMMM.jl")
using NRMMM
```

## Example

You can run the Reuse conditional sampler on the galaxy dataset modelled as a mixtures of Gaussian with a normalized generalized Gamma prior:

```bash
cd examples
julia galaxy.jl
```

## Reference

* Favaro, Teh, *MCMC for Normalized Random Measure Mixture Models*, [link](https://www.stats.ox.ac.uk/~teh/research/npbayes/FavTeh2013a.pdf), 2013.
