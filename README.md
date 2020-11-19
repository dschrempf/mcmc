
# Markov chain Monte Carlo

<p align="center"><img src="https://travis-ci.org/dschrempf/mcmc.svg?branch=master"/></p>

Sample from a posterior using Markov chain Monte Carlo (MCMC) methods.

At the moment, the library is tailored to the Metropolis-Hastings-Green
algorithm since it covers most use cases. More algorithms will be implemented
soon.


## Documentation

The [source code](https://hackage.haskell.org/package/mcmc) contains detailed documentation about general concepts as well
as specific functions.


## Examples

Have a look at the [example MCMC analyses](https://github.com/dschrempf/mcmc/tree/master/mcmc-examples). They can be built with [Stack](https://docs.haskellstack.org/en/stable/README/) and are
attached to this repository.

    git clone https://github.com/dschrempf/mcmc.git
    cd mcmc
    stack build

For example, estimate the [accuracy of an archer](https://github.com/dschrempf/mcmc/blob/master/mcmc-examples/Archery/Archery.hs) with

    stack exec archery

