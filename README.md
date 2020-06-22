
# Markov chain Monte Carlo

<p align="center"><img src="https://travis-ci.org/dschrempf/mcmc.svg?branch=master"/></p>

Sample from a posterior using Markov chain Monte Carlo.

At the moment, the library is tailored to the Metropolis-Hastings algorithm
since it covers most use cases. However, implementation of more algorithms is
planned in the future.


## Documentation

The source code contains detailed documentation about general concepts as well
as specific functions.


## Examples

Have a look at the [example MCMC analyses](https://github.com/dschrempf/mcmc/tree/master/mcmc-examples). They can be built with [Stack](https://docs.haskellstack.org/en/stable/README/) and are
attached to this repository.

    git clone https://github.com/dschrempf/mcmc.git
    stack build

For example, estimate the accuracy of an archer with

    stack exec archery

