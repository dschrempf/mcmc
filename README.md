

# Markov chain Monte Carlo sampler

<p align="center"><img src="https://travis-ci.org/dschrempf/mcmc.svg?branch=master"/></p>

Sample from a posterior using Markov chain Monte Carlo (MCMC) algorithms.

At the moment, the following algorithms are available:

-   Metropolis-Hastings-Green (Geyer, Charles J, 2011);
-   Metropolis-coupled Markov chain Monte Carlo (also known as parallel
    tempering) (Geyer, Charles J, 1991,  Altekar, Gautam and Dwarkadas, Sandhya and Huelsenbeck, John P and Ronquist, Fredrik, 2004);
-   Hamilton Monte Carlo proposal (Neal, Radford M, 2011);
-   No U-Turn Sampler (NUTS) (Matthew D. Hoffman and Andrew Gelman, 2014).


## Documentation

The [source code](https://hackage.haskell.org/package/mcmc/docs/Mcmc.html) contains detailed documentation about general concepts as well
as specific functions.


## Examples

The Git repository also includes [example MCMC analyses](https://github.com/dschrempf/mcmc/tree/master/mcmc-examples). Build them with
[cabal-install](https://cabal.readthedocs.io/en/latest/cabal-commands.html#) or [Stack](https://docs.haskellstack.org/en/stable/README/).

    git clone https://github.com/dschrempf/mcmc.git
    cd mcmc
    stack build

For example, estimate the [accuracy of an archer](https://github.com/dschrempf/mcmc/blob/master/mcmc-examples/Archery/Archery.hs) with

    stack exec archery

For a more involved example, have a look at a [phylogenetic dating project](https://github.com/dschrempf/mcmc-dating).


# References

Altekar, Gautam and Dwarkadas, Sandhya and Huelsenbeck, John P and Ronquist, Fredrik (2004). *Parallel metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference*.

Geyer, Charles J (2011). *{Introduction to Markov Chain Monte Carlo}*, CRC press.

Geyer, Charles J (1991). *Markov chain Monte Carlo maximum likelihood*.

Matthew D. Hoffman and Andrew Gelman (2014). *The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo*.

Neal, Radford M (2011). *{MCMC Using Hamiltonian Dynamics}*, CRC press.

