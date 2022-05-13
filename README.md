
# Markov chain Monte Carlo sampler

<p align="center"><img src="https://travis-ci.org/dschrempf/mcmc.svg?branch=master"/></p>

Sample from a posterior using Markov chain Monte Carlo (MCMC) algorithms.

At the moment, the following algorithms are available:

-   Metropolis-Hastings-Green <sup><a id="fnr.1" class="footref" href="#fn.1" role="doc-backlink">1</a></sup>;
-   Metropolis-coupled Markov chain Monte Carlo (also known as parallel
    tempering) <sup><a id="fnr.2" class="footref" href="#fn.2" role="doc-backlink">2</a></sup> <sup>, </sup><sup><a id="fnr.3" class="footref" href="#fn.3" role="doc-backlink">3</a></sup>.
-   Hamilton Monte Carlo proposal <sup><a id="fnr.4" class="footref" href="#fn.4" role="doc-backlink">4</a></sup>.


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


# Footnotes

<sup><a id="fn.1" href="#fnr.1">1</a></sup> Geyer, C. J., Introduction to Markov chain Monte Carlo, In Handbook of
Markov Chain Monte Carlo (pp. 45) (2011). CRC press.

<sup><a id="fn.2" href="#fnr.2">2</a></sup> Geyer, C. J., Markov chain monte carlo maximum likelihood, Computing
Science and Statistics, Proceedings of the 23rd Symposium on the Interface,
(1991).

<sup><a id="fn.3" href="#fnr.3">3</a></sup> Altekar, G., Dwarkadas, S., Huelsenbeck, J. P., & Ronquist, F., Parallel
metropolis coupled markov chain monte carlo for bayesian phylogenetic inference,
Bioinformatics, 20(3), 407–415 (2004).

<sup><a id="fn.4" href="#fnr.4">4</a></sup> Neal, R. M., Mcmc Using Hamiltonian Dynamics, In S. Brooks, A. Gelman, G.
Jones, & X. Meng (Eds.), Handbook of Markov Chain Monte Carlo (2011). CRC press.
