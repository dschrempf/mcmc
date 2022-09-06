
# Markov chain Monte Carlo sampler

<p align="center"><img src="https://travis-ci.org/dschrempf/mcmc.svg?branch=master"/></p>

Sample from a posterior using Markov chain Monte Carlo (MCMC) algorithms.

At the moment, the following algorithms are available:

-   Metropolis-Hastings-Green (<a href="#citeproc_bib_item_3">Geyer 2011</a>);
-   Metropolis-coupled Markov chain Monte Carlo (also known as parallel
    tempering) (<a href="#citeproc_bib_item_2">Geyer 1991</a>; <a href="#citeproc_bib_item_1">Altekar et al. 2004</a>);
-   Hamilton Monte Carlo proposal (<a href="#citeproc_bib_item_5">Neal 2011</a>);
-   No U-Turn Sampler (NUTS) (<a href="#citeproc_bib_item_4">Hoffman and Gelman 2014</a>).


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

<style>.csl-entry{text-indent: -1.5em; margin-left: 1.5em;}</style><div class="csl-bib-body">
  <div class="csl-entry"><a id="citeproc_bib_item_1"></a>Altekar, Gautam, Sandhya Dwarkadas, John P Huelsenbeck, and Fredrik Ronquist. 2004. “Parallel Metropolis Coupled Markov Chain Monte Carlo for Bayesian Phylogenetic Inference.” <i>Bioinformatics</i> 20 (3): 407–15.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_2"></a>Geyer, Charles J. 1991. “Markov Chain Monte Carlo Maximum Likelihood.” <i>Computing Science and Statistics, Proceedings of the 23rd Symposium on the Interface</i>.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_3"></a>———. 2011. “Introduction to Markov Chain Monte Carlo.” In <i>Handbook of Markov Chain Monte Carlo</i>, edited by Steve Brooks, Andrew Gelman, Galin Jones, and Xiao-Li Meng, 45. CRC press.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_4"></a>Hoffman, Matthew D., and Andrew Gelman. 2014. “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.” <i>Journal of Machine Learning Research</i> 15 (47): 1593–1623.</div>
  <div class="csl-entry"><a id="citeproc_bib_item_5"></a>Neal, Radford M. 2011. “MCMC Using Hamiltonian Dynamics.” In <i>Handbook of Markov Chain Monte Carlo</i>, edited by Steve Brooks, Andrew Gelman, Galin Jones, and Xiao-Li Meng. CRC press.</div>
</div>

