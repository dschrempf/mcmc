
# Markov chain Monte Carlo sampling - ChangeLog


## Unreleased changes


## 0.8.0.1

-   Improve exception handling (also during execution of monitors; also improve
    the output).
-   Increase maximum tuning parameter.
-   Minor bug fix: acceptance rates are again reported at end of run.


## 0.8.0.0

-   Use the new bytestring realfloat builder.
-   Fix a small memory leak during burn in with auto tuning.
-   Save chain after burn in.
-   New example (Poisson; see the respective [blog post](https://dschrempf.github.io/coding/2022-06-28-sample-from-a-posterior-using-markov-chain-monte-carlo-algorithms-and-haskell/)).
-   Simplify prior, likelihood, posterior, and Jacobian types.
-   Do not export `grabble`.
-   Remove chain index.
-   Remove `lengthT`.
-   Remove `ProposalOrder` default instance.
-   Use GHC 9.2.4 by default.
-   Reduce package dependencies as much as possible.


## 0.7.0.1

-   Use random-1.2. This means, that the random number generator type has changed
    from `GenIO` (`System.Random.MWC`) to `IOGenM StdGen`
    (`System.Random.Stateful`).
-   Always try to save chain, not only on `UserInterupt`.
-   Hamiltonian proposals:
    -   Split `HSettings` into `HParams` and `HStructure` (preparation for NUTS).
    -   Rename `HTune` to `HTuningConf`.
    -   No-U-Turn sampler (NUTS).
-   Significant changes to `Proposal` and `Tuner`.
    -   `PResult`.
    -   `ProposalSimple` -> `PFunction`; Return type is `(PResult, Maybe AcceptanceCounts)`.


## 0.6.2.5

-   Provide unsafe loader for MHG algorithm (useful for initializing chains with
    different prior or likelihood functions from saves).
-   Documentation and readme.
-   Log normal prior distribution.
-   Parallel computation of prior and likelihood. Speculative parallelization;
    this change is not always beneficial, we will see.


## 0.6.2.4

-   Specify covariance version bounds. Use `covariance-0.2.0.0` (specifically
    state that sigma is rescaled with `rescaleSWith`).


## 0.6.2.3

-   Allow burn in with fast proposals only (`BurnInWithCustomAutoTUning`).
    Sometimes it is advantageous to hold back slow proposals initially, especially
    when the state is so far off that it does not make sense to compute complex
    proposals.
-   Hamiltonian proposal: Use automatic differentiation specialized to `Double`
    (roughly 10 percent faster).


## 0.6.2.2

-   Remove dependency `monad-parallel`. Fix stackage build.


## 0.6.2.1

-   Improve `logGammaG`. The calculation of the gamma function involves vectors;
    the generalized version needs boxed vectors, and is slow. Using `Typeable`, we
    now check if the type is `Double`, and use the fast version in this case.


## 0.6.2.0

-   Improve leapfrog integrator.
-   Update tooling.
-   Cleanup proposals.


## 0.6.1.0

-   Revamp Hamiltonian proposal (storable vectors).
-   Use mass matrices; allow tuning of all masses (covariance estimation using
    specialized estimators).


## 0.6.0.0

-   Improve documentation.
-   Generalized priors allowing automatic differentiation.
-   Hamiltonian proposal.


### mcmc-tree

-   Moved to another repository: <https://github.com/dschrempf/mcmc-date>.


## 0.5.0.0

-   Marginal likelihood estimation using thermodynamic integration or stepping
    stone sampling.
-   Various changes of function names (e.g., metropologis -> mhg).
-   Updated examples.
-   Proper but minimal logging framework.
-   Various other changes.


## 0.4.0.0

-   Greatly improve documentation.
-   Major design change: Introduction of the `Algorithm` type class; algorithms
    are data types. See `MHG`.
-   Metropolic-coupled Markov chain Monte Carlo algorithm (parallel chains).
-   Optimal acceptance rate depends on dimension of proposal.
-   Use a circular trace with constant memory usage (big change).
-   Therefore, batch monitors use vectors now.
-   Always save chain with complete trace (but with sensible length).
-   Determine necessary trace length at initialization.
-   Rename `Item` to `Link`.
-   Rename `Status` to `Chain` and separate `Settings` and `Environment` from the
    `Chain`.
-   Many bug fixes.


## 0.3.0

-   New shorter example/test for dating trees.
-   `noData` allows running a chain without likelihood function.
-   Give proposal parameters `PName`, `PDescription`, and `PWeight` newtype
    wrappers.
-   Give `Tune` a data type.
-   Allow periodical cleansing of state (`Cleaner`).
-   Add description string to proposals, so that they can be identified in an
    easier way.
-   Add simplices and proposals on simplices.
-   `slideUniform` renamed to `slideUniformSymmetric`.
-   Merge tools into internal.
-   Do not export internal modules.


## 0.2.4

-   **Change order of arguments for proposals**.
-   &rsquo;slideStem&rsquo; was renamed to &rsquo;slideBranch&rsquo;.
-   Change ProposalSimple from newtype to type.
-   Contravariant instances of parameter and batch monitors. Use `(>$<)` instead
    of `(@.)` and `(@#)`.
-   Add `gammaDirichlet` prior for partitioned dating analyses.


## 0.2.3

-   Contrary proposals.
-   Change how monitors are lifted (use normal function, not a lens).
-   Priors.
-   Remove concurrent monitors (was slow).
-   Improve MCMC sampler output.


## 0.2.2

-   Move away from hpack.


## 0.2.1

-   Consistently use ByteString instead of Text.
-   Verbosity levels.
-   Improved handling of proposals, moves, and monitors.
-   Bactrian moves.
-   Many small changes.


## 0.1.3

Many changes; notably it is now possible to continue a Markov chain run.

