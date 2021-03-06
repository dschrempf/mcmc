
# Markov chain Monte Carlo sampling - ChangeLog


## Unreleased changes

-   Improve documentation.
-   Generalized priors allowing automatic differentiation.
-   Hamiltonian proposal (will be improved).


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
-   'slideStem' was renamed to 'slideBranch'.
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

