{
  description = "Development environment for Mcmc";

  inputs.circular.url = "github:dschrempf/circular";
  inputs.circular.inputs.flake-utils.follows = "flake-utils";
  inputs.circular.inputs.nixpkgs.follows = "nixpkgs";

  inputs.covariance.url = "github:dschrempf/covariance";
  inputs.covariance.inputs.flake-utils.follows = "flake-utils";
  inputs.covariance.inputs.nixpkgs.follows = "nixpkgs";

  inputs.dirichlet.url = "github:dschrempf/dirichlet";
  inputs.dirichlet.inputs.flake-utils.follows = "flake-utils";
  inputs.dirichlet.inputs.nixpkgs.follows = "nixpkgs";

  inputs.flake-utils.url = "github:numtide/flake-utils";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  inputs.pava.url = "github:dschrempf/pava";
  inputs.pava.inputs.flake-utils.follows = "flake-utils";
  inputs.pava.inputs.nixpkgs.follows = "nixpkgs";

  outputs =
    { self
    , circular
    , covariance
    , dirichlet
    , flake-utils
    , nixpkgs
    , pava
    }:
      flake-utils.lib.eachDefaultSystem (
        system:
          let
            lib = nixpkgs.lib;
            packageNames = [
              "mcmc"
              "mcmc-examples"
              "mcmc-statistics"
            ];
            mcmc-create-package = f: name: f name (./. + "/${name}") rec {};
            mcmc-overlay = (
              selfn: supern: {
                haskellPackages = supern.haskellPackages.override {
                  overrides = selfh: superh:
                    {
                      circular = circular.defaultPackage.${system};
                      covariance = covariance.defaultPackage.${system};
                      dirichlet = dirichlet.defaultPackage.${system};
                      pava = pava.defaultPackage.${system};
                    } // lib.genAttrs packageNames
                      (mcmc-create-package selfh.callCabal2nix);
                };
              }
            );
            overlays = [ mcmc-overlay ];
            pkgs = import nixpkgs {
              inherit system overlays;
            };
            # Set with packages.
            mcmc = lib.genAttrs packageNames (n: pkgs.haskellPackages.${n});
            # List with packages with benchmark dependencies for development
            # environment.
            mcmc-dev = builtins.mapAttrs (_: x: pkgs.haskell.lib.doBenchmark x) mcmc;
          in
            {
              packages = mcmc;

              defaultPackage = mcmc.mcmc;

              devShell = pkgs.haskellPackages.shellFor {
                packages = _: (builtins.attrValues mcmc-dev);
                buildInputs = with pkgs; [
                  bashInteractive
                  haskellPackages.cabal-install
                  haskellPackages.haskell-language-server
                  haskellPackages.stack
                ];
                doBenchmark = true;
                withHoogle = true;
              };
            }
      );
}
