{
  description = "Development environment for Mcmc";

  inputs.circular.url = "github:dschrempf/circular";

  inputs.covariance.url = "github:dschrempf/covariance";

  inputs.dirichlet.url = "github:dschrempf/dirichlet";

  inputs.flake-utils.url = "github:numtide/flake-utils";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable-small";
  # inputs.nixpkgs.url = "path:/home/dominik/Nix/Nixpkgs";

  inputs.pava.url = "github:dschrempf/pava";

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
        mcmc-create-package = f: name: f name (./. + "/${name}") rec { };
        mcmc-overlay = (
          selfn: supern: {
            haskellPackages = supern.haskell.packages.ghc923.override {
              overrides = selfh: superh:
                {
                  circular = circular.packages.${system}.default;
                  covariance = covariance.packages.${system}.default;
                  dirichlet = dirichlet.packages.${system}.default;
                  pava = pava.packages.${system}.default;
                } // lib.genAttrs packageNames
                  (mcmc-create-package selfh.callCabal2nix);
            };
          }
        );
        overlays = [ mcmc-overlay ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        hpkgs = pkgs.haskellPackages;
        # Set with packages.
        mcmc = lib.genAttrs packageNames (n: hpkgs.${n});
        # List with packages with benchmark dependencies for development
        # environment.
        mcmc-dev = builtins.mapAttrs (_: x: pkgs.haskell.lib.doBenchmark x) mcmc;
      in
      {
        packages.default = mcmc;

        devShells.default = hpkgs.shellFor {
          packages = _: (builtins.attrValues mcmc-dev);
          buildInputs = with pkgs; [
            bashInteractive
            hpkgs.cabal-install
            hpkgs.haskell-language-server
          ];
          doBenchmark = true;
          # withHoogle = true;
        };
      }
    );
}
