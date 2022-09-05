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
        ghcVersion = "ghc924";
        haskellMkPackage = f: name: f name (./. + "/${name}") rec { };
        haskellOverlay = (
          selfn: supern: {
            haskellPackages = supern.haskell.packages.${ghcVersion}.override {
              overrides = selfh: superh:
                {
                  circular = circular.packages.${system}.default;
                  covariance = covariance.packages.${system}.default;
                  dirichlet = dirichlet.packages.${system}.default;
                  pava = pava.packages.${system}.default;
                } // lib.genAttrs packageNames (haskellMkPackage selfh.callCabal2nix);
            };
          }
        );
        overlays = [ haskellOverlay ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        hpkgs = pkgs.haskellPackages;
        # Set with packages.
        mcmcPkgs = lib.genAttrs packageNames (n: hpkgs.${n});
        # List with packages with benchmark dependencies for development
        # environment.
        mcmcPkgsDev = builtins.mapAttrs (_: x: pkgs.haskell.lib.doBenchmark x) mcmcPkgs;
      in
      {
        packages = mcmcPkgs // { default = mcmcPkgs.mcmc; };

        devShells.default = hpkgs.shellFor {
          packages = _: (builtins.attrValues mcmcDev);
          buildInputs = with pkgs; [
            bashInteractive

            # TODO: `cabal-fmt` fails to build when using a newer package set.
            haskell.packages.ghc902.cabal-fmt

            hpkgs.cabal-install
            hpkgs.haskell-language-server
          ];
          doBenchmark = true;
          # withHoogle = true;
        };
      }
    );
}
