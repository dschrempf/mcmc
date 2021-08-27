{
  description = "Development environment for Mcmc";

  inputs.flake-utils.url = "github:numtide/flake-utils";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, flake-utils, nixpkgs }:
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
                  lib.genAttrs packageNames
                    (mcmc-create-package selfh.callCabal2nix);
              };
            }
          );
          overlays = [ mcmc-overlay ];
          pkgs = import nixpkgs { inherit system overlays; };
          # Set with packages.
          mcmc = lib.genAttrs packageNames (n: pkgs.haskellPackages.${n});
          # List with packages with benchmark dependencies for development
          # environment.
          mcmc-dev = builtins.mapAttrs (_: x: pkgs.haskell.lib.doBenchmark x) mcmc;
        in
          {
            packages = mcmc;

            devShell = pkgs.haskellPackages.shellFor {
              packages = _: (builtins.attrValues mcmc-dev);
              buildInputs = with pkgs.haskellPackages; [
                haskell-language-server
                cabal-install
              ];
              doBenchmark = true;
            };
          }
    );
}
