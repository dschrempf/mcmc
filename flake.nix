{
  description = "Development environment for Mcmc";

  inputs.circular.url = "github:dschrempf/circular";
  inputs.circular.inputs.nixpkgs.follows = "nixpkgs";

  inputs.covariance.url = "github:dschrempf/covariance";
  inputs.covariance.inputs.nixpkgs.follows = "nixpkgs";

  inputs.dirichlet.url = "github:dschrempf/dirichlet";
  inputs.dirichlet.inputs.nixpkgs.follows = "nixpkgs";

  inputs.flake-utils.url = "github:numtide/flake-utils";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable-small";
  # inputs.nixpkgs.url = "path:/home/dominik/Nix/Nixpkgs";

  inputs.pava.url = "github:dschrempf/pava";
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
    let
      theseHpkgNames = [
        "mcmc"
        "mcmc-examples"
        "mcmc-statistics"
      ];
      thisGhcVersion = "ghc943";
      hMkPackage = h: n: h.callCabal2nix n (./. + "/${n}") { };
      hOverlay = selfn: supern: {
        haskell = supern.haskell // {
          packageOverrides = selfh: superh:
            supern.haskell.packageOverrides selfh superh //
              nixpkgs.lib.genAttrs theseHpkgNames (hMkPackage selfh);
        };
      };
      perSystem = system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [
              hOverlay
              circular.overlays.default
              covariance.overlays.default
              dirichlet.overlays.default
              pava.overlays.default
            ];
          };
          hpkgs = pkgs.haskell.packages.${thisGhcVersion};
          hlib = pkgs.haskell.lib;
          theseHpkgs = nixpkgs.lib.genAttrs theseHpkgNames (n: hpkgs.${n});
          theseHpkgsDev = builtins.mapAttrs (_: x: hlib.doBenchmark x) theseHpkgs;
        in
        {
          packages = theseHpkgs // { default = theseHpkgs.mcmc; };

          devShells.default = hpkgs.shellFor {
            packages = _: (builtins.attrValues theseHpkgsDev);
            nativeBuildInputs = with pkgs; [
              # See https://github.com/NixOS/nixpkgs/issues/59209.
              bashInteractive

              # Haskell toolchain.
              hpkgs.cabal-fmt
              hpkgs.cabal-install
              hpkgs.haskell-language-server
            ];
            buildInputs = with pkgs; [
            ];
            doBenchmark = true;
            # withHoogle = true;
          };
        };
    in
    { overlays.default = hOverlay; } // flake-utils.lib.eachDefaultSystem perSystem;
}
