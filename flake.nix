{
  description = "Development environment for Mcmc.";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (
      system:
        let
          pkgs = import nixpkgs { inherit system; };
        in
          {
            devShell = pkgs.mkShell {
              nativeBuildInputs = with pkgs; [ cabal-install ghc ];
              buildInputs = with pkgs; [
                zlib
              ];
            };
          }
    );
}
