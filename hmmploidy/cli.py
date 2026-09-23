"""Unified command-line entry point:

    hmmploidy simulate ...
    hmmploidy genolikes ...
    hmmploidy run ...
"""
from __future__ import annotations

import sys


def main(argv=None) -> None:
    argv = list(argv) if argv is not None else sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return

    command, rest = argv[0], argv[1:]
    if command == "simulate":
        from . import simulate

        simulate.main(rest)
    elif command == "genolikes":
        from . import genotype_likelihoods

        genotype_likelihoods.main(rest)
    elif command == "run":
        from . import run

        run.main(rest)
    else:
        raise SystemExit(f"Unknown command '{command}'. Use simulate | genolikes | run")


if __name__ == "__main__":
    main()
