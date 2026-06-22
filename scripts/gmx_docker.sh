#!/usr/bin/env bash
set -euo pipefail

IMAGE="${GMX_IMAGE:-gromacs/gromacs:latest}"

# This wrapper is called like:
# scripts/gmx_docker.sh pdb2gmx -f target.pdb ...
#
# It mounts the current working directory as /work.
# gromacs.py copies minim.mdp into the working directory, so relative paths work.
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  "$IMAGE" \
  gmx "$@"
