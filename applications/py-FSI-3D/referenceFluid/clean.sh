#!/bin/sh
set -e -u

set -e -u
echo "--- Cleaning up OpenFOAM case in $(pwd)"
if [ -n "${WM_PROJECT:-}" ] || error "No OpenFOAM environment is active."; then
    # shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
    . "${WM_PROJECT_DIR}/bin/tools/CleanFunctions"
    cleanCase
    rm -rfv 0/uniform/functionObjects/functionObjectProperties
fi

