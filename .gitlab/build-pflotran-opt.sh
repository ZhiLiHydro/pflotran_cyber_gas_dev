#!/bin/sh

. $PFLOTRAN_DIR/.gitlab/skip_on_error.sh

cd $SRC_DIR

make clean
make -j4 pflotran_rxn
# prevent building of pflotran if pflotran_rxn is not built correctly
if [ ! -f pflotran_rxn ]; then
  echo "\n----- pflotran_rxn executable not properly compiled -----\n"
  echo 'failed' > $ARTIFACT_DIR/status
else
  make clean
  make -j4 gnu_code_coverage=1 gnu_runtime_checks=1 pflotran
fi
