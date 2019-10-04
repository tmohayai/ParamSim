#/bin/bash

cd ${PWD}/build
# rm -rf *
# cmake ..
# echo "done with cmake"
make install
echo "done with make install"
cd ..
bin/cafanatree_module --edepfile anatree_100.root --outfile 2.root
echo "done with re-running cafanatree module"
