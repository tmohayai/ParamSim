#/bin/bash

if [ ! -d ${PWD}/build ]; then
  mkdir -p ${PWD}/build
fi

cd ${PWD}/build
rm -rf *
cmake ..
echo "done with cmake"
make -j4 install
echo "done with make install"
cd ..
bin/cafanatree_module --infile anatree_100.root --outfile 2.root
echo "done with re-running cafanatree module"
