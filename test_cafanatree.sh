#!/bin/bash

echo "done with running cafanatree module"

#if attemting a larger-scale sample production
bin/cafanatree_module --infile anatree.root --outfile caf.root
