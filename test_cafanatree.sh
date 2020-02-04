#!/bin/bash

#if using the skimmer approach
bin/mcp_skimmer --infile anatree_100.root --outfile skimmed.root --debug 1
echo "done with running skimming module"
bin/cafanatree_module --infile skimmed.root --outfile caf.root
echo "done with running cafanatree module"

#if attemting a larger-scale sample production
bin/cafanatree_module --infile anatree.root --outfile caf.root
