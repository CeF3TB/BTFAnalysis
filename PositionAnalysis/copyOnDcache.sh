#!/bin/bash

VERSION=$1

for i in /t3/users/cmsdaqtb/data/data/BTF/CeF3/runs/analysisTrees_$VERSION/*root; do
    echo "dccp $i /pnfs/roma1.infn.it/data/cms/store/user/cmsdaqtb/analysisTrees/analysisTrees_V00"
done
