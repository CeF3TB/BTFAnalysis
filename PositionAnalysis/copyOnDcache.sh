#!/bin/bash

VERSION=$1

mkdir /pnfs/roma1.infn.it/data/cms/store/user/cmsdaqtb/analysisTrees/analysisTrees_$VERSION

for i in /t3/users/cmsdaqtb/data/data/BTF/CeF3/runs/analysisTrees_$VERSION/*root; do
    dccp $i /pnfs/roma1.infn.it/data/cms/store/user/cmsdaqtb/analysisTrees/analysisTrees_$VERSION
done
