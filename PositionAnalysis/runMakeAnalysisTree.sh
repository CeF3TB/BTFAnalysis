#!/bin/bash                                                                                                                                                  

INPUT_DIR=/t3/users/cmsdaqtb/data/data/BTF/CeF3/runs/dataTree
TAG=$1

#ln -s $INPUT_DIR data_2

for file in $(ls $INPUT_DIR/run_BTF_*| awk '{FS="n_"}{print $2}'|awk '{FS=".r"}{print $1}'); do
     ./makeAnalysisTree $file $TAG
done