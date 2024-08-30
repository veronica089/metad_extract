#!/bin/bash

#you get filenotfound error not recognized without activating conda
#for sys in "wt_Glu0" "mut_Glu0" "wt_Glu-" "mut_Glu-"
for sys in "mut_Glu0" "wt_Glu0" "mut_Glu-" "wt_Glu-"
do
	for struct in 0 1
	do
        python get_dist.py $sys $struct 
        done
done
