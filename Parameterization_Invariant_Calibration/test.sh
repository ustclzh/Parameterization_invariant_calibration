#!/bin/bash
/opt/matlab/latest/bin/matlab -nodesktop -nosplash -r "inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7; design_criterion = $8; run('simulation_cluster.m');quit"