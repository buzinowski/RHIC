#!/bin/bash

echo "Processing cbmsim_reduce.C..."
root.exe -l <<EOF
.L cbmsim_reduced.C
cbmsim_reduced t("mpddst_reduced_9971543_998.root")
t.Loop("foutput.root")
EOF
root.exe -l foutput.root 

