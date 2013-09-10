#!/bin/sh
source /vols/cms02/bainbrid/qcd/trigger/SUSY2/setup.sh
cd /vols/cms02/bainbrid/qcd/trigger/SUSY2/bryn/python
./Data.py -j ./Data.py_20110421_10_37_30.json -J ${SGE_TASK_ID}