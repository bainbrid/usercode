#!/usr/bin/env python

import setupSUSY
from libFrameworkSUSY import *
from libbainbrid import *
from icf.core import PSet

toy = ToyMC(
    PSet(
    OutputFile="output",
    Directory="Default",
    nAT = [ 0.55 ],
    nPT = [ 50. ],
    nHT = [ 150., 200., 250., 300., 350. ],
#    nAT = [ 0.50, 0.51, 0.52, 0.53, 0.54, 0.55 ],
#    nPT = [ 50., 40., 30., 20. ],
#    nHT = [ 150., 200., 250., 300., 350., 400., 450 ],
    nBins = [ 20000 ],
    ).ps()
    )

