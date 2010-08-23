#!/usr/bin/env python


# To use, in framework/src/common/EventElement.cc:
# Change "print_branches_to_keep" to true...

# -----------------------------------------------------------------------------
# Includes

from golden_cff import *
from samples_cff import *

# -----------------------------------------------------------------------------
# Configuration

# Available configurations
conf_gic = deepcopy(conf_ic5_gen)
conf_gak = deepcopy(conf_ak5_gen)
conf_ic5 = deepcopy(conf_ic5_calo)
conf_ak5 = deepcopy(conf_ak5_calo)
conf_jpt = deepcopy(conf_ak5_jpt)
conf_pf  = deepcopy(conf_ak5_pf)

# Chosen configuration
conf = conf_ak5

# -----------------------------------------------------------------------------
# Analyses

dataset = qcd6[0]

dataset.FirstEntry = 1
dataset.LastEntry  = 2

a=Analysis("EventContent")
a.Run("../results",conf,[dataset])

