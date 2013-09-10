#!/usr/bin/env python

import setupSUSY
from icf.core import PSet
from libbainbrid import *

def createRatioPlots( ratio,
                      truth = False,
                      pt_default = 50.,
                      ht_default = 350.,
                      pt1_default = 125.,
                      pt2_default = 100. ) :

      # Define pT bins to be used
      pt_min = 20
      pt_max = 50
      pt_step = 10
      
      # Define x3 from pT and HT defaults for analysis 
      x3_default = (2.*pt_default)/(ht_default+pt_default)
      factor = (2.-x3_default)/x3_default

      # Determine maximum values for HT and Meff (from max pT value)
      ht_max = pt_max*1. * factor
      meff_max = pt_max*1. + ht_max

      # Iterate through AlphaT values
      for ii in range(0,11) :
            alphat = 500 + ii * 5

            # Iterate through pT values
            for pt in range(pt_min,pt_max+pt_step,pt_step) :

                  # Define HT and Meff from pT value
                  ht = pt*1. * factor
                  meff = pt*1. + ht

                  # Define scaled pT values appropriate for given Meff
                  pt1 = pt1_default * meff / meff_max
                  pt2 = pt2_default * meff / meff_max

                  # Iterate through scaled pT values
                  for jj in range(0.8,1.2,0.1) :

                        # Define scaled pT values appropriate for given Meff
                        pt1 *= jj
                        pt2 *= jj
                        
                        # Define histograms 
                        dir = "Ratio"+str(alphat)+"Pt"+str(pt)+"Scale"+str(jj)
                        ratio.append( RobPlottingOps( PSet(DirName = dir,
                                                           MinObjects=2,
                                                           MaxObjects=8,
                                                           Ratio = True,
                                                           MinJetPt = pt*1.,
                                                           MinJetPt1 = pt1,
                                                           MinJetPt2 = pt2,
                                                           UseGen = truth,
                                                           MinGenPt = pt_min*1.,
                                                           MaxGenMET = 10.,
                                                           AlphaTcut = (alphat/1000.),
                                                           ).ps() ) )

def attachRatioPlots( ratio,
                      tree,
                      last_operation ) :
      
      if len(ratio) > 0 :
            for kk in range(0,len(ratio)) :
                  if kk == 0 :
                        tree.TAttach(last_operation,ratio[kk])
                  else :
                        tree.TAttach(ratio[kk-1],ratio[kk])

      
