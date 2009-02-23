import FWCore.ParameterSet.Config as cms

# $Id: j2tParametersVX_cfi.py,v 1.3 2009/02/23 13:12:13 bainbrid Exp $
j2tParametersVX = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    trackQuality = cms.string("goodIterative"),
    coneSize = cms.double(0.5)
)

