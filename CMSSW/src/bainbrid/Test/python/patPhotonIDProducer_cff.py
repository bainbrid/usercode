import FWCore.ParameterSet.Config as cms

from RecoEgamma.PhotonIdentification.photonId_cfi import PhotonIDProd

patPhotonIDProducer = cms.EDProducer(
    "PATPhotonIDProducer",
    Photons  = cms.InputTag("selectedLayer1Photons")
    )

for i in dir(PhotonIDProd):
    if isinstance(getattr(PhotonIDProd,i),cms._ParameterTypeBase):
        setattr(patPhotonIDProducer,i,getattr(PhotonIDProd,i))
