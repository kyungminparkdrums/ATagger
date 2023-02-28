import FWCore.ParameterSet.Config as cms

calibratedPhotonProducer = cms.EDProducer('CalibratedPhotonProducer',
  src = cms.InputTag('gedGsfElectrons'),
  recHitCollectionEB = cms.InputTag('reducedEcalRecHitsEB'),
  recHitCollectionEE = cms.InputTag('reducedEcalRecHitsEE'),
  correctionFile = cms.string(''),
  minEtToCalibrate = cms.double(5),
  produceCalibratedObjs = cms.bool(True),
  semiDeterministic = cms.bool(True),
  valueMapsStored = cms.vstring(
    'energyScaleStatUp',
    'energyScaleStatDown',
    'energyScaleSystUp',
    'energyScaleSystDown',
    'energyScaleGainUp',
    'energyScaleGainDown',
    'energySigmaRhoUp',
    'energySigmaRhoDown',
    'energySigmaPhiUp',
    'energySigmaPhiDown',
    'energyScaleUp',
    'energyScaleDown',
    'energySigmaUp',
    'energySigmaDown',
    'energyScaleValue',
    'energySigmaValue',
    'energySmearNrSigma',
    'ecalEnergyPreCorr',
    'ecalEnergyErrPreCorr',
    'ecalEnergyPostCorr',
    'ecalEnergyErrPostCorr'
  )
)
