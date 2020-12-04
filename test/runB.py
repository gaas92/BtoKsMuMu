import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootupler")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')



process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        #MiniAOD UltraLegacy2018 we need parked data
        '/store/data/Run2018A/Charmonium/MINIAOD/12Nov2019_UL2018_rsb-v1/10000/08F41CB9-8F1F-D44F-A5FC-D17E38328C4C.root',

 )
)

process.TFileService = cms.Service("TFileService",
       fileName = cms.string('testParcked-MiniAOD.root'),
)

#### B Fitter-Tupler 
process.Btupler = cms.EDAnalyzer('BtoKsMuMu', 
                          muons = cms.InputTag("slimmedMuons"),
                          #packedPFcand = cms.InputTag("packedPFCandidates"),
                          #primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          #secundaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                          #bslabel = cms.InputTag("offlineBeamSpot"),
                          #TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          #OnlyBest = cms.bool(False),
                          #isMC = cms.bool(False),
                          #OnlyGen = cms.bool(False)                               
    )

process.p = cms.Path(process.Btupler)
