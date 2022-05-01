import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")
from FWCore.ParameterSet.VarParsing import VarParsing


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

globaltag_ = '102X_upgrade2018_realistic_v12'
tagNprobe_ = 'TagAndProbeProducer_MC' 
testData   = '/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/08A85CD9-6209-9A41-AFDC-648C3688EA3D.root'

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag_, '')

process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/BdToK0sMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/17621299-2380-064D-A5E9-B78072D95A7D.root'
 )
)


outname = 'MC_onlyGen.root'

process.load("myAnalyzers.BtoKsMuMu.Psiks0_OnlyGen_cfi")

process.rootuple.isMC = cms.bool(True) # this is only for test
process.rootuple.isRes = cms.bool(False)
process.rootuple.OnlyGen = cms.bool(True)
process.rootuple.GenParticles = cms.InputTag("prunedGenParticles") 


process.TFileService = cms.Service("TFileService",

       fileName = cms.string(outname),
)

process.p = cms.Path(process.rootuple)
 
