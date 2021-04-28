import FWCore.ParameterSet.Config as cms
from glob import glob

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring( open('myGenFiles/Res_noProbeFilterDecayFilter_Gen.txt').readlines())
        #fileNames = cms.untracked.vstring(
        #    ['/store/mc/RunIIAutumn18DR/BdToK0sMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/GEN-SIM-DIGI-RAW/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/100000/001BF693-83CF-4648-AC6B-43B3EB4EAE30.root']
        #)
		)



process.TFileService = cms.Service("TFileService",
         fileName = cms.string('test_B0ResGen_noProbeFilterDecayFilter.root'),                                  
)

#Analyzer
process.demo = cms.EDAnalyzer('McGenAnalyzer',
                              genParticles = cms.InputTag('genParticles'))

process.p = cms.Path(process.demo)