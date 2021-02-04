import FWCore.ParameterSet.Config as cms
from glob import glob

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
        #fileNames = cms.untracked.vstring( open('my_gen-NF_files.txt').readlines())
        fileNames = cms.untracked.vstring(
            ['file:/eos/user/g/gayalasa/Step0-BPHnoFilters_B0toK0MM/crab_2018-PrivateMC-Step0-BPHnoFilters_B0toK0MM-2021-01-27-19-15/210127_181518/0000/BPH_B0toK0MM_GenNF_1.root']
        )
		)



process.TFileService = cms.Service("TFileService",
         fileName = cms.string('test_B0Gen_NoFilter.root'),                                  
)

#Analyzer
process.demo = cms.EDAnalyzer('McGenAnalyzer',
                              genParticles = cms.InputTag('genParticles'))

process.p = cms.Path(process.demo)