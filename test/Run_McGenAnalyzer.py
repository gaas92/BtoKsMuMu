import FWCore.ParameterSet.Config as cms
from glob import glob

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring( open('my_gen-NF_files.txt').readlines())
		)



process.TFileService = cms.Service("TFileService",
         fileName = cms.string('test_B0Gen_NoFilter.root'),                                  
)

#Analyzer
process.demo = cms.EDAnalyzer('McGenAnalyzer',
                              genParticles = cms.InputTag('genParticles'))

process.p = cms.Path(process.demo)