import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")
from FWCore.ParameterSet.VarParsing import VarParsing

'''
############ Command line args ################
'''

options = VarParsing('python')

options.register('isRes', False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Run Resonant or non resonant MC")

options.register('onlyGen', False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Run MC onlyGen level")

options.register('globalTag', '102X_upgrade2018_realistic_v12',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "Set global tag"
) 

options.register('maxE', -1,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.int,
                "Maximum number of events"
)

options.register('tg', 'PrivateMC_Condor_Result',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "tag for outputfile"
)

options.register('inputFile', 'noProbeFilterDecayFilter_MiniAOD1_0.txt',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "input txt file"
)

options.register('saveInSync', False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "save result file in sync")

options.parseArguments()


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v24', '')# for 2018
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v20', '')# for 2017
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')# for 2016

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxE))
files_my_gen = [string for string in open('myGenFiles/'+options.inputFile).readlines() if len(string) > 10]
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(files_my_gen
 )
)


process.load("myAnalyzers.BtoKsMuMu.Psiks0_BestPA_V0Ext_Rootupler_cfi")
process.rootuple.isMC = cms.bool(True) # this is only for test
process.rootuple.isRes = cms.bool(options.isRes)
process.rootuple.OnlyGen = cms.bool(options.onlyGen)
process.rootuple.GenParticles = cms.InputTag("prunedGenParticles") 


outname = 'PrivateMC_{}_{}_part{}.root'.format('Res' if options.isRes  else 'notRes', 'Gen' if options.onlyGen else 'Reco' , str(options.inputFile)[len(options.inputFile)-7: len(options.inputFile)-4])


process.TFileService = cms.Service("TFileService",

       fileName = cms.string(outname),
)

process.p = cms.Path(process.rootuple)
 
