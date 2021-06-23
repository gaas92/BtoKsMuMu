#Taken from: https://github.com/ocerri/BPH_RDntuplizer/blob/master/config/cmssw_privateMC_TagAndProbeTrigger.py
import os, sys
import FWCore.ParameterSet.Config as cms
#import FWCore.ParameterSet.VarParsing as VarParsing
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BPHRDntuplizer', eras.Run2_2018)
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')


# Needed for transient track builder
# process.load('Configuration.StandardSequences.Services_cff')
# process.load('Configuration.EventContent.EventContent_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


'''
############ Command line args ################
'''

options = VarParsing('python')

options.register('isMC', False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Run MC or RD")

options.register('onlyGen', False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Run MC onlyGen level")

options.register('globalTag', '102X_dataRun2_v11',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "Set global tag"
) 

options.register('maxE', -1,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.int,
                "Maximum number of events"
)

options.register('tg', 'TnP_Condor_Result',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "tag for outputfile"
)

options.register('inputFile', 'noProbeFilterDecayFilter_MiniAOD1_0.txt',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "input txt file"
)

options.parseArguments()

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v12', '')
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag , '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

'''
#####################   Input    ###################
'''
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxE)
)

files_my_gen = [string for string in open('myGenFiles/'+options.inputFile).readlines() if len(string) > 10]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #Parked Data
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root'
        files_my_gen

 ),
    inputCommands=cms.untracked.vstring('keep *',
                                        'drop GenLumiInfoHeader_generator__SIM'),
    skipBadFiles=cms.untracked.bool(True)
)

'''
#####################   Output   ###################
'''


#outname = f'TagAndProbeTrigger_{ "MC" if options.isMC esle "RD"}_part{options.input_file[len(options.input_file)-7: len(options.input_file)-4]}.root'
#outname = 'TagAndProbeTrigger_{}_part{}.root'.format( "MC" if options.isMC esle "RD", options.input_file[len(options.input_file)-7: len(options.input_file)-4])
outname = 'TagAndProbeTrigger_part.root'

process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outname),
      #closeFileFast = cms.untracked.bool(True)
      )



'''
#################   Sequence    ####################
'''

process.l1bits=cms.EDProducer("L1TriggerResultsConverter",
                              src=cms.InputTag("gtStage2Digis"),
                              # src_ext=cms.InputTag("gtStage2Digis"),
                              # storeUnprefireableBit=cms.bool(True),
                              legacyL1=cms.bool(False),
)

#process.TnP = cms.EDFilter("TagAndProbeProducer_MC",
process.TnP = cms.EDFilter("TagAndProbeProducer",
        muonIDScaleFactors = cms.int32(0),
        requireTag = cms.int32(1),
        verbose = cms.int32(0)
)


process.p = cms.Path(
                    process.l1bits +
                    process.TnP
                    )



'''
#############   Overall settings    ################
'''

process.MessageLogger.cerr.FwkReport.reportEvery = 10000