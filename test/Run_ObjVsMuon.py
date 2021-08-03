#Taken from: https://github.com/ocerri/BPH_RDntuplizer/blob/master/config/cmssw_privateMC_TagAndProbeTrigger.py
import os, sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

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

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v12', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')



'''
#####################   Input    ###################
'''
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(args.maxEvents)
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #MC-RES /BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM
        #'/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/08A85CD9-6209-9A41-AFDC-648C3688EA3D.root',
        #'/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/260000/51C3402C-206E-FA4E-9BCF-C1CD7DD79786.root' 
        #'/store/mc/RunIIAutumn18MiniAOD/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/N1_102X_upgrade2018_realistic_v15-v1/00000/0923C65B-A0B1-6F4B-A7FB-54A2ECCFE4B3.root'
        #Parked Data
        '/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root'

 ),
    inputCommands=cms.untracked.vstring('keep *',
                                        'drop GenLumiInfoHeader_generator__SIM'),
    skipBadFiles=cms.untracked.bool(True)
)

'''
#####################   Output   ###################
'''

outname = 'ObjVsMuon_tuple.root'

process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outname),
      #closeFileFast = cms.untracked.bool(True)
      )



'''
#################   Sequence    ####################
'''

process.ObjVsMu = cms.EDProducer("ObjVsMuon_tupler",
        triggerBits = cms.InputTag("TriggerResults","","HLT"),
        objects = cms.InputTag("slimmedPatTrigger"),
        verbose = cms.int32(0)
)

process.p = cms.Path(
                    process.ObjVsMu
                    )


'''
#############   Overall settings    ################
'''

process.MessageLogger.cerr.FwkReport.reportEvery = 10000