import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

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
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
files_my_gen = [string for string in open('noProbeFilterDecayFilter_MiniAOD.txt').readlines() if len(string) > 10]
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(files_my_gen

        #Parked Data
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2129F080-982D-A649-8D20-4944620E99A6.root',
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/12C659FD-D961-6640-9C3D-48CD12A0D033.root', 
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root'

        #MiniAOD UltraLegacy2018
        #'/store/data/Run2018A/Charmonium/MINIAOD/12Nov2019_UL2018_rsb-v1/10000/08F41CB9-8F1F-D44F-A5FC-D17E38328C4C.root',

        #MiniAOD UltraLegacy2017
        #'/store/data/Run2017F/Charmonium/MINIAOD/09Aug2019_UL2017-v1/20000/00BACB48-9B0F-8F48-A68B-2F08A3E9E681.root',

        #MiniAOD UltraLegacy2016
        #'/store/data/Run2016G/Charmonium/MINIAOD/21Feb2020_UL2016-v1/30000/00013A18-278D-5B48-9BEF-1083A8F5C9D7.root',
        
 )
)

#process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
#                                        triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
#                                                                        'HLT_Dimuon25_Jpsi_v*',
#                                                                        'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
#                                                                        'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
#                                                                        'HLT_DoubleMu4_Jpsi_Displaced_v*'                                   
#                                                                       ),
#                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
#                                        l1tResults = cms.InputTag( "" ),
#                                        throw = cms.bool(False)
#                                        )

process.load("myAnalyzers.BtoKsMuMu.Psiks0_BestPA_V0Ext_Rootupler_cfi")
process.rootuple.isMC = cms.bool(True) # this is only for test
process.rootuple.isRes = cms.bool(False)
process.rootuple.OnlyGen = cms.bool(True)
process.rootuple.GenParticles = cms.InputTag("prunedGenParticles") 

process.TFileService = cms.Service("TFileService",

       fileName = cms.string('Rootuple_BdtoMuMuks0_PARKED_Bpa_PVExt_MiniAODGen.root'),
)

#process.mySequence = cms.Sequence(
#                                   process.triggerSelection *
#                                   process.rootuple
#				   )
#
#process.p = cms.Path(process.mySequence)
process.p = cms.Path(process.rootuple)

