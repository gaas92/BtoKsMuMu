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
#files_my_gen = [string for string in open('myGenFiles/noProbeFilterDecayFilter_MiniAOD2_4.txt').readlines() if len(string) > 10]
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#files_my_gen

        #Parked Data
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/2129F080-982D-A649-8D20-4944620E99A6.root',
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/12C659FD-D961-6640-9C3D-48CD12A0D033.root', 
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root'
        #Failed Datafiles
        '/store/data/Run2018B/ParkingBPH1/MINIAOD/05May2019-v2/230000/BDE201E3-92B3-0640-A9B0-9FA6D4A64DB2.root', 
        '/store/data/Run2018B/ParkingBPH1/MINIAOD/05May2019-v2/230000/BEF7C437-2D8F-9E47-9387-2CBF82813237.root', 

        '/store/data/Run2018B/ParkingBPH1/MINIAOD/05May2019-v2/230000/DF476FF5-ABC1-8C42-85EB-CA30E74209F7.root', 

        '/store/data/Run2018D/ParkingBPH2/MINIAOD/05May2019promptD-v1/130007/97C9FFF0-E03C-8A48-A8C3-F06CD63E65EC.root', 
        '/store/data/Run2018D/ParkingBPH2/MINIAOD/05May2019promptD-v1/130007/994F7213-2598-4843-BB41-CE653908B894.root',
        '/store/data/Run2018D/ParkingBPH2/MINIAOD/05May2019promptD-v1/130007/99910C92-D70D-EC45-B108-AE1D6BD2096B.root',

        '/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/50004/44C185C3-CF69-3E41-8877-3F70821771DB.root', 
        '/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/50004/44EDA4E4-0807-8D4C-9ED5-6159DA0D5C2C.root', 
        '/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/50004/456859F1-50ED-C34F-8C20-EA784E3729A7.root', 
        '/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/50004/456859F1-50ED-C34F-8C20-EA784E3729A7.root' 

        # Test on MC 
        #'/store/mc/RunIIAutumn18MiniAOD/BdToK0sMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/17621299-2380-064D-A5E9-B78072D95A7D.root'   

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
process.rootuple.isMC = cms.bool(False) # this is only for test
process.rootuple.isRes = cms.bool(False)
process.rootuple.OnlyGen = cms.bool(False)
process.rootuple.GenParticles = cms.InputTag("prunedGenParticles") 

process.TFileService = cms.Service("TFileService",

       fileName = cms.string('Rootuple_BdtoMuMuks0_PARKED_Bpa_PVExt_MiniAOD.root'),
       #fileName = cms.string('Rootuple_BdtoMuMuks0_PARKED_Bpa_PVExt_MiniAOD_GEN_24.root'),
)

#process.mySequence = cms.Sequence(
#                                   process.triggerSelection *
#                                   process.rootuple
#				   )
#
#process.p = cms.Path(process.mySequence)
process.p = cms.Path(process.rootuple)
 
