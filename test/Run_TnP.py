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


isMC = True

if isMC:
    globaltag_ = '102X_upgrade2018_realistic_v12'
    tagNprobe_ = 'TagAndProbeProducer_MC' 
    testData   = '/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/08A85CD9-6209-9A41-AFDC-648C3688EA3D.root'
else :
    globaltag_ = '102X_dataRun2_v11'
    tagNprobe_ = 'TagAndProbeProducer' 
    testData   = '/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root' 

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag_, '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

'''
############ Command line args ################
'''

#args = VarParsing.VarParsing('analysis')
#args.register('inputFile', '', args.multiplicity.list, args.varType.string, "Input file or template for glob")
#args.outputFile = ''
#args.parseArguments()

'''
#####################   Input    ###################
'''
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(args.maxEvents)
    input = cms.untracked.int32(-1)
)

#from glob import glob
#if args.inputFile:
#    if len(args.inputFile) == 1 and '*' in args.inputFile[0]:
#        flist = glob(args.inputFile[0])
#    else:
#        flist = args.inputFile
#elif args.inputFiles:
#    if len(args.inputFiles) == 1 and args.inputFiles[0].endswith('.txt'):
#        with open(args.inputFiles[0]) as f:
#            flist = [l[:-1] for l in f.readlines()]
#    else:
#        flist = args.inputFiles
#else:
#    fdefault = os.environ['CMSSW_BASE'] + '/src/ntuplizer/BPH_RDntuplizer/production/'
#    # fdefault += 'inputFiles_BP_Tag_B0_MuNuDmst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3.txt'
#    fdefault += 'inputFiles_BP_Tag-Probe_B0_JpsiKst_Hardbbbar_evtgen_HELAMP_PUc0_10-2-3.txt'
#    with open(fdefault) as f:
#        flist = [l[:-1] for l in f.readlines()]
#    flist = flist[:10]
#
#for i in range(len(flist)):
#    if os.path.isfile(flist[i]):
#        flist[i] = 'file:' + flist[i]

#process.source = cms.Source("PoolSource",
#                            #fileNames = cms.untracked.vstring(tuple(flist)),
#                            fileNames = cms.untracked.vstring(
#                                '/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/08A85CD9-6209-9A41-AFDC-648C3688EA3D.root' 
#                            ),
#                            inputCommands=cms.untracked.vstring('keep *',
#                                                                'drop GenLumiInfoHeader_generator__SIM'),
#                            skipBadFiles=cms.untracked.bool(True)
#                           )
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        testData
        #MC-RES /BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM
        #'/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/08A85CD9-6209-9A41-AFDC-648C3688EA3D.root',
        #'/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/260000/51C3402C-206E-FA4E-9BCF-C1CD7DD79786.root' 
        #'/store/mc/RunIIAutumn18MiniAOD/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/N1_102X_upgrade2018_realistic_v15-v1/00000/0923C65B-A0B1-6F4B-A7FB-54A2ECCFE4B3.root'
        #Parked Data
        #'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root'

        #Files with problems in crab:
        #'/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/50004/44EDA4E4-0807-8D4C-9ED5-6159DA0D5C2C.root',
        #'/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/50004/45F5AB3C-34F0-5949-A7AA-1B21643731DA.root'
 ),
    inputCommands=cms.untracked.vstring('keep *',
                                        'drop GenLumiInfoHeader_generator__SIM'),
    skipBadFiles=cms.untracked.bool(True)
)

'''
#####################   Output   ###################
'''
#if args.outputFile == '.root':
#    outname = 'TagAndProbeTrigger_CAND.root'
#elif args.outputFile.startswith('_numEvent'):
#    outname = 'TagAndProbeTrigger_CAND' + args.outputFile
#else:
#    outname = args.outputFile

outname = 'TagAndProbeTrigger_CAND.root'

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
process.TnP = cms.EDFilter(tagNprobe_,
        muonIDScaleFactors = cms.int32(0),
        requireTag = cms.int32(1),
        verbose = cms.int32(1)
)


process.p = cms.Path(
                    process.l1bits +
                    process.TnP
                    )


# DEBUG -- dump the event content
# process.output = cms.OutputModule(
#                 "PoolOutputModule",
#                       fileName = cms.untracked.string('edm_output.root'),
#                       )
# process.output_step = cms.EndPath(process.output)
#
# process.schedule = cms.Schedule(
# 		process.p,
# 		process.output_step)


'''
#############   Overall settings    ################
'''

process.MessageLogger.cerr.FwkReport.reportEvery = 10000