import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")
from FWCore.ParameterSet.VarParsing import VarParsing
 
'''
############ Command line args ################
'''
 
options = VarParsing('python') 

options.register('isMC', True,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Run MC or RD")

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
'''new'''
options.register('inputMINIAODFile', '/store/data/Run2018B/ParkingBPH3/MINIAOD/05May2019-v2/40002/3294F295-DBA4-CD40-9F7E-C8DF4B144DF1.root',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "Single MiniAOD file"
)
options.register('singleFile', False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "if Single MiniAOD file"
)
options.register('njob', 0,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.int,
                "Id of the job"
)
'''end new'''

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

if options.isMC:
    globaltag_ = '102X_upgrade2018_realistic_v12'
    tagNprobe_ = 'TagAndProbeProducer_MC' 
    testData   = '/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/08A85CD9-6209-9A41-AFDC-648C3688EA3D.root'
else :
    globaltag_ = '102X_dataRun2_v11'
    tagNprobe_ = 'TagAndProbeProducer' 
    testData   = '/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/05030414-6C93-DC46-AD09-68D76E2FB466.root' 


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag_, '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v24', '')# for 2018
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v20', '')# for 2017
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')# for 2016

process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

'''
#####################   Input    ###################
'''

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxE))
if options.singleFile:
    files_to_run = str(options.inputMINIAODFile)
else:
    files_to_run = [string for string in open('myGenFiles/'+options.inputFile).readlines() if len(string) > 10]    

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(files_to_run
    fileNames = cms.untracked.vstring(['/store/mc/RunIIAutumn18MiniAOD/BdToK0sJPsi_ToMuMu_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/100000/08A85CD9-6209-9A41-AFDC-648C3688EA3D.root']
 )
)

'''
#####################   Output   ###################
'''

if options.singleFile:
    container = str(int(options.njob/1000))
    if 'Run2018B' in files_to_run and 'ParkingBPH3' in files_to_run:
        folder = 'B{}000'.format(container)
        pattern = ''
    else :
        folder = 'failed_crab_jobs'   
        pattern = files_to_run.split('/')[4] + '_' +files_to_run.split('/')[3]

    outname = 'Data_Reco/{}/SpecialReco_{}_part{}.root'.format(folder,'MC' if options.isMC  else 'RD'+pattern,  options.njob)

else:
    outname = 'PrivateMC_{}_{}_part{}.root'.format('JPsi' if options.isRes  else 'MuMu', 'Gen' if options.onlyGen else 'Reco' , str(options.inputFile)[len(options.inputFile)-7: len(options.inputFile)-4])

if not options.onlyGen:
    process.load("myAnalyzers.BtoKsMuMu.Psiks0_BestPA_V0Ext_Rootupler_cfi")
else:
    process.load("myAnalyzers.BtoKsMuMu.Psiks0_OnlyGen_cfi")

process.rootuple.isMC = cms.bool(True) # this is only for test
process.rootuple.isRes = cms.bool(options.isRes)
process.rootuple.OnlyGen = cms.bool(options.onlyGen)
process.rootuple.GenParticles = cms.InputTag("prunedGenParticles") 


if options.saveInSync:
    outname = '/eos/user/g/gayalasa/Sync/CondorResults/' + outname

process.TFileService = cms.Service("TFileService",

       fileName = cms.string(outname),
)

process.p = cms.Path(process.rootuple)
 
