#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
from CRABClient.UserUtilities import config

def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = '',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.requestName = None
        #config.General.workArea = 'ZMuondecay'
        config.General.workArea = 'V5_Data'
	config.General.transferOutputs = True
	config.General.transferLogs = False

        config.JobType.pluginName = 'Analysis'
	config.JobType.psetName = '/afs/cern.ch/work/g/gayalasa/public/B0Analysis/CMSSW_10_6_12/src/myAnalyzers/BtoKsMuMu/test/MuMuks0_BestPA_V0Ext_Rootupler.py' #MC Parked configfile
	config.JobType.allowUndistributedCMSSW = True

        config.Data.inputDataset = None
	config.Data.inputDBS = 'global'
   #     config.Data.splitting = 'Automatic'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 5
   #     config.Data.totalUnits = 30
	#config.Data.lumiMask = '' # no idea 
	config.Data.publication = True
        config.Data.outputDatasetTag = None
	config.Data.outLFNDirBase = '/store/user/gayalasa/V5_UL/'
	#config.Site.storageSite = 'T3_US_FNALLPC'
	config.Site.storageSite = 'T3_CH_CERNBOX'
        #config.Site.whitelist = ['T2_US*']
        #config.Data.ignoreLocality = True
        #config.Site.storageSite = None # Choose your site. 
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        # Data taken from here: https://indico.cern.ch/event/1094697/contributions/4608916/attachments/2346139/4000698/Bstojpsiks_15_11_2021.pdf
        # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis
        inputDatasets = [ 
                          '/Charmonium/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD', # UL 2016 35.92 fb-1
                          '/Charmonium/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD', 
                          '/Charmonium/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD',
                          '/Charmonium/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD',
                          '/Charmonium/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD',
                          '/Charmonium/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD',
                          '/Charmonium/Run2016F-21Feb2020_UL2016-v1/MINIAOD',
                          '/Charmonium/Run2016G-21Feb2020_UL2016-v1/MINIAOD',
                          '/Charmonium/Run2016H-21Feb2020_UL2016-v1/MINIAOD', 

                          '/Charmonium/Run2017B-09Aug2019_UL2017-v1/MINIAOD', # UL 2017 42.42 fb-1
                          '/Charmonium/Run2017C-09Aug2019_UL2017-v1/MINIAOD',
                          '/Charmonium/Run2017D-09Aug2019_UL2017-v1/MINIAOD',
                          '/Charmonium/Run2017E-09Aug2019_UL2017-v1/MINIAOD',
                          '/Charmonium/Run2017F-09Aug2019_UL2017-v1/MINIAOD',

                          '/Charmonium/Run2018A-12Nov2019_UL2018_rsb-v1/MINIAOD',  # UL 2018 58.97 fb-1
                          '/Charmonium/Run2018B-12Nov2019_UL2018-v1/MINIAOD',
                          '/Charmonium/Run2018C-12Nov2019_UL2018_rsb_v2-v2/MINIAOD',
                          '/Charmonium/Run2018D-12Nov2019_UL2018-v1/MINIAOD'
                 	]
 
        for inDS in inputDatasets:
             # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            #config.General.requestName = inDS.split('/')[1]+'-'+inDS.split('/')[2]
            config.General.requestName = inDS.split('/')[1]+'-'+inDS.split('/')[2]
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
            # Submit.
            try:
                print "Submitting for input dataset %s" % (inDS)
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Submission for input dataset %s failed: %s" % (inDS, hte.headers)
            except ClientException as cle:
                print "Submission for input dataset %s failed: %s" % (inDS, cle)

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print "-"*len(msg)
            print msg
            print "-"*len(msg)
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
            except ClientException as cle:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)


if __name__ == '__main__':
    main()
