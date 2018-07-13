from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_RECOSIM.py'
config.JobType.allowUndistributedCMSSW = False
#config.JobType.outputFiles=['step2.root']

config.section_("Data")
# MC example
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.inputDataset = '/Step1/hatake-CMSSW_10_2_0_pre6_NewGeom_Step2_v2-50c842c03a2b8d8257110386097e3d98/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#KH config.Data.totalUnits = 100
# MC example ends
# Data example
#config.Data.inputDataset = '/ZeroBias/Run2018A-v1/RAW'
#config.Data.inputDataset = '/HLTPhysics/Run2018A-v1/RAW'
#config.Data.runRange = '315361-315690'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 10
# Data example ends

config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' # Parameter Data.publishDbsUrl has been renamed to Data.publishDBS
config.Data.outputDatasetTag = 'CMSSW_10_2_0_pre6_NewGeom_Step3_v3' # <== Check!!!

config.Data.outLFNDirBase = '/store/user/hatake/crab_outputs'  # Data.outLFN has been renamed to Data.outLFNDirBase
config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T3_US_Baylor'
#KH (this whitelisting below is not really necessary. we can use any T2/T3 for running jobs. we can still send output to Baylor)
config.Site.whitelist = ['T3_US_Baylor','T2_US_*']
