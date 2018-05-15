# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase2_realistic -n 10 --era Phase2 --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry Extended2023D17 --filein file:step2.root --fileout file:step3.root
#------------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
import FWCore.ParameterSet.VarParsing as VarParsing

#------------------------------------------------------------------------------------
# Declare the process and input variables
#------------------------------------------------------------------------------------
#process = cms.Process('NOISE',eras.Run2_50ns)#for 50ns 13 TeV data
#process = cms.Process('NOISE',eras.Run2_25ns)#for 25ns 13 TeV data
options = VarParsing.VarParsing ('analysis')
process = cms.Process('RECO',eras.Phase2)

options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")
#
# Dataset
# /RelValSinglePiPt25Eta1p7_2p7/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/GEN-SIM-RECO
# /SinglePiPt*Eta1p6_2p8/PhaseIITDRFall17*93X_upgrade2023_realistic_v2*/GEN-SIM-RECO
#
# pt=5 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt5Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/18D5CC99-22AA-E711-8020-90B11C08CDC7.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt5Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/B8AEECCD-22AA-E711-82FF-1866DA7F8E98.root'
#options.outputFile = 'results_pt5.root'
#
# pt=10 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt10Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/0844EA3B-14A9-E711-9F1D-44A842CFD633.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt10Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/8CF1BCF2-3AA9-E711-B0BE-A4BF0112BD74.root'
#options.outputFile = 'results_pt10.root'
#
# pt=15 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt15Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/1E9C9D26-4EA9-E711-B775-6CC2173DA2F0.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt15Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/5C395DE1-ACA9-E711-B391-0CC47A7EEE32.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt15Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/6AD5611D-81A9-E711-92A7-3417EBE7009F.root'
#options.outputFile = 'results_pt15.root'
#
# pt=25 GeV sample *relval*
options.inputFiles = '/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2CC28E03-BDA6-E711-815F-0CC47A4D767E.root','/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/7A0E458B-C1A6-E711-9ED9-0025905B8604.root','/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/F880EB0D-C1A6-E711-AC1A-0CC47A4C8F18.root'
options.outputFile = 'step3_pt25.root'
#
# pt=50 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17GS/SinglePiPt50Eta1p6_2p8/GEN-SIM/93X_upgrade2023_realistic_v2-v1/50000/3C70EDD8-2CA4-E711-AFF1-7CD30AD091F0.root','/store/mc/PhaseIITDRFall17GS/SinglePiPt50Eta1p6_2p8/GEN-SIM/93X_upgrade2023_realistic_v2-v1/50000/60E9D7F3-71A3-E711-937B-001E67E6F7E7.root','/store/mc/PhaseIITDRFall17GS/SinglePiPt50Eta1p6_2p8/GEN-SIM/93X_upgrade2023_realistic_v2-v1/50000/9C2B2E60-4AA3-E711-AB55-6C3BE5B59200.root','/store/mc/PhaseIITDRFall17GS/SinglePiPt50Eta1p6_2p8/GEN-SIM/93X_upgrade2023_realistic_v2-v1/50000/C434E268-4DA3-E711-844F-0CC47A57CC42.root'
#options.outputFile = 'results_pt50.root'
#
# pt=100 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt100Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/BE82B842-31AE-E711-9EA4-0026B94DBDA2.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt100Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/E20ADE15-A5AE-E711-9434-0023AEEEB55F.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt100Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/FC41ACC6-2FAE-E711-B431-F04DA2747854.root'
#options.outputFile = 'results_pt100.root'
#
#'/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2637F672-C7A6-E711-B4EF-0025905A612A.root'
options.maxEvents = -1 # -1 means all events
options.skipEvents = 0 # default is 0.

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
options.parseArguments()
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(options.secondaryInputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents) # default is 0.
)

# process.TFileService = cms.Service("TFileService", 
#                                    fileName = cms.string(options.outputFile)
# )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
    fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)

#------------------------------------------------------------------------------------
# import of standard configurations
#------------------------------------------------------------------------------------
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#KH
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# #------------------------------------------------------------------------------------
# # Set up our analyzer
# #------------------------------------------------------------------------------------
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_Tree_cfi")
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_Event_cfi")
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_GenParticles_cfi")
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_HBHERecHits_cfi")
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_HGCRecHits_cfi")
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_HGCSimHits_cfi")
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_SimTracks_cfi")
# process.load("HGCalAnalysis.HGCalTreeMaker.HGCalTupleMaker_RecoTracks_cfi")

# process.load("Validation.HGCalValidation.hgcalHitValidation_cfi")

#------------------------------------------------------------------------------------
# Output definition
#------------------------------------------------------------------------------------

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)
#KH keep Uncalibrated Rechits
#process.FEVTDEBUGHLToutput.outputCommands = cms.untracked.vstring('keep *_HGCalUncalibRecHit_*_*')
process.FEVTDEBUGHLToutput.outputCommands.extend( ['keep *_HGCalUncalibRecHit_*_*'] )

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# #------------------------------------------------------------------------------------
# # HGCalTupleMaker sequence definition
# #------------------------------------------------------------------------------------
# process.tuple_step = cms.Sequence(
#     # Make HCAL tuples: Event, run, ls number
#     process.hgcalTupleTree*
#     process.hgcalTupleEvent
#     # process.hgcalTupleHBHERecHits*
#     # process.hgcalTupleHGCRecHits*
#     # process.hgcalTupleGenParticles*
#     # process.hgcalTupleHGCSimHits*
#     # process.hgcalTupleSimTracks*
#     # process.hgcalTupleGeneralTracks
# )

# #-----------------------------------------------------------------------------------
# # Path and EndPath definitions
# #-----------------------------------------------------------------------------------
# process.myTask = cms.Task()
# #process.myTask.add(*[getattr(process,prod) for prod in process.producers_()])
# #process.myTask.add(*[getattr(process,filt) for filt in process.filters_()])
# #process.myTask.add(process.hgcalTupleTree)
# process.myTask.add(process.hgcalTupleEvent)
# process.preparation = cms.Path(
#     #process.hgcalHitValidation*
#     #process.digireco_step
#     process.tuple_step
# )
# process.preparation.associate(process.myTask)

#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)

process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

#-----------------------------------------------------------------------------------
# Schedule definition
#-----------------------------------------------------------------------------------
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    #process.preparation,
    process.FEVTDEBUGHLToutput_step
)

#KH process.schedule.associate(process.myTask)

#-----------------------------------------------------------------------------------

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

