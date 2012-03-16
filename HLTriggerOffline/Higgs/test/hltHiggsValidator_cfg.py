import FWCore.ParameterSet.Config as cms

process = cms.Process("HLTHiggsOfflineAnalysis")

process.load("HLTriggerOffline.Higgs.HiggsValidation_cff")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")

##############################################################################
##### Templates to change parameters in hltMuonValidator #####################
# process.hltMuonValidator.hltPathsToCheck = ["HLT_IsoMu3"]
# process.hltMuonValidator.genMuonCut = "abs(mother.pdgId) == 24"
# process.hltMuonValidator.recMuonCut = "isGlobalMuon && eta < 1.2"
##############################################################################

hltProcessName = "HLT"
process.hltHiggsValidator.hltProcessName = hltProcessName
process.hltHiggsValidator.HWW.hltPathsToCheck = cms.vstring(
		"HLT_Mu30_etqa2p1",
		"HLT_IsoMu24_eta2p1",
		"HLT_Ele27_WP80",
		)

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string(autoCond['startup'])

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	    'file:/afs/cern.ch/user/d/duarte/scratch0/step2_RAW2DIGI_RECO.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/C8D75A92-3942-E111-B75E-0018F3D09708.root', 
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/BCD347BA-3942-E111-912A-003048FFCB74.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/D65E2E4C-3A42-E111-959B-003048FFD79C.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/BC5D3B13-3B42-E111-BBD9-003048678FF8.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/8E82C028-3B42-E111-B350-0026189438D6.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/EA772D2E-3B42-E111-A566-002618943913.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/263E0428-3B42-E111-A636-0026189438C9.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/90675304-3B42-E111-BAF0-003048678E94.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/B001B713-3B42-E111-80CF-001A92971ADC.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/C27FD31B-3B42-E111-86FC-003048678C26.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/963D6675-3A42-E111-B90B-003048FFD7BE.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/E4415916-3B42-E111-8156-002618FDA279.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/B2109704-3B42-E111-A1B3-002618943911.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/6690BD27-3B42-E111-BA4B-003048678B1A.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/9A46D22D-3B42-E111-82CF-0026189437F2.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/881E757E-3B42-E111-8A95-00261894396F.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/A67D5D46-3B42-E111-91CB-003048FFCC18.root',
	    #        '/store/relval/CMSSW_5_2_0_pre1/RelValH130GGgluonfusion/GEN-SIM-DIGI-RECO/START50_V9_FastSim-v1/0010/2299587A-3D42-E111-AF97-002618943867.root',
    ),
    secondaryFileNames = cms.untracked.vstring(
	    'file:/afs/cern.ch/user/d/duarte/scratch0/H130GGgluonfusion_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root',
    )
)

process.DQMStore = cms.Service("DQMStore")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.MessageLogger.destinations += ['HiggsValidationMessages']
process.MessageLogger.categories   += ['HiggsValidation']
process.MessageLogger.debugModules += ['*']#HLTHiggsValidator','HLTHiggsSubAnalysis','HLTHiggsPlotter']
process.MessageLogger.HiggsValidationMessages = cms.untracked.PSet(
    threshold       = cms.untracked.string('DEBUG'),
    default         = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    HiggsValidation = cms.untracked.PSet(limit = cms.untracked.int32(1000))
    )

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *', 
        'keep *_MEtoEDMConverter_*_HLTMuonOfflineAnalysis'),
    fileName = cms.untracked.string('hltHiggsValidator.root')
)

process.analyzerpath = cms.Path(
    process.hltHiggsValidator *
    process.MEtoEDMConverter
)

process.outpath = cms.EndPath(process.out)
