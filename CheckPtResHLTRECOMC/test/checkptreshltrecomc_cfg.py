import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# Global tag for data
process.GlobalTag.globaltag = 'GR_R_42_V14::All'


#- Running over AOD
from PhysicsTools.PatAlgos.tools.coreTools import *
restrictInputToAOD(process)

#- Running over data
removeMCMatching(process)

#- Remove Jets (it seems there are not present)
#- It needs to remove explicitaly the relative
#- path of Jets (it seems that the removeSpecific...
#- function does not do it)
#removeSpecificPATObjects(process,['Jets'])
#process.patDefaultSequence.remove( process.patJets )
#- Remove all but muons 
removeAllPATObjectsBut(process,['Muons'])
# -- The above function does not do the remove 
# -- of the path automaticaly for the Jets and METs
process.patDefaultSequence.remove( process.patJets )
process.patDefaultSequence.remove( process.patMETs )

#- Trigger information
# load the PAT trigger Python tools
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process) 

#================================================================================
# De moment comentat, copiat d'un notes...
# require physics declared
#process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
#process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'
#
## require scraping filter
#process.scrapingVeto = cms.EDFilter("FilterOutScraping",
#                                    applyfilter = cms.untracked.bool(True),
#                                    debugOn = cms.untracked.bool(False),
#                                    numtrack = cms.untracked.uint32(10),
#                                    thresh = cms.untracked.double(0.2)
#                                    )
## HB + HE noise filtering
#process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
#
#
#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#if mytrigs is not None :
#    process.hltSelection = hltHighLevel.clone(TriggerResultsTag = 'TriggerResults::HLT', HLTPaths = mytrigs)
#    process.hltSelection.throw = False
#
#
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                           minimumNDOF = cms.uint32(4) ,
#                                           maxAbsZ = cms.double(24),
#                                           maxd0 = cms.double(2)
#                                           )
#================================================================================



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options.wantSummary = True
#process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.source = cms.Source("PoolSource",
		# replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/7A73BEC0-4188-E011-A150-001F296BE5FA.root',
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/7C2F6CCF-4488-E011-9AE7-002481A8A782.root',
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/86AE8DCB-8087-E011-A507-00237DA1AC2A.root',
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/86BC9D61-C485-E011-BFB5-001E0B482944.root',
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/8822D251-1887-E011-A1F3-0017A4770034.root',
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/9637E316-DE86-E011-B3B3-78E7D1651098.root',
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/98870B28-6C86-E011-9348-1CC1DE047FD8.root',
			'file:/gpfs/gaes/cms/store/data/Run2011A/SingleMu/AOD/HWW-PromptSkim-v4/0000/9A202612-B985-E011-96BB-1CC1DE1CE026.root',
			)
		)


#-- Do not need the PAT output
process.outpath.remove(process.out)
#process.out.outputCommands.append( 'keep *_offlinePrimaryVertices*_*_*' )

process.ptrestuple = cms.EDAnalyzer('CheckPtResHLTRECOMC'
		, outputfile      = cms.string("testkk.root")
		, isData          = cms.bool(True)
		, vertex          = cms.string("offlinePrimaryVertices")
		, patTriggerEvents= cms.string("patTriggerEvent")
		, patMuons        = cms.string("selectedPatMuons")
#                , triggerEvent = cms.string("hltTriggerSummaryAOD")
		)
#
#
#process.p = cms.Path(process.ptrestuple)
# -- Using the default sequence from patTemplate

#--- Extract trigger information (paths, objects, ...)
#process.load("HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi")
#process.hltEventAnalyzerAOD.triggerName=cms.string("@")  # @ --> all 


process.p = cms.Path(process.patDefaultSequence*process.ptrestuple )
