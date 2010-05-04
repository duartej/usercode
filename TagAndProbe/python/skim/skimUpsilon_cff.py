import FWCore.ParameterSet.Config as cms

ptMinCut       = 'pt > 2 || (abs(eta) > 1 && p > 2)';
triggerPath1Mu = 'HLT_Mu3'                 # Single-muon Trigger name
triggerFilt1Mu = 'hltSingleMu3L3Filtered3' # Single-muon Trigger last filter name
triggerPath2Mu = 'HLT_L1DoubleMuOpen'                  # Double-muon Trigger name
triggerFilt2Mu = 'hltDoubleMuLevel1PathL1OpenFiltered' # Double-muon last filter name
massRangeMu  = (8.0, 12.0)
massRangeTk  = (8.0, 12.0) 
massRangeSta = (5.0, 14.0)


##     ____  _ _                     _                  _    ___  ____  
##    / ___|| (_)_ __ ___  _ __ ___ (_)_ __   __ _     / \  / _ \|  _ \ 
##    \___ \| | | '_ ` _ \| '_ ` _ \| | '_ \ / _` |   / _ \| | | | | | |
##     ___) | | | | | | | | | | | | | | | | | (_| |  / ___ \ |_| | |_| |
##    |____/|_|_|_| |_| |_|_| |_| |_|_|_| |_|\__, | /_/   \_\___/|____/ 
##                                           |___/                      
## ==== pruner of GenParticles ====
genMuons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *  ",                     # this is the default
        "++keep abs(pdgId) = 13",        # keep muons and their parents
        "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
    )
)
## ==== Just ignore all the too low pt stuff ====
goodTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string(ptMinCut),
)

slimAOD = cms.Sequence(
    genMuons +
    goodTracks 
)


##    __  __       _          ____   _  _____   __  __                       
##   |  \/  | __ _| | _____  |  _ \ / \|_   _| |  \/  |_   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \ | |_) / _ \ | |   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ |  __/ ___ \| |   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___| |_| /_/   \_\_|   |_|  |_|\__,_|\___/|_| |_|___/
##
##                                                                           
##  We take them from Onia2MuMuPAT but we make a few changes
##
from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import *

## Trigger matches
## put a PT cut on the muons
patMuons.cut = ptMinCut;
## use the genMuons as MC source, so that we keep them and have the correct mother refs
muonMatch.matched = 'genMuons'
## switch off genParticle embedding, as we keep the genMuons collection
## also switch off standalone muon track embedding, as we keep it separately
patMuonsWithoutTrigger.embedGenMatch = False
patMuonsWithoutTrigger.embedStandAloneMuon = False

##    __  __       _                _   _                                                 
##   |  \/  | __ _| | _____    ___ | |_| |__   ___ _ __   _ __ ___  _   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \  / _ \| __| '_ \ / _ \ '__| | '_ ` _ \| | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ | (_) | |_| | | |  __/ |    | | | | | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___|  \___/ \__|_| |_|\___|_|    |_| |_| |_|\__,_|\___/|_| |_|___/
##                                                                                        
##   
## ==== For bare tracks, make candidates assuming the muon mass hypothesis ====
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("goodTracks"),      
    particleType = cms.string("mu+"),
) 
## ==== For skimming with standalone muons, use the raw standalone track, so that it doesn't take the tracker P4 if it's also a global or tracker muon ====
staTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("standAloneMuons","UpdatedAtVtx"), 
    particleType = cms.string("mu+"),
)
## ==== Muons that are not only standalone (for skimming only) ====
nonStaOnlyMuon = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull"),
)

## ==== Golden muons to be used for tags. use PAT ones, so I can check HLT =====
PASS_HLT_1 = "!triggerObjectMatchesByFilter('%s').empty()" % (triggerFilt1Mu,);
PASS_HLT_2 = "!triggerObjectMatchesByFilter('%s').empty()" % (triggerFilt2Mu,);
tags1Mu = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && "+ PASS_HLT_1), 
)
tags2Mu = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && " + PASS_HLT_2), 
)

otherMuons = cms.Sequence(
    nonStaOnlyMuon +
    tags1Mu +
    tags2Mu +
    tkTracks +
    staTracks
)

allMuons = cms.Sequence(
    patMuonSequence *
    otherMuons
)

##    __  __       _           _   _ ___       _ _             
##   |  \/  | __ _| | _____   | | | | _ \ ___ (_| | ___  _  _  
##   | |\/| |/ _` | |/ / _ \   \ | /||_) / __|| | |/ _ \| \| |
##   | |  | | (_| |   <  __/    | | | _ /\__ \| | | (_| |    | 
##   |_|  |_|\__,_|_|\_\___|    |_| |_|  |__ /|_|_|\__ /|_|\_| 
##                                                        
##   
upsilonMu  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tags1Mu@+ nonStaOnlyMuon@-"),
    cut = cms.string("%f < mass < %f" % massRangeMu),
)
upsilonTk  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tags1Mu@+ tkTracks@-"),
    cut = cms.string("%f < mass < %f" % massRangeTk),
)
upsilonSta = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tags2Mu@+ staTracks@-"),
    cut = cms.string("%f < mass < %f" % massRangeSta),
)
allUpsilons = cms.Sequence(
    upsilonMu + upsilonTk + upsilonSta
)

##    __  __       _                     _   _     
##   |  \/  | __ _(_)_ __    _ __   __ _| |_| |__  
##   | |\/| |/ _` | | '_ \  | '_ \ / _` | __| '_ \ 
##   | |  | | (_| | | | | | | |_) | (_| | |_| | | |
##   |_|  |_|\__,_|_|_| |_| | .__/ \__,_|\__|_| |_|
##                          |_|                    
##   
upsilonSkimMainSequence = cms.Sequence(
    slimAOD  *    
    allMuons *
    allUpsilons
)

##    ____  _    _             ____       _   _         
##   / ___|| | _(_)_ __ ___   |  _ \ __ _| |_| |__  ___ 
##   \___ \| |/ / | '_ ` _ \  | |_) / _` | __| '_ \/ __|
##    ___) |   <| | | | | | | |  __/ (_| | |_| | | \__ \
##   |____/|_|\_\_|_| |_| |_| |_|   \__,_|\__|_| |_|___/
##                                                      
##   
upsilonMuFilter  = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("upsilonMu"),
    minNumber = cms.uint32(1),
)
upsilonTkFilter  = upsilonMuFilter.clone(src = 'upsilonTk')
upsilonStaFilter = upsilonMuFilter.clone(src = 'upsilonSta')


##     ___        _               _   
##    / _ \ _   _| |_ _ __  _   _| |_ 
##   | | | | | | | __| '_ \| | | | __|
##   | |_| | |_| | |_| |_) | |_| | |_ 
##    \___/ \__,_|\__| .__/ \__,_|\__|
##                   |_|              
##   
upsilonSkimOut = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("skimUpsilon.root"),
    #fileName = cms.untracked.string("/tmp/gpetrucc/skimUpsilon_ppMuX.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_patMuons_*_Skim",
        "keep *_goodTracks_*_Skim",
        "keep *_genMuons_*_Skim",
        "keep recoTrackExtras_standAloneMuons_*_*",          ## track states at the muon system, used both by patMuons and standAloneMuons
        "keep recoTracks_standAloneMuons_UpdatedAtVtx_*",    ## bare standalone muon tracks, using standalone muon momentum (with BS constraint)
        "keep edmTriggerResults_*_*_Skim",                   ## to know which kind of skim channel got us the event   
        "keep edmTriggerResults_*_*_HLT",                    ## to know what got us on tape
        "keep l1extraL1MuonParticles_l1extraParticles_*_*",  ## if we ever want to do L1 efficiency too ## <<--- Not in 3.1.X AODSIM
        "keep *_offlinePrimaryVertices__*",                  ## vertices and BS are not very useful on MC
        "keep *_offlineBeamSpot__*",                         ## but they can be important on data
       #"keep *_upsilonMu_*_Skim", "keep *_upsilonTk_*_Skim", "keep *_upsilonSta_*_Skim",                       ## <<--- keep these for monitoring
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring(
        "Skim_upsilonMu",
        "Skim_upsilonTk",
        "Skim_upsilonSta",
    )),
)

#def Summer09_Trigger(process):
#    process.patTrigger.processName = 'HLT8E29'
#    if hasattr(process, 'filterHLT1Mu'): process.filterHLT1Mu.TriggerResultsTag = cms.InputTag('TriggerResults','','HLT8E29')
#    if hasattr(process, 'filterHLT2Mu'): process.filterHLT2Mu.TriggerResultsTag = cms.InputTag('TriggerResults','','HLT8E29')
#    ## Remove matches to l1Extra which is not in 3.1.X samples and that we woudln't use anyway
#    for n in process.patTrigger.parameters_(): 
#        if n.startswith("l1Extra"): delattr(process.patTrigger, n)
#
#def Spring10ReDigi_Trigger(process):
#    process.patTrigger.processName = 'REDIGI'
#    if hasattr(process, 'muonMatchHLTCtfTrack'):
#        process.muonMatchHLTCtfTrack.collectionTags[0] = process.muonMatchHLTCtfTrack.collectionTags[0].replace('::HLT','::REDIGI')
#    if hasattr(process, 'filterHLT1Mu'): process.filterHLT1Mu.TriggerResultsTag = cms.InputTag('TriggerResults','','REDIGI')
#    if hasattr(process, 'filterHLT2Mu'): process.filterHLT2Mu.TriggerResultsTag = cms.InputTag('TriggerResults','','REDIGI')
#    ## Remove matches to l1Extra which is not in 3.1.X samples and that we woudln't use anyway
#    for n in process.patTrigger.parameters_():
#        if n.startswith("l1Extra"): delattr(process.patTrigger, n)

