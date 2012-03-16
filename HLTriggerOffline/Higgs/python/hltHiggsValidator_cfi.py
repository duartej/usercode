import FWCore.ParameterSet.Config as cms


hltHiggsValidator = cms.EDAnalyzer("HLTHiggsValidator",
		
    hltProcessName = cms.string("HLT"),
    analysis       = cms.vstring("HWW"), #,"H2tau", "HZZ", "Hgg", "Htaunu",
    
    genParticleLabel = cms.string("genParticles"),
    parametersTurnOn = cms.vdouble(0,
                                   1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                   22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
                                   45, 50, 55, 60, 65, 70,
                                   80, 100, 200, 500, 1000, 2000,
                                   ), 

    parametersEta      = cms.vdouble(48, -2.400, 2.400),
    parametersPhi      = cms.vdouble(50, -3.142, 3.142),

    # set criteria for matching at L1, L2, L3
    cutsDr = cms.vdouble(0.4, 0.4, 0.015),
    # parameters for attempting an L1 match using a propagator
    maxDeltaPhi = cms.double(0.4),
    maxDeltaR   = cms.double(0.4),

    # set cuts on generated and reconstructed objects 
    # --- Muons
    Mu_genCut     = cms.string("abs(pdgId) == 13 && status == 1"),
    Mu_recCut     = cms.string("isGlobalMuon"),
    Mu_cutMinPt   = cms.double(10),
    Mu_cutMaxEta  = cms.double(2.4),
    Mu_cutMotherId= cms.uint32(23), # CUIDDADO
    Mu_cutDr      = cms.vdouble(0.4,0.4,0.015),
    
    # --- Electrons
    Ele_genCut      = cms.string("abs(pdgId) == 11 && status == 1"),
    Ele_recCut      = cms.string("hadronicOverEm < 0.05 && eSuperClusterOverP > 0.5 && eSuperClusterOverP < 2.5"),
    Ele_cutMinPt    = cms.double(10),
    Ele_cutMaxEta   = cms.double(2.5),
    Ele_cutMotherId = cms.uint32(23), # CUIDDADO
    Ele_cutDr       = cms.vdouble(0.4,0.4,0.015),


    # The specific parameters per analysis
    HWW = cms.PSet( 
	    hltPathsToCheck = cms.vstring(
		    "HLT_(L[12])?(Double)?(Iso)?Mu[0-9]*(Open)?(_NoVertex)?(_eta2p1)?(_v[0-9]*)?$",
		    "HLT_Dimuon0_Jpsi_v10",
		    "HLT_Dimuon13_Jpsi_Barrel_v5",
		    ),
	    recMuonLabel = cms.string("muons"),
	    recElecLabel = cms.string("gsfElectrons"),
	    ),
)
