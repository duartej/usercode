import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('CheckPtResHLTRECOMC'
		, outputfile = cms.string("output.root")
		, isData     = cms.bool(True)
		)
