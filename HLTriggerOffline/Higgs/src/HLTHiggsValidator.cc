// -*- C++ -*-
//
// Package:     HLTHiggsValidator
// Class:       HLTHiggsValidator
// 

//
// Jordi Duarte Campderros (based on the Jason Slaunwhite 
// and Jeff Klukas coded from the HLTriggerOffline/Muon package
//
// $Id: HLTHiggsValidator.cc,v 1.1 2012/03/12 16:34:41 duarte Exp $
//
//

// system include files
//#include<memory>
#include<iostream>

#include "HLTriggerOffline/Higgs/interface/HLTHiggsValidator.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TPRegexp.h"
//////// Class Methods ///////////////////////////////////////////////////////
// Constructor
HLTHiggsValidator::HLTHiggsValidator(const edm::ParameterSet& pset) :
      	_pset(pset),
	_hltProcessName(pset.getParameter<string>("hltProcessName")),
	_analysisnames(pset.getParameter<vstring>("analysis")),
	_collections(0),
	_dbe(0)
{
}


void HLTHiggsValidator::beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup) 
{
      	// Initialize hltConfig
      	bool changedConfig;
	if(!_hltConfig.init(iRun, iSetup, _hltProcessName, changedConfig)) 
	{
		edm::LogError("HLTHiggsValidator") 
			<< "Initialization of HLTConfigProvider failed!!"; 
	    	return;
      	}
	
      	// Get the set of trigger paths we want to make plots for
	/*std::set<std::string> hltPaths;
      	for(size_t i = 0; i < _hltPathsToCheck.size(); ++i) 
	{
	     	TPRegexp pattern(_hltPathsToCheck[i]);
	    	for(size_t j = 0; j < _hltConfig.triggerNames().size(); ++j)
		{
			if(TString(_hltConfig.triggerNames()[j]).Contains(pattern))
			{
				hltPaths.insert(_hltConfig.triggerNames()[j]);
			}
		}
      	}*/

	//std::vector<std::string> hltPathsV(hltPaths.begin(),hltPaths.end());
       
	for(size_t i = 0; i < _analysisnames.size() < ++i)
	{
		HLTHiggsSubAnalysis analyzer(pset_, _analysisnames.at(i));
		_analyzers.push_back(analyzer);
	}
	// Call the Plotter beginRun (which stores the triggers paths..:)
      	for(std::vector<HLTHiggsSubAnalysis>::iterator iter = _analyzers.begin(); 
			iter != _analyzers.end(); ++iter) 
	{
	    	iter->beginRun(iRun, iSetup);
	
	}
}
	

void HLTHiggsValidator::analyze(const Event& iEvent, const EventSetup& iSetup)
{
      	for(std::vector<HLTHiggsSubAnalysis>::iterator iter = _analyzers.begin(); 
			iter != _analyzers.end(); ++iter) 
	{
	     	iter->analyze(iEvent, iSetup,this->_collections);
      	}
	
	// Deleting after analysis
	if( this->_collections != 0)
	{
		delete this->_collections;
		this->_collections = 0;
	}
}



void HLTHiggsValidator::beginJob()
{
}

void HLTHiggsValidator::endRun(const edm::Run & iRun, const edm::EventSetup& iSetup)
{
      	// vector<HLTMuonPlotter>::iterator iter;
      	// for(std::vector<HLTHiggsPlotter>::iterator iter = _analyzers.begin(); 
	//                 iter != analyzers_.end(); ++iter) 
	// {
      	//         iter->endRun(iRun, iSetup);
      	// }
}


void HLTHiggsValidator::endJob()
{
}



//define this as a plug-in
DEFINE_FWK_MODULE(HLTHiggsValidator);
