
/** \file HLTHiggsSubAnalysis.cc
 *  $Date: 2012/03/16 01:55:33 $
 *  $Revision: 1.2 $
 */





#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"

#include "HLTriggerOffline/Higgs/interface/HLTHiggsSubAnalysis.h"
#include "HLTriggerOffline/Higgs/src/EVTColContainer.cc"

#include "TPRegexp.h"

#include "TString.h"

#include<set>


HLTHiggsSubAnalysis::HLTHiggsSubAnalysis(const edm::ParameterSet & pset,
		const std::string & analysisname) :
	_pset(pset),
	_analysisname(analysisname),
	_hltProcessName(pset.getParameter<std::string>("hltProcessName")),
	_genParticleLabel(pset.getParameter<std::string>("genParticleLabel")),
	_dbe(0)
{
	edm::ParameterSet anpset = pset.getParameter<edm::ParameterSet>(analysisname);
	//FIXME: CHECK ERRORS

	// Collections labels (but genparticles already initialized)
	this->bookobjects( anpset );

	_hltPathsToCheck = anpset.getParameter<std::vector<std::string> >("hltPathsToCheck");

	_dbe = edm::Service<DQMStore>().operator->();
      	_dbe->setVerbose(0);
}

HLTHiggsSubAnalysis::~HLTHiggsSubAnalysis()
{
}


void HLTHiggsSubAnalysis::beginJob() 
{
}



void HLTHiggsSubAnalysis::beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup)
{
	std::string baseDir = "HLT/Higgs/"+_analysisname+"/";
      	_dbe->setCurrentFolder(baseDir);

	// Initialize the confighlt
	bool changedConfig;
	if(!_hltConfig.init(iRun,iSetup,_hltProcessName,changedConfig))
	{
		edm::LogError("HLTHiggsSubAnalysis") 
			<< "Initializtion of HLTConfigProvider failed!!";
	}


	// Parse the paths and initi
	_hltPaths.clear();
	for(size_t i = 0; i < _hltPathsToCheck.size(); ++i)
	{
		TPRegexp pattern(_hltPathsToCheck[i]);
		for(size_t j = 0 ; j < _hltConfig.triggerNames().size(); ++j)
		{
			std::string thetriggername = _hltConfig.triggerNames()[j];
			if(TString(thetriggername).Contains(pattern))
			{
				_hltPaths.insert(thetriggername);
			}
		}
	}
      	// Initialize the plotters (analysers for each trigger path)
	_analyzers.clear();
  	for(std::set<std::string>::iterator iPath = _hltPaths.begin(); 
			iPath != _hltPaths.end(); iPath++) 
	{
		std::string path = * iPath;
		std::string shortpath = path;
	    	if(path.rfind("_v") < path.length())
		{
			shortpath = path.substr(0, path.rfind("_v"));
		}
      				
		HLTHiggsPlotter analyzer(_pset, shortpath, this->getObjectsType(shortpath), _dbe);
		_analyzers.push_back(analyzer);
    	}
      	// Call the beginRun (which books all the histograms)
      	for(std::vector<HLTHiggsPlotter>::iterator it = _analyzers.begin(); 
			it != _analyzers.end(); ++it) 
	{
	    	it->beginRun(iRun, iSetup);
	}
}



void HLTHiggsSubAnalysis::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup, 
		EVTColContainer * collections)
{
      	static int eventNumber = 0;
      	eventNumber++;
      	LogTrace("HLTHiggsValidation") << "In HLTHiggsSubAnalysis::analyze,  " 
		<< "Event: " << eventNumber;

	// Initialize the collection (the ones which hasn't been initiliazed yet)
 	this->initobjects(iEvent,collections);
	
	for(std::vector<HLTHiggsPlotter>::iterator it = _analyzers.begin();
			it != _analyzers.end(); ++it)
	{
		it->analyze(iEvent,iSetup,collections);
	}  
}

const std::vector<unsigned int> HLTHiggsSubAnalysis::getObjectsType(const std::string & hltPath) const
{
	std::vector<unsigned int> objsType;
	for(std::map<unsigned int,std::string>::const_iterator it = _recLabels.begin();
			it != _recLabels.end(); ++it)
	{
		std::string objTypeStr = this->getTypeString( it->first );
		// Check if it is needed this object for this trigger
		if( ! TString(hltPath).Contains(objTypeStr) )
		{
			continue;
		}
		objsType.push_back(it->first);
	}
	return objsType;
}

const std::vector<unsigned int> HLTHiggsSubAnalysis::getObjectsType() const // TO BE DEPRECATED
{
	std::vector<unsigned int> objsType;
	for(std::map<unsigned int,std::string>::const_iterator it = _recLabels.begin();
			it != _recLabels.end(); ++it)
	{
		objsType.push_back(it->first);
	}
	return objsType;
}



void HLTHiggsSubAnalysis::bookobjects( const edm::ParameterSet & anpset )
{
	if( anpset.exists("recMuonLabel") )
	{
		_recLabels[MUON] = anpset.getParameter<std::string>("recMuonLabel");
	}
	if( anpset.exists("recElecLabel") )
	{
		_recLabels[ELEC] = anpset.getParameter<std::string>("recElecLabel");
	}
	if( anpset.exists("recPhotonLabel") )
	{
		_recLabels[PHOTON] = anpset.getParameter<std::string>("recPhotonLabel");
	}
	if( anpset.exists("recMETLabel") )
	{
		_recLabels[MET] = anpset.getParameter<std::string>("recMETLabel");
	}
	if( anpset.exists("recPFMETLabel") )
	{
		_recLabels[PFMET] = anpset.getParameter<std::string>("recPFMETLabel");
	}
	if( anpset.exists("recPFTauLabel") )
	{
		_recLabels[PFTAU] = anpset.getParameter<std::string>("recPFTauLabel");
	}
	if( anpset.exists("recPFJetLabel") )
	{
		_recLabels[PFJET] = anpset.getParameter<std::string>("recPFJetLabel");
	}
	if( anpset.exists("recMHTLabel") )
	{
		_recLabels[MHT] = anpset.getParameter<std::string>("recMHTLabel");
	}

	if( _recLabels.size() < 1 )
	{
		edm::LogError("HLTHiggsVal") << "In HLTHiggsSubAnalysis::bookobjects, " 
		<< "Not included any object (recMuonLabel, recElecLabel, ...)  "
	       	<< "in the analysis " << _analysisname;
		return;
	}
}

void HLTHiggsSubAnalysis::initobjects(const edm::Event & iEvent, EVTColContainer * col)
{
	if( col != 0 && col->isAllInit() )
	{
		// Already init, not needed to do nothing
		return;
	}
	if( ! col->isCommonInit() )
	{
		edm::Handle<trigger::TriggerEventWithRefs> rawTEH;
		iEvent.getByLabel("hltTriggerSummaryRAW",rawTEH);
		if(rawTEH.failedToGet())
		{
			edm::LogError("HLTMuonVal") << "No trigger summary found"; 
			return;
		}
		col->rawTriggerEvent = rawTEH.product();

		edm::Handle<reco::GenParticleCollection> genPart;
		iEvent.getByLabel(_genParticleLabel,genPart);
		if( genPart.isValid() )
		{
			col->genParticles = genPart.product();
		}
	}
		
		
	for(std::map<unsigned int,std::string>::iterator it = _recLabels.begin(); 
			it != _recLabels.end(); ++it)
	{
		if( it->first == MUON )
		{
			edm::Handle<reco::MuonCollection> theHandle;
			iEvent.getByLabel(it->second, theHandle);
			col->set(theHandle.product());
		}
		else if( it->first == ELEC )
		{
			edm::Handle<reco::GsfElectronCollection> theHandle;
			iEvent.getByLabel(it->second, theHandle);
			col->set(theHandle.product());
		}
		else
		{
			edm::LogError("HLTHiggsVal") << "HLTHiggsSubAnalysis::initobjects " 
				<< " NOT IMPLEMENTED (yet) ERROR: '" << it->second << "'";
			//return; ??
		}
	}
}

const std::string HLTHiggsSubAnalysis::getTypeString(const unsigned int & objtype) const
{
	std::string objTypestr("Mu");

	if( objtype == HLTHiggsSubAnalysis::ELEC )
	{
		objTypestr = "Ele";
	}
	else if( objtype == HLTHiggsSubAnalysis::PHOTON )
	{
		objTypestr = "Photon";
	}
	else if( objtype == HLTHiggsSubAnalysis::JET )
	{
		objTypestr = "Jet";
	}
	else if( objtype == HLTHiggsSubAnalysis::PFJET )
	{
		objTypestr = "PFJet";
	}
	else if( objtype == HLTHiggsSubAnalysis::MET )
	{
		objTypestr = "ET";
	}
	/*else
	{ ERROR FIXME
	}*/

	return objTypestr;
}
