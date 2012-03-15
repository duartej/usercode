
/** \file HLTHiggsSubAnalysis.cc
 *  $Date: 2011/09/07 16:31:47 $
 *  $Revision: 1.1 $
 */



#include "HLTriggerOffline/Higgs/interface/HLTHiggsSubAnalysis.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates//interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates//interface/GsfElectronFwd.h"

#include<set>

HLTHiggsSubAnalysis::HLTHiggsSubAnalysis(const edm::ParameterSet & pset,
		const std::string & analysisname) :
	_genParticleLabel(pset.getParameter<vstring>("genParticleLabel")),
	_analysisname(analysisname)
{
	edm::ParameterSet anpset = pset.getParameter<edm::ParameterSet>(analysisname);
	//FIXME: CHECK ERRORS

	// Collections labels (but genparticles already initialized)
	this->bookobjects( anpset );

	_hltPathsToCheck = anpset.getParameter<vstring>("hltPathsToCheck");
	_dbe = Service<DQMStore>().operator->();
      	_dbe->setVerbose(0);
}

HLTHiggsSubAnalysis::~HLTHiggsSubAnalysis()
{
}


void HLTHiggsSubAnalysis::beginJob() 
{
}



void HLTHiggsSubAnalysis::beginRun(const Run & iRun, const EventSetup & iSetup)
{
	std::string baseDir = "HLT/Higgs/"+_analysisname+"/";
      	_dbe->setCurrentFolder(baseDir);

      	// Initialize the plotters (analysers for each trigger path)
	_analyzers.clear();
  	for(std::vector<std::string>::iterator iPath = _hltPaths.begin(); 
			iPath != _hltPaths.end(); iPath++) 
	{
		std::string path = * iPath;
		std::string shortpath = path;
	    	if(path.rfind("_v") < path.length())
		{
			shortpath = path.substr(0, path.rfind("_v"));
		}
      				
		HLTHiggsPlotter analyzer(pset_, shortpath, this->getObjectsType(), _dbe);
		_analyzers.push_back(analyzer);
    	}
      	// Call the beginRun (which books all the histograms)
      	for (std::vector<HLTHiggsPlotter>::iterator it = _analyzers.begin(); 
			it != _analyzers.end(); ++iter) 
	{
	    	it->beginRun(iRun, iSetup);
	}
}



void HLTHiggsSubAnalysis::analyze(const Event & iEvent, const EventSetup & iSetup, 
		EVTColContainer * collections)
{
      	static int eventNumber = 0;
      	eventNumber++;
      	LogTrace("HLTHiggsVal") << "In HLTHiggsSubAnalysis::analyze,  " 
		<< "Event: " << eventNumber;

	// Initialize the collection (the ones which hasn't been initiliazed yet)
 	this->initobjects(collections);
	
	for(std::vector<HLTHiggsPlotter>::iterator it = _analyzers.begin();
			it != _analyzers.end(); ++it)
	{
		it->analyze(iEvent,iSetup,collections);
	}  
}


const std::vector<unsigned int> getObjectsType() const
{
	std::vector<unsigned int> objsType;
	for(std::map<unsigned int,std::string>::iterator it = _recLabel.begin();
			it != _recLabel.end(); ++it)
	{
		objsType.push_back(it->first);
	}
	return objsType;
}



void HLTHiggsSubAnalysis::bookobjects( const edm::ParameterSet & anpset )
{
	if( anpset.exists("recMuonLabel") )
	{
		_recLabels[MUON] = pset.getParameter<vstring>("recMuonLabel");
	}
	if( anpset.exists("recElecLabel") )
	{
		_recLabels[ELEC] = pset.getParameter<vstring>("recElecLabel");
	}
	if( anpset.exists("recPhotonLabel") )
	{
		_recLabel[PHOTHON] = pset.getParameter<vstring>("recPhotonLabel");
	}
	if( anpset.exists("recMETLabel") )
	{
		_recLabel[MET] = pset.getParameter<vstring>("recMETLabel");
	}
	if( anpset.exists("recPFMETLabel") )
	{
		_recLabels[PFMET] = pset.getParameter<vstring>("recPFMETLabel");
	}
	if( anpset.exists("recPFTauLabel") )
	{
		_recLabels[PFTAU] = pset.getParameter<vstring>("recPFTauLabel");
	}
	if( anpset.exists("recPFJetLabel") )
	{
		_recLabels[PFJET] = pset.getParameter<vstring>("recPFJetLabel");
	}
	if( anpset.exists("recMHTLabel") )
	{
		_recLabels[MHT] = pset.getParameter<vstring>("recMHTLabel");
	}

	if( _recLabels.size() < 1 )
	{
		LogError("HLTHiggsVal") << "In HLTHiggsSubAnalysis::bookobjects, " 
		<< "Not included any object (recMuonLabel, recElecLabel, ...)  "
	       	<< "in the analysis " << _analysisname;
		return;
	}
}

void HLTHiggsSubAnalysis::initobjects(EVTColContainer * col)
{
	if( col != 0 && col->isAllInit() )
	{
		// Already init, not needed to do nothing
		return;
	}
	if( ! col.isCommonInit() )
	{
		edm::Handle<trigger::TriggerEventWithRefs> rawTEH;
		iEvent.getByLabel("hltTriggerSummaryRAW",rawTEH);
		if(rawTriggerEvent.failedToGet())
		{
			LogError("HLTMuonVal") << "No trigger summary found"; 
			return;
		}
		(*col).rawTriggerEvent = rawTEH.getproduct();

		edm::Handle<reco::GenParticleCollection> genPart;
		iEvent.getByLabel(_genParticleLabel,genPart);
		if( genPart.isValid() )
		{
			(*col).genParticles = genPart.product()
		}
	}
		
		
	for(std::map<unsigned int,std::string>::iterator it = _recLabels.begin(); 
			it != _recLabels.end(); ++it)
	{
		if( it->first == MUON )
		{
			edm::Handle<reco::MuonCollection> theHandle;
			iEvent.getByLabel(it->second, theHandle);
		}
		else if( it->first == ELEC )
		{
			edm::Handle<reco::GsfElectronCollection> theHandle;
			iEvent.getByLabel(it->second, theHandle);
		}
		else
		{
			LogError("HLTHiggsVal") << "HLTHiggsSubAnalysis::initobjects " 
				<< " NOT IMPLEMENTED (yet) ERROR: '" << it->second << "'";
			//return; ??
		}

		col->set(theHandle.product());
	}
}
