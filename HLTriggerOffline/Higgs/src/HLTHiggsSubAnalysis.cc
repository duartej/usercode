
/** \file HLTHiggsSubAnalysis.cc
 *  $Date: 2012/03/16 01:55:33 $
 *  $Revision: 1.2 $
 */

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"

#include "HLTriggerOffline/Higgs/interface/HLTHiggsSubAnalysis.h"
#include "HLTriggerOffline/Higgs/src/EVTColContainer.cc"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTriggerOffline/Higgs/interface/HLTHiggsSubAnalysis.h"
#include "HLTriggerOffline/Higgs/src/MatchStruct.cc"

#include "TPRegexp.h"
#include "TString.h"

#include<set>
#include<algorithm>


HLTHiggsSubAnalysis::HLTHiggsSubAnalysis(const edm::ParameterSet & pset,
		const std::string & analysisname) :
	_pset(pset),
	_analysisname(analysisname),
	_minCandidates(0),
	_hltProcessName(pset.getParameter<std::string>("hltProcessName")),
	_genParticleLabel(pset.getParameter<std::string>("genParticleLabel")),
      	_parametersEta(pset.getParameter<std::vector<double> >("parametersEta")),
  	_parametersPhi(pset.getParameter<std::vector<double> >("parametersPhi")),
  	_parametersTurnOn(pset.getParameter<std::vector<double> >("parametersTurnOn")),
	_genSelector(0),
	_recMuonSelector(0),
	_recElecSelector(0),
	/*_recMETSelector(0),
	_recPFMETSelector(0),
	_recJetSelector(0),
	_recPFJetSelector(0),*/
	_recPhotonSelector(0),
	_dbe(0)
{
	for(std::map<unsigned int,std::string>::const_iterator it = _recLabels.begin();
			it != _recLabels.end(); ++it)
	{
		const std::string objStr = this->getTypeString(it->first);
		_genCut[it->first] = pset.getParameter<std::string>( std::string(objStr+"_genCut").c_str() );
		_recCut[it->first] = pset.getParameter<std::string>( std::string(objStr+"_recCut").c_str() );
		_cutMinPt[it->first] = pset.getParameter<double>( std::string(objStr+"_cutMinPt").c_str() );
		_cutMaxEta[it->first] = pset.getParameter<double>( std::string(objStr+"_cutMaxEta").c_str() );
		_cutMotherId[it->first] = pset.getParameter<unsigned int>( std::string(objStr+"_cutMotherId").c_str() );
		_cutsDr[it->first] = pset.getParameter<std::vector<double> >( std::string(objStr+"_cutDr").c_str() );
	}
	// Specific parameters for this analysis
	edm::ParameterSet anpset = pset.getParameter<edm::ParameterSet>(analysisname);

	// Collections labels (but genparticles already initialized)
	this->bookobjects( anpset );

	//FIXME: CHECK ERRORS
	_hltPathsToCheck = anpset.getParameter<std::vector<std::string> >("hltPathsToCheck");
	_minCandidates = anpset.getParameter<unsigned int>("minCandidates");

	_dbe = edm::Service<DQMStore>().operator->();
      	_dbe->setVerbose(0);
}

HLTHiggsSubAnalysis::~HLTHiggsSubAnalysis()
{
	if( _genSelector != 0)
	{
		delete _genSelector;
		_genSelector =0;
	}
	if( _recMuonSelector != 0)
	{
		delete _recMuonSelector;
		_recMuonSelector =0;
	}
	if( _recElecSelector != 0)
	{
		delete _recElecSelector;
		_recElecSelector =0;
	}
	if( _recPhotonSelector != 0)
	{
		delete _recPhotonSelector;
		_recPhotonSelector =0;
	}
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


	// Parse the input paths to get them if there are in the table 
	// and associate them the last filter of the path (in order to extract the
	// trigger objects to match with the reco/gen objects)
	_hltPaths.clear();
	std::map<std::string,std::string> pathLastFilter;
	for(size_t i = 0; i < _hltPathsToCheck.size(); ++i)
	{
		TPRegexp pattern(_hltPathsToCheck[i]);
		for(size_t j = 0 ; j < _hltConfig.triggerNames().size(); ++j)
		{
			std::string thetriggername = _hltConfig.triggerNames()[j];
			if(TString(thetriggername).Contains(pattern))
			{
				_hltPaths.insert(thetriggername);
				std::vector<std::string> mods = _hltConfig.moduleLabels(thetriggername);
				// Extract the last filter (to get the trigger objects 
				// which passing the filter)
				pathLastFilter[thetriggername] = mods.at(mods.size()-2); 
			}
		}
	}

      	// Initialize the plotters (analysers for each trigger path)
	_analyzers.clear();
  	for(std::set<std::string>::iterator iPath = _hltPaths.begin(); 
			iPath != _hltPaths.end(); ++iPath) 
	{
		std::string path = * iPath;
		std::string shortpath = path;
	    	if(path.rfind("_v") < path.length())
		{
			shortpath = path.substr(0, path.rfind("_v"));
		}
		_shortpath2long[shortpath] = path;

		// Sanity check: the object needed by a trigger path  was
		// introduced by the user via config python (_recLabels datamember)
		const std::vector<unsigned int> objsNeedHLT = this->getObjectsType(shortpath);
		if( objsNeedHLT.size() == 0 )
		{
			edm::LogError("HiggsValidation") << "In HLTHiggsSubAnalysis::beginRun, " 
				<< "Incoherence found in the python configuration file!!\nThe SubAnalysis '" 
				<< _analysisname << "' has been asked to evaluate the trigger path '"
				<< shortpath << "' (found it in 'hltPathsToCheck') BUT this path"
				<< " needs a variable which has not been instantiate ('recVariableLabels'" 
				<< ")" ;
			exit(-1);
		}

		// the hlt path, its last filter, the objects (elec,muons,photons,...)
		// needed to evaluate the path are the argumens of the plotter
		HLTHiggsPlotter analyzer(_pset, shortpath, pathLastFilter[path],
				objsNeedHLT, _dbe);
		_analyzers.push_back(analyzer);
    	}

      	// Call the beginRun (which books all the path dependent histograms)
      	for(std::vector<HLTHiggsPlotter>::iterator it = _analyzers.begin(); 
			it != _analyzers.end(); ++it) 
	{
	    	it->beginRun(iRun, iSetup);
	}

	// Book the gen/reco analysis-dependent histograms (denominators)
	for(std::map<unsigned int,std::string>::const_iterator it = _recLabels.begin();
			it != _recLabels.end(); ++it)
	{
		const std::string objStr = this->getTypeString(it->first);
		std::vector<std::string> sources(2);
		sources[0] = "gen";
		sources[1] = "rec";
	  
		for(size_t i = 0; i < sources.size(); i++) 
		{
			std::string source = sources[i];
			bookHist(source, objStr, "Eta");
			bookHist(source, objStr, "Phi");
			bookHist(source, objStr, "MaxPt1");
			bookHist(source, objStr, "MaxPt2");
		}
	}
}



void HLTHiggsSubAnalysis::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup, 
		EVTColContainer * cols)
{
      	static int eventNumber = 0;
      	eventNumber++;
      	LogTrace("HLTHiggsValidation") << "In HLTHiggsSubAnalysis::analyze,  " 
		<< "Event: " << eventNumber;

	// Initialize the collection (the ones which hasn't been initialiazed yet)
 	this->initobjects(iEvent,cols);

	//! Map to reference the object type with its selector
	std::map<unsigned int, void *> recObjSelRef;
	// Map to reference the source (gen/reco) with the recoCandidates
	std::map<unsigned int,std::vector<MatchStruct> > sourceMatchMap;
	// utility map
	std::map<unsigned int,std::string> u2str;
	u2str[GEN]="gen";
	u2str[RECO]="rec";
	// Extract the match structure containing the gen/reco candidates (electron, muons,...)
	// common to all the SubAnalysis
	for(std::map<unsigned int,std::string>::iterator it = _recLabels.begin(); 
			it != _recLabels.end(); ++it)
	{
		std::vector<unsigned int> sources;
		bool hasGen = false;
		if( cols->genParticles != 0 )
		{
			// using the mc particles to evaluate the efficiencies
			sources.push_back(GEN);
			hasGen = true;
		}

		bool hasReco = false;
		if(cols->get(it->first) != 0)
		{
			// Using the reco particles to evaluate the efficiencies
			sources.push_back(RECO);
			hasReco = true;
		}

		// Extraction of the gen/reco candidates 
		for(size_t sourceNo = 0; sourceNo < sources.size(); ++sourceNo)
		{
			const unsigned int source = sources[sourceNo];
			// Initialize selectors when first event
			if(!_genSelector) 
 			{
 				_genSelector      = new StringCutObjectSelector<reco::GenParticle>(_genCut[it->first]);
			}
			if(!recObjSelRef[it->first])
			{
				recObjSelRef[it->first]= this->InitSelector(it->first);
			}
 		    	// Make each good gen/rec object into the base cand for a MatchStruct
 			std::vector<MatchStruct> matches;
 		    	if(source == GEN && hasGen)
 			{
 			  	for(size_t i = 0; i < cols->genParticles->size(); ++i)
 				{
 					if(_genSelector->operator()(cols->genParticles->at(i)))
 					{
 				       		matches.push_back(MatchStruct(&cols->genParticles->at(i),it->first));
 					}
 				}
 			}
 			if(source == RECO && hasReco)
 			{
				const std::vector<reco::Candidate> * recCol = cols->get(it->first);
				StringCutObjectSelector<reco::Candidate> * recSelector = 
					static_cast<StringCutObjectSelector<reco::Candidate>* >(recObjSelRef[it->first]);
 			  	for(size_t i = 0; i < recCol->size(); i++)
 				{
 					if(recSelector->operator()(recCol->at(i)))
 					{
 				      		matches.push_back(MatchStruct(&recCol->at(i),it->first));
 					}
 				}
 			}
			// Sort the MatchStructs by pT for later filling of turn-on curve
			std::sort(matches.begin(), matches.end(), matchesByDescendingPt());
			sourceMatchMap[source] = matches;
		}
	}
	
	// Filling the histograms if pass the minimum amount of candidates needed by the analysis
	for(std::map<unsigned int,std::vector<MatchStruct> >::iterator it = sourceMatchMap.begin(); 
			it != sourceMatchMap.end(); ++it)
	{
		// it->first: gen/reco   it->second: matches (std::vector<MatchStruc>)
		if( it->second.size() < _minCandidates )
		{
			continue;
		}
		
		// Filling the gen/reco objects (eff-denominators): 
		// Just the first two if there are more
		for(size_t j = 0; j < it->second.size(); ++j)
		{
			std::string objTypeStr = this->getTypeString((it->second)[j].objType);

			float pt  = (it->second)[j].candBase->pt();
			float eta = (it->second)[j].candBase->eta();
			float phi = (it->second)[j].candBase->phi();
			this->fillHist(u2str[it->first],objTypeStr,"Eta",eta);
			this->fillHist(u2str[it->first],objTypeStr,"Phi",phi);
			if( j == 0 )
			{
				this->fillHist(u2str[it->first],objTypeStr,"MaxPt1",pt);
			}
			else if( j == 1 )
			{
				this->fillHist(u2str[it->first],objTypeStr,"MaxPt2",pt);
				break;
			}
		}
	}
	
	// Setting up the trigger: EVTColContainer tiene que cogerlo...
	edm::Handle<edm::TriggerResults> trigResults;
	edm::InputTag trigResultsTag("TriggerResults","",cols->rawTriggerEvent->usedProcessName());
        iEvent.getByLabel(trigResultsTag,trigResults);

	const edm::TriggerNames trigNames = iEvent.triggerNames(*trigResults);
	//const edm::TriggerNames trigNames = iEvent.triggerNames(cols->triggerResults);


	// Calling to the plotters analysis (where the evaluation of the different trigger paths are done)
	for(std::map<unsigned int,std::vector<MatchStruct> >::iterator um = sourceMatchMap.begin();
			um != sourceMatchMap.end(); ++um)
	{
		// gen/rec
		std::string source = u2str[um->first];
		for(std::vector<HLTHiggsPlotter>::iterator it = _analyzers.begin();
					it != _analyzers.end(); ++it)
		{
			const std::string hltPath = _shortpath2long[it->gethltpath()];
			//const bool ispassTrigger =  cols->triggerResults->accept(trigNames.triggerIndex(hltPath));
			const bool ispassTrigger =  trigResults->accept(trigNames.triggerIndex(hltPath));
			it->analyze(ispassTrigger,source,um->second);
		}
	}
}

// Return the objects (muons,electrons,photons,...) needed by a hlt path. Note that it 
// returns a vector which can contain repeated elements if it is a double path for instance
const std::vector<unsigned int> HLTHiggsSubAnalysis::getObjectsType(const std::string & hltPath) const
{
	std::vector<unsigned int> objsType;
	// The object to deal has to be entered via the config .py
	for(std::map<unsigned int,std::string>::const_iterator it = _recLabels.begin();
			it != _recLabels.end(); ++it)
	{
		std::string objTypeStr = this->getTypeString( it->first );
		// Check if it is needed this object for this trigger
		if( ! TString(hltPath).Contains(objTypeStr) )
		{
			continue;
		}

		// Number of ocurrences of the object (so how many objects we need)
		int n = 0;
		std::string::size_type pos = 0;
		while( (pos = hltPath.find( objTypeStr, pos)) != std::string::npos )
		{
			++n;
			pos += objTypeStr.size();
		} // End algorithm to extract n-ocurrences
		int nocur = n;
		for(int i = 0; i < nocur; ++i)
		{
			objsType.push_back(it->first);
		}
		// Or (and) if has the double label in the name
		if( TString(hltPath).Contains("Double") )
		{
			objsType.push_back(it->first);
		}
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
	
		//edm::Handle<edm::TriggerResults> trigResults;
		//edm::InputTag trigResultsTag("TriggerResults","",col->rawTriggerEvent->usedProcessName());
		//iEvent.getByLabel(trigResultsTag,trigResults);

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
		else if( it->first == PHOTON )
		{
			edm::Handle<reco::PhotonCollection> theHandle;
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

void HLTHiggsSubAnalysis::bookHist(const std::string & source, 
		const std::string & objType, const std::string & variable)
{
	std::string sourceUpper = source; 
      	sourceUpper[0] = std::toupper(sourceUpper[0]);
	std::string name = source + objType + variable ;
      	TH1F * h = 0;

      	if(variable.find("MaxPt") != std::string::npos) 
	{
		std::string desc = (variable == "MaxPt1") ? "Leading" : "Next-to-Leading";
		std::string title = "pT of " + desc + " " + sourceUpper + " " + objType;
	    	const size_t nBins = _parametersTurnOn.size() - 1;
	    	float * edges = new float[nBins + 1];
	    	for(size_t i = 0; i < nBins + 1; i++)
		{
			edges[i] = _parametersTurnOn[i];
		}
	    	h = new TH1F(name.c_str(), title.c_str(), nBins, edges);
		delete edges;
      	}
      	else 
	{
		std::string symbol = (variable == "Eta") ? "#eta" : "#phi";
		std::string title  = symbol + " of " + sourceUpper + " " + objType;
		std::vector<double> params = (variable == "Eta") ? _parametersEta : _parametersPhi;

	    	int    nBins = (int)params[0];
	    	double min   = params[1];
	    	double max   = params[2];
	    	h = new TH1F(name.c_str(), title.c_str(), nBins, min, max);
      	}
      	h->Sumw2();
      	_elements[name] = _dbe->book1D(name, h);
      	delete h;
}

void HLTHiggsSubAnalysis::fillHist(const std::string & source, 
		const std::string & objType, const std::string & variable, const float & value )
{
	std::string sourceUpper = source; 
      	sourceUpper[0] = toupper(sourceUpper[0]);
	std::string name = source + objType + variable ;

	_elements[name]->Fill(value);
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



void * HLTHiggsSubAnalysis::InitSelector(const unsigned int & objtype)
{	
	void * selector = 0;

	if( objtype == HLTHiggsSubAnalysis::MUON )
	{
		_recMuonSelector = new StringCutObjectSelector<reco::Muon>(_recCut[objtype]);
		selector = _recMuonSelector;
	}
	else if( objtype == HLTHiggsSubAnalysis::ELEC )
	{
		_recElecSelector = new StringCutObjectSelector<reco::GsfElectron>(_recCut[objtype]);
		selector = _recElecSelector;
	}
	else if( objtype == HLTHiggsSubAnalysis::PHOTON )
	{
		_recPhotonSelector = new StringCutObjectSelector<reco::Photon>(_recCut[objtype]);
		selector = _recPhotonSelector;
	}
/*	else if( objtype == HLTHiggsSubAnalysis::JET )
	{
		_recJetSelector = new StringCutObjectSelector<reco::Jet>(_recCut[objtype]);
		selector = _recMuonSelector;
	}
	else if( objtype == HLTHiggsSubAnalysis::PFJET )
	{
		_recPFJetSelector = new StringCutObjectSelector<reco::pfJet>(_recCut[objtype]);
		selector = _recMuonSelector;
	}
	else if( objtype == HLTHiggsSubAnalysis::MET )
	{
		_recMETSelector = new StringCutObjectSelector<reco::caloMET>(_recCut[objtype]);
		selector = _recMuonSelector;
	}
	else if( objtype == HLTHiggsSubAnalysis::PFMET )
	{
		_recPFMETSelector = new StringCutObjectSelector<reco::pfMET>(_recCut[objtype]);
		selector = _recMuonSelector;
	}*/
/*	else
	{
FIXME: ERROR NO IMPLEMENTADO
	}*/

	return selector;
}
