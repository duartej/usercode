
/** \file HLTHiggsPlotter.cc
 *  $Date: 2011/09/07 16:31:47 $
 *  $Revision: 1.1 $
 */

#include "HLTriggerOffline/Higgs/interface/HLTHiggsPlotter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"

class DQMStore;

HLTHiggsPlotter::HLTHiggsPlotter(const ParameterSet & pset,
		const std::string & hltPath,
		const std::vector<unsigned int> objectsType,
	       	DQMStore * dbe) :
	_hltPath(hltPath),
	_objectsType(objectsType),
      	_cutsDr(pset.getParameter<std::vector<double> >("cutsDr")),
      	_parametersEta(pset.getParameter<std::vector<double> >("parametersEta")),
  	_parametersPhi(pset.getParameter<std::vector<double> >("parametersPhi")),
  	_parametersTurnOn(pset.getParameter<std::vector<double> >("parametersTurnOn")),
  	_genCut1(pset.getParameter<string>("genCut1")),
      	_recCut1(pset.getParameter<string>("recCut1")),
  	_genCut2(pset.getParameter<string>("genCut2")),
      	_recCut2(pset.getParameter<string>("recCut2")),
	_genSelector(0),
	_recSelector(0),
	_dbe(dbe)
{
}


void HLTHiggsPlotter::beginJob() 
{
}



void HLTHiggsPlotter::beginRun(const Run & iRun, const EventSetup & iSetup)
{
	static int runNumber = 0;
      	runNumber++;

	for(std::vector<unsigned int>::iterator it = _objectsType.begin(); 
			it != _objectsType.end(); ++it)
	{
		if( *it == HLTHiggsSubAnalysis::ELEC )
		{
			_cutMaxEta[*it] = 2.5; 
		}
		else
		{
			_cutMaxEta[*it] = 2.4;

			if( *it == HLTHiggsSubAnalysis::MUON &&_hltPath.find("eta2p1") != string::npos) 
			{
				_cutMaxEta[*it] = 2.1;
			}
		}

		// Choose a pT cut for gen/rec muons based on the pT cut in the _hltPath
		unsigned int threshold = 0;
		std::string objTypestr = this->getTypeString( *it );
		TPRegexp ptRegexp(std::string(objTypestr+"([0-9]+)").c_str());  
		TObjArray * regexArray = ptRegexp.MatchS(_hltPath);
		if(regexArray->GetEntriesFast() == 2) 
		{	
			threshold = atoi(((TObjString *)regexArray->At(1))->GetString());
		}
		delete regexArray;
		// We select a whole number min pT cut slightly above the _hltPath's final 
		// pt threshold, then subtract a bit to let through particle gun muons with
		// exact integer pT:
		_cutMinPt[*it] = ceil(threshold * 1.1) - 0.01;
		if(_cutMinPt[*it] < 0.) 
		{
			_cutMinPt[*it] = 0.;
		}
	 
		/*std::string baseDir = "HLT/Higgs/Distributions/";
	  	  _dbe->setCurrentFolder(baseDir + _hltPath);*/
	
		std::vector<std::string> sources(2);
		sources[0] = "gen";
		sources[1] = "rec";
	  
		for(size_t i = 0; i < sources.size(); i++) 
		{
			std::string source = sources[i];
			bookHist(source, objTypeStr, "Eta");
			bookHist(source, objTypeStr, "Phi");
			bookHist(source, objTypeStr, "MaxPt1");
			bookHist(source, objTypeStr, "MaxPt2");
		}
	}
}

void HLTHiggsPlotter::analyze(const Event & iEvent, const EventSetup & iSetup, 
		EVTColCollection * cols)
{
      	static int eventNumber = 0;
      	eventNumber++;
      	LogTrace("HLTMuonVal") << "In HLTHiggsPlotter::analyze,  " 
		<< "Event: " << eventNumber;

	
	for(std::vector<unsigned int>::iterator it = _objectsType.begin(); 
			it != _objectsType.end(); ++it)
	{
		std::vector<std::string> sources;
		bool hasGen = false;
 	      	if( (*cols).genParticles != 0 )
 		{
 			sources.push_back("gen");
 			hasGen = true;
 		}
 		bool hasReco = false;
 	      	if( col->get(*it) != 0 )
 		{
 			sources.push_back("rec");
 			hasReco = true;
 		}
 	      	for(size_t sourceNo = 0; sourceNo < sources.size(); ++sourceNo) 
 		{
 			std::string source = sources[sourceNo];
 		     	// If this is the first event, initialize selectors --> Por que no ponerlo en el beginRun
			if(!_genSelector) 
 			{
 				_genSelector = new StringCutObjectSelector<reco::GenParticle>(_genCut[*it]);
 			}
			/* FIXME:!!!!
			if( _objtypeSelRef[*it] == 0 )
			{
				this->initSelector( *it );
			}*/
			_recSelector = new StringCutObjectSelector<reco::Candidate>(_recCut[objtype]);
 		    	// Make each good gen/rec object into the base cand for a MatchStruct
 			std::vector<MatchStruct> matches;
 		    	if(source == "gen" && hasGen)
 			{
 			  	for(size_t i = 0; i < cols->genParticles->size(); ++i)
 				{
 					if(_genSelector->operator()(cols->genParticles->at(i)))
 					{
 				       		matches.push_back(MatchStruct(&cols->genParticles->at(i)));
 					}
 				}
 			}
 			if(source == "rec" && hasReco)
 			{
				std::vector<reco::Candidate> * recCol = cols->get(*it);
 			  	for(size_t i = 0; i < recCol->size(); i++)
 				{
 					if(_recSelector->operator()(recCol->at(i)))
 					{
 				      		matches.push_back(MatchStruct(&recMuons->at(i)));
 					}
 				}
 			}
 		     
 		     	// Sort the MatchStructs by pT for later filling of turn-on curve
 			std::sort(matches.begin(), matches.end(), matchesByDescendingPt());
 		    
 			const bool isDoublePath = (hltPath_.find("Double") != string::npos);
 		    	const int nObjectsToPassPath = (isDoublePath) ? 2 : 1;
 			std::vector<reco::RecoChargedCandidateRef> refsHlt;
 			std::vector<const reco::RecoChargedCandidate*> candsHlt;
 	    
 			const int hltStep = i - 1;
 			InputTag tag = InputTag("L3", "", _hltProcessName); // Ultimo????
 		  	size_t iFilter = rawTriggerEvent->filterIndex(tag);
 		  	if(iFilter < rawTriggerEvent->size()) 
 			{
				// FIXME A funcion needed to decide which trigger type is trigger::TriggerWHATEVER
 		      		rawTriggerEvent->getObjects(iFilter, trigger::TriggerMuon,refsHlt);
 		    	}
 		  	else
 			{
 				LogTrace("HLTMuonVal") << "No collection with label " << tag;
 			}
 		  
 			for(size_t j = 0; j < refsHlt.size(); ++j)
 			{
 				if(refsHlt[j].isAvailable()) 
 				{
 			      		candsHlt[i].push_back(& * refsHlt[i][j]);
 				}
 			       	else 
 				{
 					LogWarning("HLTHiggsPlotter")
 					    	<< "Ref refsHlt[j]: product not available " << j;
 				}
 		 	}
 		     
 		    	// Add trigger objects to the MatchStructs
 		    	findMatches(matches, candsHlt);
 		     
 			std::vector<size_t> matchesInEtaRange;
 			std::vector<bool> hasMatch(matches.size(), true);
 		   
 			for(size_t j=0; j < matches.size(); ++j)
 			{
 				if(matches[j].candHlt == 0)
 				{
 			    		hasMatch[j] = false;
 				}
 		      		else if(!hasMatch[j]) 
 				{
 			    		LogTrace("HLTMuonVal") 
 						<< "Match found for HLT without previous match!";
 				}
 		    		break;
 	       		}
 			
 			if(std::count(hasMatch.begin(), hasMatch.end(), true) < nObjectsToPassPath)
 			{
 				break;
 			}
 		  
 			std::string pre  = source + "Pass";
 			std::string post = "_HLT";
 		  
 		   	for(size_t j = 0; j < matches.size(); j++) 
 			{
 				float pt  = matches[j].candBase->pt();
 				float eta = matches[j].candBase->eta();
 				float phi = matches[j].candBase->phi();
 				if(hasMatch[j]) 
 				{ 
 			       		if(matchesInEtaRange.size() >= 1 && j == matchesInEtaRange[0])
 					{
 				    		elements_[pre + "MaxPt1" + post]->Fill(pt);
 					}
 			      		if(matchesInEtaRange.size() >= 2 && j == matchesInEtaRange[1])
 					{
 						elements_[pre + "MaxPt2" + post]->Fill(pt);
 					}
 			      		if(pt > cutMinPt_) 
 					{
 				    		elements_[pre + "Eta" + post]->Fill(eta);
 					}
 			    		if(fabs(eta) < cutMaxEta_)
 					{
 				  		elements_[pre + "Phi" + post]->Fill(phi);
 					}
 			 	}
 		  	}
 	      	} // End loop over sources
	} // End loop over objects
}



void HLTHiggsPlotter::findMatches(std::vector<MatchStruct> & matches,
	    	std::vector<const RecoChargedCandidate *> candsHlt)
{
  	
	std::set<size_t> indicesHlt(candsHlt.size());
    	for (size_t j = 0; j < candsHlt.size(); j++)
	{
		indicesHlt.insert(j);
	}
	
       	for(size_t i = 0; i < matches.size(); i++) 
	{
	    	const Candidate * cand = matches[i].candBase;
		
	    	double bestDeltaR = _cutsDr[0];
		size_t bestMatch = kNull;
	    	matches[i].candHlt.assign(candsHlt.size(), 0);
	    	for (size_t j = 0; j < candsHlt.size(); j++) 
		{
		   	bestDeltaR = _cutsDr[level - 2];
			bestMatch = kNull;
		  	for(std::set<size_t>::iterator it = indicesHlt.begin();
				       	it != indicesHlt[j].end(); ++it) 
			{
				double dR = deltaR(cand->eta(), cand->phi(),
		       				candsHlt[*it]->eta(), candsHlt[*it]->phi());
				if(dR < bestDeltaR) 
				{
			       		bestMatch = *it;
			      		bestDeltaR = dR;
				}
		  	}
		  	if(bestMatch != kNull)
			{
				matches[i].candHlt = candsHlt[bestMatch];
			}
		  	indicesHlt.erase(bestMatch);
	    	}
//     cout << "    Muon: " << cand->eta() << ", ";
//     if (matches[i].candL1) cout << matches[i].candL1->eta() << ", ";
//     else cout << "none, ";
//     for (size_t j = 0; j < candsHlt.size(); j++) 
//       if (matches[i].candHlt[j]) cout << matches[i].candHlt[j]->eta() << ", ";
//       else cout << "none, ";
//     cout << endl;
	}
}




void HLTHiggsPlotter::bookHist(const std:string & source, 
		const std::string & objType, const std::string & type)
{
	std::string sourceUpper = source; 
      	sourceUpper[0] = std::toupper(sourceUpper[0]);
	string name = source + objType + "Pass" + type + "_" + _hltPath;
      	TH1F * h = 0;

      	if(type.find("MaxPt") != string::npos) 
	{
		std::string desc = (type == "MaxPt1") ? "Leading" : "Next-to-Leading";
		std::string title = "pT of " + desc + " " + sourceUpper + " " + objType + " "
                   "matched to HLT";
	    	const size_t nBins = parametersTurnOn_.size() - 1;
	    	float * edges = new float[nBins + 1];
	    	for(size_t i = 0; i < nBins + 1; i++)
		{
			edges[i] = parametersTurnOn_[i];
		}
	    	h = new TH1F(name.c_str(), title.c_str(), nBins, edges);
		delete edges;
      	}
      	else 
	{
		std::string symbol = (type == "Eta") ? "#eta" : "#phi";
		std::string title  = symbol + " of " + sourceUpper + " " + objType + " "+
    			"matched to HLT";
		std::vector<double> params = (type == "Eta") ? _parametersEta : _parametersPhi;

	    	int    nBins = (int)params[0];
	    	double min   = params[1];
	    	double max   = params[2];
	    	h = new TH1F(name.c_str(), title.c_str(), nBins, min, max);
      	}
      	h->Sumw2();
      	_elements[name] = _dbe->book1D(name, h);
      	delete h;
}

void HLTHiggsPlotter::fillHist(const std::string & source, const std::string & objType, 
		const std::string & type, const float & value )
{
	std::string sourceUpper = source; 
      	sourceUpper[0] = std::toupper(sourceUpper[0]);
	string name = source + objType + "Pass" + type + "_" + _hltPath;

	_elements[name]->Fill(value);
}


void HTLHiggsPlotter::initSelector(const unsigned int & objtype)
{	
	if( objtype == HLTHiggsSubAnalyzer::MUON )
	{
		_recMuonSelector = new StringCutObjectSelector<reco::Muon>(_recCut[objtype]);
	}
	else if( objtype == HLTHiggsSubAnalyzer::ELEC )
	{
		_recElecSelector = new StringCutObjectSelector<reco::GsfElectron>(_recCut[objtype]);
	}
	else if( objtype == HLTHiggsSubAnalyzer::Photon )
	{
		_recPhotonSelector = new StringCutObjectSelector<reco::Photon>(_recCut[objtype]);
	}
	else if( objtype == HLTHiggsSubAnalyzer::JET )
	{
		_recJetSelector = new StringCutObjectSelector<reco::Jet>(_recCut[objtype]);
	}
	else if( objtype == HLTHiggsSubAnalyzer::PFJET )
	{
		_recPFJetSelector = new StringCutObjectSelector<reco::pfJet>(_recCut[objtype]);
	}
	else if( objtype == HLTHiggsSubAnalyzer::MET )
	{
		_recMETSelector = new StringCutObjectSelector<reco::caloMET>(_recCut[objtype]);
	}
	else if( objtype == HLTHiggsSubAnalyzer::PFMET )
	{
		_recPFMETSelector = new StringCutObjectSelector<reco::pfMET>(_recCut[objtype]);
	}
	/*else
	{
FIXME: ERROR NO IMPLEMENTADO
	}*/
}



//! 
const std::string HLTHiggsPlotter::getString(const unsigned int & objtype) const
{
	std::string objTypestr("Mu");

	if( objtype == HLTHiggsSubAnalyzer::ELEC )
	{
		objTypestr = "Ele";
	}
	else if( objtype == HLTHiggsSubAnalyzer::Photon )
	{
		objTypestr = "Photon";
	}
	else if( objtype == HLTHiggsSubAnalyzer::JET )
	{
		objTypestr = "Jet";
	}
	else if( objtype == HLTHiggsSubAnalyzer::PFJET )
	{
		objTypestr = "PFJet";
	}
	else if( objtype == HLTHiggsSubAnalyzer::MET )
	{
		objTypestr = "ET";
	}
	/*else
	{ ERROR FIXME
	}*/

	return objTypestr;
}
