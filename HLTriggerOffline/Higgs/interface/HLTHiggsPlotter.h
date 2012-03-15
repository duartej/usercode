#ifndef HLTriggerOffline_Higgs_HLTHiggsPlotter_H
#define HLTriggerOffline_Higgs_HLTHiggsPlotter_H

/** \class HLTHiggsPlotter
 *  Generate histograms for trigger efficiencies Higgs related
 *  Documentation available on the CMS TWiki:
 *  https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWGHLTValidate
 *
 *  $Date: 2012/03/12 16:10:32 $
 *  $Revision: 1.0 $
 *  \author  J. Duarte Campderros (based and adapted on J. Klukas,
 *           M. Vander Donckt and J. Alcaraz code from the 
 *           HLTriggerOffline/Muon package)
 *  \author  J. Klukas, M. Vander Donckt, J. Alcaraz
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"


#include <vector>
#include <cstring>
#include <map>
#include <iostream>


class HLTHiggsPlotter 
{
       	public:
	      	HLTHiggsPlotter(const edm::ParameterSet & pset, const std::string & hltPath,
	       			const std::vector<std::string> & moduleLabels,
			       	const std::vector<std::string> & stepLabels);
	      	void beginJob();
	      	void beginRun(const edm::Run &, const edm::EventSetup &);
	      	void analyze(const edm::Event &, const edm::EventSetup &);
		
       	private:
	      	
	      	struct MatchStruct 
		{
		    	const reco::Candidate            * candBase;
		    	std::vector<const reco::RecoChargedCandidate *> candHlt;
		    	MatchStruct() 
			{
			  	candBase   = 0;
		    	}
		    	MatchStruct(const reco::Candidate * cand) 
			{
			  	candBase = cand;
		    	}
		    	bool operator<(MatchStruct match) 
			{      
				return candBase->pt() < match.candBase->pt();
		    	}
		    	bool operator>(MatchStruct match) 
			{
			  	return candBase->pt() > match.candBase->pt();
		    	}
	      	};
		
	      	struct matchesByDescendingPt 
		{
		     	bool operator() (MatchStruct a, MatchStruct b) 
			{     
			       	return a.candBase->pt() > b.candBase->pt();
		    	}
	      	};

	      	void analyzePath(const edm::Event &,
	     			const std::string &, const std::string &,
	     			const std::vector<MatchStruct>, 
	     			edm::Handle<trigger::TriggerEventWithRefs>);
	      	void findMatches(std::vector<MatchStruct> &, 
			  	std::vector< std::vector< const reco::RecoChargedCandidate *> >
			  	);
	      	void bookHist(std::string, std::string, std::string, std::string);
	      	void fillHist(std::string, std::string, std::string, std::string);

		std::string _subanalysisname;
		
	      	std::string  _hltPath;
		
	      	std::vector<double> _parametersEta;
	      	std::vector<double> _parametersPhi;
	      	std::vector<double> _parametersTurnOn;
		
		std::map<unsigned int,double> _cutMinPt;
		std::map<unsigned int,double> _cutMaxEta;
		std::map<unsigned int,unsigned int> _cutMotherId;
		std::map<unsigned int,std::vector<double> > _cutsDr;
		std::map<unsigned int,std::string> _genCut;
		std::map<unsigned int,std::string> _recCut;
		
	      	StringCutObjectSelector<reco::GenParticle> * _genSelector;
	      	StringCutObjectSelector<reco::Candidate>   * _recSelector;
	      	/*StringCutObjectSelector<reco::Muon>        * _recMuonSelector;
	      	StringCutObjectSelector<reco::GsfElectron> * _recElecSelector;
	      	StringCutObjectSelector<reco::caloMET>     * _recMETSelector;
	      	StringCutObjectSelector<reco::pfMET>       * _recPFMETSelector;
	      	StringCutObjectSelector<reco::Photon>      * _recPhotonSelector;*/
	      	

		//! Map to reference the object type with its selector
		std::map<unsigned int, void *> _objtypeSelRef;
		
	      	HLTConfigProvider _hltConfig;
		
	      	DQMStore* _dbe;
	      	std::map<std::string, MonitorElement *> _elements;		
};
#endif
