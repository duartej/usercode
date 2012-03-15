#ifndef HLTriggerOffline_Higgs_HLTHiggsSubAnalysis_H
#define HLTriggerOffline_Higgs_HLTHiggsSubAnalysis_H

/** \class HLTHiggsSubAnalysis
 *  Generate histograms for trigger efficiencies Higgs related
 *  Documentation available on the CMS TWiki:
 *  https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWGHLTValidate
 *
 *  $Date: 2012/03/12 16:10:32 $
 *  $Revision: 1.0 $
 *  \author  J. Duarte Campderros (based and adapted on J. Klukas,
 *           M. Vander Donckt and J. Alcaraz code from the 
 *           HLTriggerOffline/Muon package)
 *
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "DQMServices/Core/interface/DQMStore.h"


#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>
#include <cctype>

#include "TPRegexp.h"

//! container with all the objects needed
struct EVTColContainer
{
	const reco::GenParticleCollection * genParticles;
	/*const reco::MuonCollection * muons;
	const reco::GsfElectronCollection * electrons;
	const reco::PhotonCollection * photons;
	const reco::CaloMETCollection * caloMETs;
	const std::vector<reco::PFMET> * pfMETs;*/
	const std::vector<reco::Candidate> * muons;
	const std::vector<reco::Candidate> * electrons;
	const std::vector<reco::Candidate> * photons;
	const std::vector<reco::Candidate> * caloMETs;
	const std::vector<reco::Candidate> * pfMETs;
	const trigger::TriggerEventWithRefs rawTriggerEvent;
	const int nOfCollections = 5;
	int nInitialized;
	//std::map<unsigned int,void *> _objectsPtr;
	std::map<unsigned int,std::vector<reco::Candidate> *> _objectsPtr;
	EVTColContainer():
		nInitialized(0)
	{
		genParticles = 0;
		muons = 0; electrons = 0; photons = 0; caloMETs=0; pfMETs=0;
		rawTriggerEvent = 0;
		_objectsPtr[HLTHiggsSubAnalysis::MUON] = muons;
		_objectsPtr[HLTHiggsSubAnalysis::ELEC] = electrons;
		_objectsPtr[HLTHiggsSubAnalysis::PHOTONS] = photons;
		_objectsPtr[HLTHiggsSubAnalysis::CALOMET] = caloMETs;
		_objectsPtr[HLTHiggsSubAnalysis::PFMET] = pfMETs;
	}
	//! 
	bool isAllInit()
	{
		return (nInitialized == nOfCollections);
	}
	bool isCommonInit()
	{
		return (rawTriggerEvent == 0);
	}
	//! 
	void reset()
	{
		nInitialized = 0;
		genParticles = 0;
		muons = 0; electrons = 0; photons = 0; caloMETs=0; pfMETs=0;
		rawTriggerEvent = 0;
	}
	//! Setter: multiple overloaded function
	void set(const reco::MuonCollection * v)
	{
		//muons = v;
		muons = static_cast<std::vector<reco::Candidate> *>(v);
		++nInitialized;
	}
	void set(const reco::GsfElectronCollection * v)
	{
		//electrons = v;
		electrons = static_cast<std::vector<reco::Candidate> *>(v);
		++nInitialized;
	}
	void set(const reco::PhotonCollection * v)
	{
		photons = static_cast<std::vector<reco::Candidate> *>(v);
		//photons = v;
		++nInitialized;
	}
	void set(const reco::CaloMETCollection * v)
	{
		caloMETs = static_cast<std::vector<reco::Candidate> *>(v);
		//caloMETs = v;
		++nInitialized;
	}
	void set(const std::vector<reco::PFMET> * v)
	{
		pfMETs = static_cast<std::vector<reco::Candidate> *>(v);
		//pfMETs = v;
		++nInitialized;
	}
	//! Getter: just to get the pointer
	/*void * get(const unsigned int & objtype) const
	{
		return _objectsPtr[objtype];
	}*/
	//std::vector<reco::Candidate> * getcollection(const unsigned int & objtype) const
	std::vector<reco::Candidate> * get(const unsigned int & objtype) const
	{
		return static_cast<std::vector<reco::Candidate> *>(_objectsPtr[objtype]);
	}
};

class HLTHiggsSubAnalysis 
{
       	public:
		enum 
		{
			MUON,
			ELEC,
			PHOTON,
			MET,
			PFMET,
			PFTAU,
			PFJET,
			MHT,
			_nMAX
		};

		HLTHiggsSubAnalysis(const edm::ParameterSet & pset, 
				const std::string & analysisname,
				const std::vector<std::string> & hltPaths);
		~HLTHiggsSubAnalysis();
	      	void beginJob();
	      	void beginRun(const edm::Run &, const edm::EventSetup &);
	      	void analyze(const edm::Event &, const edm::EventSetup &, EVTColContainer * colls);

		//! Extract what objects need this analysis
		const std::vector<unsigned int> getObjectsType() const;

		
       	private:
		void bookobjects(EVTColCollections * col);
		void initobjects(EVTColCollections * col);

		std::string _analysisname;
		
		std::vector<std::string> _hltPathsToCheck;

		// The name of the object collections to be used in this
		// analysis. 
	      	std::string _genParticleLabel;
		std::map<unsigned int,std::string> _recLabels;

		// The plotters: managers of each hlt path where the plots are done
		std::vector<HLTHiggsPlotter> _analyzers;
		
	      	DQMStore* _dbe;
};
#endif
