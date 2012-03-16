#ifndef HLTRIGGEROFFLINE_HIGGS_EVTCOLCONTAINER
#define HLTRIGGEROFFLINE_HIGGS_EVTCOLCONTAINER

/** \class EVTColContainer
 *  Generate histograms for trigger efficiencies Higgs related
 *  Documentation available on the CMS TWiki:
 *  https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWGHLTValidate
 *
 *  $Date: 2012/03/15 17:53:00 $
 *  $Revision: 1.1 $
 *  \author  J. Duarte Campderros
 *
 */

#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HLTriggerOffline/Higgs/interface/HLTHiggsSubAnalysis.h"

#include<vector>
#include<map>

//! container with all the objects needed
struct EVTColContainer
{
	int nOfCollections;
	int nInitialized;
	const reco::GenParticleCollection * genParticles;
	const std::vector<reco::Candidate> * muons;
	const std::vector<reco::Candidate> * electrons;
	const std::vector<reco::Candidate> * photons;
	//const std::vector<reco::Candidate> * caloMETs;
	const std::vector<reco::Candidate> * pfMETs;
	const trigger::TriggerEventWithRefs * rawTriggerEvent;
	std::map<unsigned int,const std::vector<reco::Candidate> *> _objectsPtr;
	EVTColContainer():
		nOfCollections(4),
		nInitialized(0),
		genParticles(0),
		muons(0),
		electrons(0),
		photons(0),
		pfMETs(0),
		rawTriggerEvent(0)
	{
		_objectsPtr[HLTHiggsSubAnalysis::MUON] = muons;
		_objectsPtr[HLTHiggsSubAnalysis::ELEC] = electrons;
		_objectsPtr[HLTHiggsSubAnalysis::PHOTON] = photons;
	//	_objectsPtr[HLTHiggsSubAnalysis::CALOMET] = caloMETs;
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
		muons = 0; electrons = 0; photons = 0; pfMETs=0; //caloMETs=0; 
		rawTriggerEvent = 0;
	}
	//! Setter: multiple overloaded function
	void set(const reco::MuonCollection * v)
	{
		const void * tmp = v;
		muons = static_cast<const std::vector<reco::Candidate> *>(tmp);
		++nInitialized;
	}
	void set(const reco::GsfElectronCollection * v)
	{
		const void * tmp = v;
		electrons = static_cast<const std::vector<reco::Candidate> *>(tmp);
		++nInitialized;
	}
	void set(const reco::PhotonCollection * v)
	{
		const void * tmp = v;
		photons = static_cast<const std::vector<reco::Candidate> *>(tmp);
		++nInitialized;
	}
	/*void set(const reco::CaloMETCollection * v)
	{
		void * tmp = v;
		caloMETs = static_cast<const std::vector<reco::Candidate> *>(v);
		//caloMETs = v;
		++nInitialized;
	}*/
	void set(const std::vector<reco::PFMET> * v)
	{
		const void * tmp = v;
		pfMETs = static_cast<const std::vector<reco::Candidate> *>(tmp);
		//pfMETs = v;
		++nInitialized;
	}
	const std::vector<reco::Candidate> * get(const unsigned int & objtype)
	{
		return _objectsPtr[objtype];
	}
};
#endif
