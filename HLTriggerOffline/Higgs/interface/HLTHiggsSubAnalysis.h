#ifndef HLTriggerOffline_Higgs_HLTHiggsSubAnalysis_H
#define HLTriggerOffline_Higgs_HLTHiggsSubAnalysis_H

/** \class HLTHiggsSubAnalysis
 *  Generate histograms for trigger efficiencies Higgs related
 *  Documentation available on the CMS TWiki:
 *  https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWGHLTValidate
 *
 *  $Date: 2012/03/15 17:53:00 $
 *  $Revision: 1.1 $
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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DQMServices/Core/interface/DQMStore.h"

#include "HLTriggerOffline/Higgs/interface/HLTHiggsPlotter.h"


#include<vector>
#include<set>
#include<map>
#include<cstring>

class EVTColContainer;

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
			JET,
			PFJET,
			MHT,
			_nMAX
		};

		HLTHiggsSubAnalysis(const edm::ParameterSet & pset, 
				const std::string & analysisname );
		~HLTHiggsSubAnalysis();
	      	void beginJob();
	      	void beginRun(const edm::Run & iRun, const edm::EventSetup & iEventSetup);
	      	void analyze(const edm::Event & iEvent, const edm::EventSetup & iEventSetup, EVTColContainer * cols);

		//! Extract what objects need this analysis
		const std::vector<unsigned int> getObjectsType() const;

		
       	private:
		void bookobjects(const edm::ParameterSet & anpset);
		void initobjects(const edm::Event & iEvent, EVTColContainer * col);

		std::string _analysisname;

		edm::ParameterSet _pset;
		
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
