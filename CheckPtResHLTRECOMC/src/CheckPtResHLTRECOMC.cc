// -*- C++ -*-
//
// Package:    CheckPtResHLTRECOMC
// Class:      CheckPtResHLTRECOMC
// 
/**\class CheckPtResHLTRECOMC CheckPtResHLTRECOMC.cc PtResolutions/CheckPtResHLTRECOMC/src/CheckPtResHLTRECOMC.cc

 Description: Just extract the resolution

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jordi Duarte Campderros
//         Created:  Sun May 29 15:13:24 CEST 2011
// $Idp
//
//


// system include files
#include <memory>
#include<vector>
#include<list>
#include<map>
#include<utility>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
//#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

// Ojo si migras esto a interface/.cc --> debes llevarte esta function
template<class T>
const std::list<int> getindordlept( const T * collection, const int & id=13 )
{
	unsigned int maxsize =20;
	std::map<float,int> ordered;
	
	bool isMuon = true;

	for(unsigned int k = 0; k < collection->size(); k++)
	{
		// Only for MC particles
		if( typeid( T ) == typeid( reco::GenParticleCollection ) )
		{
			if( collection->at(k).status() == 1 && abs(collection->at(k).pdgId()) == id )
			{
				isMuon = true;
			}
			else
			{
				isMuon = false;
			}
		}
		//---- 
		if( isMuon )
		{
			ordered[collection->at(k).pt()] = k;
		}
		
		if( ordered.size() >= maxsize)
		{
			break;
		}
	}

	std::list<int> indexs;
	for(std::map<float,int>::reverse_iterator it = ordered.rbegin(); it != ordered.rend(); it++)
	{
		indexs.push_back( it->second );
	}

	return indexs;
}

//
// class declaration
//
class CheckPtResHLTRECOMC : public edm::EDAnalyzer 
{
	public:
		explicit CheckPtResHLTRECOMC(const edm::ParameterSet&);
		~CheckPtResHLTRECOMC();
		
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		void cleanDM();
		void fillhistoandmatch(const int & restype);
		std::vector<std::pair<int,int > > getindexsdRmatch( 
				const std::vector<double> & etameasv, const std::vector<double> & etagenv,
				const std::vector<double> & phimeasv, const std::vector<double> & phigenv,
				const std::vector<double> & ptmeasv, const std::vector<double> & ptgenv
				);

		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;
		
		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		
		// ----------member data ---------------------------
		enum { RECOHLT, MCRECO, MCHLT };  // GeneratedMeasured --> (pt^Gen-pt^Meas)/pt^Meas
		                                  // Equivalent to     --> (1/pt^meas-1/pt^Gen)/1/pt^Gen
		// MC particles
		std::vector<double> * _MCMupt;
		std::vector<double> * _MCMueta;
		std::vector<double> * _MCMuphi;
		// HLT particles
		std::vector<double> * _HLTMupt;
		std::vector<double> * _HLTMueta;
		std::vector<double> * _HLTMuphi;
		// PAT-RECO particles
		std::vector<double> * _RECOMupt;
		std::vector<double> * _RECOMueta;
		std::vector<double> * _RECOMuphi;
		// MATCHING INFO
		// Each key contains a vector which:
		//     * number of element: index of the Measured object
		//     * element vector content: index of the Generated
		//                     matched 
		std::map<int,std::vector<std::pair<int,int> > *> _indxMeasGen;

		// All vectors references (to facilitate the destrucction)
		std::list<std::vector<double>** > _registry;
		// Reco-Vertex and Event Content
		int _Nvrt;

		int _RunNumber;
		int _EventNumber;
		int _LuminosityBlock;
		int _Bx;
		int _Orbit;


		// Storage
		TFile * _file;
		TTree * _tree;

		std::map<int,TH2F *> _ptreshistos;

		// configuration
		bool _isData;
		std::string _tagPatTriggerEvents;
		//std::string _tagTriggerEvent;
		std::string _tagVertex;
		std::string _tagPatMuon;

};

//
// constructors and destructor
//
CheckPtResHLTRECOMC::CheckPtResHLTRECOMC(const edm::ParameterSet& iConfig)
{

	_file = new TFile( iConfig.getParameter<std::string>("outputfile").c_str(), "RECREATE" );
	_tree = new TTree( "ptres", "" );

	// Configuration
	_isData = iConfig.getParameter<bool>("isData");
	_tagPatTriggerEvents = iConfig.getParameter<std::string>("patTriggerEvents");
	//_tagTriggerEvent = iConfig.getParameter<std::string>("triggerEvents");
	_tagPatMuon = iConfig.getParameter<std::string>("patMuons");
	_tagVertex = iConfig.getParameter<std::string>("vertex");
	
	// Initialiting datamembers
	_MCMupt   =0;
	_MCMueta  =0;
        _MCMuphi  =0; 
       
        _HLTMupt  =0;
        _HLTMueta =0;
        _HLTMuphi =0;
       
        _RECOMupt =0;
        _RECOMueta=0;
        _RECOMuphi=0;
	// And registring
	_registry.push_back( &_MCMupt );
	_registry.push_back( &_MCMueta );
	_registry.push_back( &_MCMuphi );
	_registry.push_back( &_HLTMupt );
	_registry.push_back( &_HLTMueta );
	_registry.push_back( &_HLTMuphi );
	_registry.push_back( &_RECOMupt );
	_registry.push_back( &_RECOMueta );
	_registry.push_back( &_RECOMuphi );
	
	_Nvrt           = 0;
	_RunNumber      = 0;
        _EventNumber    = 0;
        _LuminosityBlock= 0;
        _Bx             = 0; 
        _Orbit          = 0;
	// Branches initialization
	_tree->Branch("MCMuPt",  &_MCMupt); 
	_tree->Branch("MCMuEta", &_MCMueta); 
	_tree->Branch("MCMuPhi", &_MCMuphi); 

	_tree->Branch("ohMuL3Pt", &_HLTMupt); 
	_tree->Branch("ohMuL3Eta",&_HLTMueta); 
	_tree->Branch("ohMuL3Phi",&_HLTMuphi); 

	_tree->Branch("recoMuonPt", &_RECOMupt); 
	_tree->Branch("recoMuonEta",&_RECOMueta); 
	_tree->Branch("recoMuonPhi",&_RECOMuphi); 

	_tree->Branch("Nvrt",&_Nvrt);


	_tree->Branch("RunNumber",&_RunNumber);
	_tree->Branch("EventNumber",&_EventNumber);
	_tree->Branch("LuminosityBlock",&_LuminosityBlock);
	_tree->Branch("BunchCrossing",&_Bx);
	_tree->Branch("OrbitNumber",&_Orbit);
	
	// Histograms and Matching
	_ptreshistos[RECOHLT] = new TH2F("recohlt", "", 250, -1.5,1.5, 100, 0,16 );
	_indxMeasGen[RECOHLT] = 0;
	_tree->Branch("idx_recohlt", &(_indxMeasGen[RECOHLT]) );
	if( ! _isData )
	{
		_ptreshistos[MCRECO] = new TH2F("mcreco","", 250, -1.5,1.5, 100, 0.0,16.0 );
		_ptreshistos[MCHLT] = new TH2F("mchlt", "", 250, -1.5,1.5, 100, 0.0,16.0 );

		_indxMeasGen[MCRECO] = 0;
		_tree->Branch("idx_mcreco", &(_indxMeasGen[MCRECO]) );
		_indxMeasGen[MCHLT]  = 0;
		_tree->Branch("idx_mchlt", &(_indxMeasGen[MCHLT]) );
	}
	// Associated to the file
	for(std::map<int,TH2F*>::iterator it = _ptreshistos.begin(); it != _ptreshistos.end(); it++)
	{
		it->second->SetDirectory(_file);
	}

}


CheckPtResHLTRECOMC::~CheckPtResHLTRECOMC()
{
	//for(std::map<int,TH2F *>::iterator it = _ptreshistos.begin(); it != _ptreshistos.end(); it++)
	//{
	//	it->second->Write();
	//}
	_file->Write();
	_file->Purge();
	_file->Close();
	_file->Delete();
	delete _file;
}


//
// member functions
//
void CheckPtResHLTRECOMC::cleanDM()
{
	for(std::list<std::vector<double>** >::iterator it = _registry.begin(); it != _registry.end(); it++)
	{
		if( *(*it) )
		{
			delete *(*it);
			*(*it) = 0;
		}
	}
	// Freeing indx
	for(std::map<int,std::vector<std::pair<int,int> > *>::iterator it = _indxMeasGen.begin();
			it != _indxMeasGen.end(); it++)
	{
		if( it->second )
		{
			delete it->second;
			it->second = 0;
		}
	}

	_Nvrt           = 0;
	_RunNumber      = 0;
        _EventNumber    = 0;
        _LuminosityBlock= 0;
        _Bx             = 0; 
        _Orbit          = 0;

}

std::vector<std::pair<int,int > > CheckPtResHLTRECOMC::getindexsdRmatch( 
		const std::vector<double> & etameasv, const std::vector<double> & etagenv,
		const std::vector<double> & phimeasv, const std::vector<double> & phigenv,
		const std::vector<double> & ptmeasv, const std::vector<double> & ptgenv
		)
{
	// FIXME: entrar-ho amb el python config
	const std::pair<double,double> range = std::pair<double,double>(6.0,28.0);
	const double dRcut = 0.05;

	std::map<int,int> igenmeasmap;
	std::map<int,double> drpairs;
	
	for(unsigned int imeas = 0; imeas < ptmeasv.size(); imeas++)
	{
		const double ptmeas = ptmeasv.at(imeas);
		if( ptmeas < range.first or ptmeas > range.second )
		{
			continue;
		}
		for(unsigned int igen = 0; igen < ptgenv.size(); igen++)
		{
			const double deta = etagenv.at(igen)-etameasv.at(imeas);
			const double dphi = phigenv.at(igen)-phimeasv.at(imeas);
			const double dR = sqrt(deta*deta+dphi*dphi);
			if( dR > dRcut )
			{
				continue;
			}
			// Desambiguation: using lesser dR
			if( igenmeasmap.count( igen ) == 0 )
			{
				igenmeasmap[igen] = imeas;
				drpairs[igen] = dR;
			}
			else
			{
				// Storing the lesser dR
				if( dR < drpairs[igen] )
				{
					igenmeasmap[igen] = imeas;
					drpairs[igen] = dR;
				}
			}
		}
	}
	
	std::vector<std::pair<int,int> > indexs;
	for(std::map<int,int>::iterator it = igenmeasmap.begin(); it != igenmeasmap.end(); it++)
	{
		// output: <measured index, generated index>
		indexs.push_back(std::pair<int,int>(it->second,it->first));
	}

	return indexs;
}

void CheckPtResHLTRECOMC::fillhistoandmatch(const int & restype)
{

	std::vector<double> * ptmeasv =0; 
	std::vector<double> * etameasv=0;
	std::vector<double> * phimeasv=0;

	std::vector<double> * ptgenv  =0; 
	std::vector<double> * etagenv =0;
	std::vector<double> * phigenv =0;

	if( restype == MCRECO )
	{
		ptgenv  = _MCMupt;
		etagenv = _MCMueta;
		phigenv = _MCMuphi;

		ptmeasv = _RECOMupt;
		etameasv= _RECOMueta;
		phimeasv= _RECOMuphi;
	}
	else if( restype == MCHLT )
	{
		ptgenv  = _MCMupt;
		etagenv = _MCMueta;
		phigenv = _MCMuphi;

		ptmeasv = _HLTMupt;
		etameasv= _HLTMueta;
		phimeasv= _HLTMuphi;
	}
	else if( restype == RECOHLT )
	{
		ptgenv  = _RECOMupt;
		etagenv = _RECOMueta;
		phigenv = _RECOMuphi;

		ptmeasv = _HLTMupt;
		etameasv= _HLTMueta;
		phimeasv= _HLTMuphi;
	}
	else
	{
		// FIXME: cms::Exception
		std::cout << "Error: Not a known resolution type '" << restype << "'" << std::endl;
		exit(-1);
	}
	// Init index ---> FIXED BUG!! (Affecting MC case)
	_indxMeasGen[restype]= new std::vector<std::pair<int,int> >;

	// meas, gen
	const std::vector<std::pair<int, int> > matchedindx = getindexsdRmatch(*etameasv,*etagenv,
			*phimeasv, *phigenv, *ptmeasv, *ptgenv );

	for(std::vector<std::pair<int,int> >::const_iterator it = matchedindx.begin();
			it != matchedindx.end(); it++)
	{
		const double ptmeas= ptmeasv->at(it->first);
		const double ptgen = ptgenv->at(it->second);

		_ptreshistos[restype]->Fill( (ptgen-ptmeas)/ptmeas, _Nvrt );

		// Matching
		_indxMeasGen[restype]->push_back(std::pair<int,int>(it->first,it->second));
	}
}

// ------------ method called for each event  ------------
void
CheckPtResHLTRECOMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	// Vertex number
	edm::Handle<std::vector<reco::Vertex> > handleVertex;
	iEvent.getByLabel( _tagVertex, handleVertex );

	_Nvrt = handleVertex->size();

	// Event Stuff
	_RunNumber      = iEvent.id().run();
	_EventNumber    = iEvent.id().event();
	_LuminosityBlock= iEvent.luminosityBlock();
	_Bx             = iEvent.bunchCrossing();
	_Orbit          = iEvent.orbitNumber();

	// HLT Stuff
	_HLTMupt = new std::vector<double>;
	_HLTMueta = new std::vector<double>;
	_HLTMuphi = new std::vector<double>;
	edm::Handle<pat::TriggerEvent > handlePatTriggerEvent;
	iEvent.getByLabel( _tagPatTriggerEvents, handlePatTriggerEvent );
	const pat::TriggerObjectCollection * triggerObj = handlePatTriggerEvent->objects() ;
	pat::TriggerObjectCollection * trMuons = new pat::TriggerObjectCollection;
	for(pat::TriggerObjectCollection::const_iterator tm = triggerObj->begin(); tm != triggerObj->end(); tm++)
	{
		if( abs(tm->pdgId()) != 13 )
		{
			continue;
		}
		trMuons->push_back( (*tm) );
		//_HLTMupt->push_back( (*tm).pt() );
		//_HLTMueta->push_back( (*tm).eta() );
		//_HLTMuphi->push_back( (*tm).phi() );
	}
	// Ordered
	const std::list<int> muhltindx = getindordlept<pat::TriggerObjectCollection>( trMuons );
	for(std::list<int>::const_iterator it = muhltindx.begin(); it != muhltindx.end(); it++)
	{
		_HLTMupt->push_back( trMuons->at(*it).pt() );
		_HLTMueta->push_back( trMuons->at(*it).eta() );
		_HLTMuphi->push_back( trMuons->at(*it).phi() );
	}
	/*edm::Handle<trigger::TriggerEvent> handleTrEvt;
	iEvent.getByLabel( _tagTriggerEvent, handleTrEvt );*/

	// RECO Stuff
	_RECOMupt = new std::vector<double>;
	_RECOMueta = new std::vector<double>;
	_RECOMuphi = new std::vector<double>;
	edm::Handle<pat::MuonCollection > handlePatMuons;
	iEvent.getByLabel( _tagPatMuon, handlePatMuons );
	const std::vector<pat::Muon> * patMuons = handlePatMuons.product();
	// ordered
	const std::list<int> muonrecoindx = getindordlept<std::vector<pat::Muon> >( patMuons );
	for(std::list<int>::const_iterator it = muonrecoindx.begin(); it != muonrecoindx.end(); it++)
	{
		_RECOMupt->push_back( patMuons->at(*it).pt() );
		_RECOMueta->push_back( patMuons->at(*it).eta() );
		_RECOMuphi->push_back( patMuons->at(*it).phi() );
	}

	// Filling the histos
	fillhistoandmatch( RECOHLT );

	// MC stuff
	if( ! _isData )
	{
		edm::Handle<reco::GenParticleCollection> genParticles;
		iEvent.getByLabel("genParticles", genParticles );
		if( ! genParticles.isValid() )
		{
			// Throw an exception better!! (find out where
			// are the exception codes??!!)
			std::cout << "Error: genParticles not found. Is this a MC sample?" << std::endl;
			exit(-1);
		}
		_MCMupt = new std::vector<double>;
		_MCMueta = new std::vector<double>;
		_MCMuphi = new std::vector<double>;
		// Get indexs of ordered gen-muons
		const std::vector<reco::GenParticle> * gpv = genParticles.product();
		const std::list<int> muonindx = getindordlept<reco::GenParticleCollection>( gpv, 13 );
		for(std::list<int>::const_iterator it = muonindx.begin(); it != muonindx.end(); it++)
		{
			_MCMupt->push_back( gpv->at(*it).pt() );
			_MCMueta->push_back( gpv->at(*it).eta() );
			_MCMuphi->push_back( gpv->at(*it).phi() );
		}

		fillhistoandmatch( MCRECO );
		fillhistoandmatch( MCHLT );
	}

	// Storing 
	_tree->Fill();
	// Freeing memory and cleaning
	cleanDM();
}

// ------------ method called once each job just before starting event loop  ------------
void 
CheckPtResHLTRECOMC::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CheckPtResHLTRECOMC::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
CheckPtResHLTRECOMC::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CheckPtResHLTRECOMC::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CheckPtResHLTRECOMC::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CheckPtResHLTRECOMC::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CheckPtResHLTRECOMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CheckPtResHLTRECOMC);
