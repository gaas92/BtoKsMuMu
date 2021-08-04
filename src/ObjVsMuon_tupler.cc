// Taken from https://github.com/ocerri/BPH_RDntuplizer/blob/master/plugins/ObjVsMuon_tupler.cc
// Original Author: Olmo Cerri 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <iostream>
#include <string>
#include <regex>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1I.h"

using namespace std;

class ObjVsMuon_tupler : public edm::EDAnalyzer {
   public:
      explicit ObjVsMuon_tupler(const edm::ParameterSet& iConfig);
      ~ObjVsMuon_tupler() {};


   private:
      //void beginJob(const edm::EventSetup&) {};
      //void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
      //void addToTree();

      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
      edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;

      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
      edm::EDGetTokenT<reco::BeamSpot> BSLabel_;

      edm::Service<TFileService> fs;

      TTree* tree;
      map<string, float> outMap;
      bool treeDeclared = false;

      bool isRealData;
      unsigned int runNum;
      unsigned int lumiNum;
      unsigned long long eventNum;

      vector<string> triggerTags = {"Mu12_IP6", 
	                                  "Mu10p5_IP3p5", "Mu9_IP0", "Mu9_IP3", "Mu9_IP4", "Mu9_IP5", "Mu9_IP6", 
                                    "Mu8p5_IP3p5",  "Mu8_IP3", "Mu8_IP5", "Mu8_IP6", 
                                    "Mu7_IP4"};

      std::vector<float>       *TriggerObj_pt, *TriggerObj_eta, *TriggerObj_phi, *TriggerObj_ch, *TriggerObj_ip;
      std::vector<int>         *obj_HLT_Mu7_IP4, *obj_HLT_Mu8_IP3, *obj_HLT_Mu8_IP5, *obj_HLT_Mu8_IP6, *obj_HLT_Mu8p5_IP3p5;
      std::vector<int>         *obj_HLT_Mu9_IP0, *obj_HLT_Mu9_IP3, *obj_HLT_Mu9_IP4, *obj_HLT_Mu9_IP5, *obj_HLT_Mu9_IP6, *obj_HLT_Mu10p5_IP3p5, *obj_HLT_Mu12_IP6;

      std::vector<float>       *TriggerMu_pt, *TriggerMu_eta, *TriggerMu_phi, *TriggerMu_ch, *TriggerMu_ip;
      std::vector<int>         *mu_HLT_Mu7_IP4, *mu_HLT_Mu8_IP3, *mu_HLT_Mu8_IP5, *mu_HLT_Mu8_IP6, *mu_HLT_Mu8p5_IP3p5;
      std::vector<int>         *mu_HLT_Mu9_IP0, *mu_HLT_Mu9_IP3, *mu_HLT_Mu9_IP4, *mu_HLT_Mu9_IP5, *mu_HLT_Mu9_IP6, *mu_HLT_Mu10p5_IP3p5, *mu_HLT_Mu12_IP6;

      int verbose = 0;
};


ObjVsMuon_tupler::ObjVsMuon_tupler(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( iConfig.getParameter<edm::InputTag>("triggerBits") ) ),
  triggerPrescalesSrc_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),

  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),

  //Trigger Muon Selecctor 
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),

  TriggerObj_pt(0), TriggerObj_eta(0), TriggerObj_phi(0), TriggerObj_ch(0), TriggerObj_ip(0),
  obj_HLT_Mu7_IP4(0), obj_HLT_Mu8_IP3(0), obj_HLT_Mu8_IP5(0), obj_HLT_Mu8_IP6(0), obj_HLT_Mu8p5_IP3p5(0),
  obj_HLT_Mu9_IP0(0), obj_HLT_Mu9_IP3(0), obj_HLT_Mu9_IP4(0), obj_HLT_Mu9_IP5(0), obj_HLT_Mu9_IP6(0), obj_HLT_Mu10p5_IP3p5(0), obj_HLT_Mu12_IP6(0),
  
  TriggerMu_pt(0), TriggerMu_eta(0), TriggerMu_phi(0), TriggerMu_ch(0), TriggerMu_ip(0),
  mu_HLT_Mu7_IP4(0), mu_HLT_Mu8_IP3(0), mu_HLT_Mu8_IP5(0), mu_HLT_Mu8_IP6(0), mu_HLT_Mu8p5_IP3p5(0),
  mu_HLT_Mu9_IP0(0), mu_HLT_Mu9_IP3(0), mu_HLT_Mu9_IP4(0), mu_HLT_Mu9_IP5(0), mu_HLT_Mu9_IP6(0), mu_HLT_Mu10p5_IP3p5(0), mu_HLT_Mu12_IP6(0),

  verbose( iConfig.getParameter<int>( "verbose" ) )

{
  
  //tree = fs->make<TTree>( "T", "Trigger Objects and Trigger Muons TTree ");
}


// ------------ method called to for each event  ------------
void ObjVsMuon_tupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
//void ObjVsMuon_tupler::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescalesSrc_, triggerPrescales);

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];

  //Old used Trigger info
  edm::Handle<edm::TriggerResults> triggerResults_handle;
  iEvent.getByToken(triggerBitsSrc_, triggerResults_handle);

  //Trigger Object Selector
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  //Muons
  edm::Handle<edm::View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(muon_Label,thePATMuonHandle);

  //Beam spot handle & position (for displaced Tracks ...)
  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByToken(BSLabel_, theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
  math::XYZPoint referencePos(theBeamSpot->position());


  isRealData = iEvent.isRealData() ? 1 : 0 ;
  runNum     = iEvent.id().run();
  lumiNum    = iEvent.luminosityBlock();
  eventNum   = iEvent.id().event();

  //BPH trigger footprint
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  if (verbose) {cout << "\n ==== TRIGGER PATHS ==== " << endl;}
  //look and fill the BParked Trigger Objects ... 
  
  //emulate BParking MuonTriggerSelector 
  int int_obj = 0;
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames   
    int_obj++; 
    obj.unpackFilterLabels(iEvent, *triggerResults_handle);
    obj.unpackPathNames(names);
    bool isTriggerMuon = false;

    // checa que al menos un elemento sea un muon (ID=83)
    // que pasa con los demas?
    //if(debug) std::cout << "\tfilterIds size:   " << obj.filterIds().size()<< "\n";
    for (unsigned h = 0; h < obj.filterIds().size(); ++h)
    	if(obj.filterIds()[h] == 83){ 
        	isTriggerMuon = true; 
         	// std::cout << "\t\tType IDs:   " << 83 <<"\n";  //83 = muon
        	//break;
      	} else {
         	// std::cout << "\t\tXXXXXXXXXX   Not Muon:   " << obj.filterIds()[h] <<"\n";
         	isTriggerMuon = false;
    }
    if(!isTriggerMuon) continue;
    // Ahora checa que dentro de los filterlabes en al menos uno
    // exista hltL3 y Park
    isTriggerMuon = false;
    //std::cout << "\tfilterLabels size:  " << obj.filterLabels().size()<< "\n";
	  std::string filterLabel = "";
	  std::string filterName_ = "";
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){   
        std::string filterName = obj.filterLabels()[h];
		
        //if (debug) std::cout << "\t\tfilterlabel:  " << h << filterName << "\n";
        if(filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
          isTriggerMuon = true;
			    filterLabel = filterName;
        }
    }
	  if(!isTriggerMuon) continue;

	  float obj_pt = obj.pt();
	  float obj_eta = obj.eta();
	  float obj_phi = obj.phi();
	  TriggerObj_pt->push_back(obj_pt);
	  TriggerObj_eta->push_back(obj_eta);
	  TriggerObj_phi->push_back(obj_phi);
	  float obj_ch = obj.charge();
	  TriggerObj_ch->push_back(obj_ch);
    TriggerObj_ip->push_back(0.0);

    if(verbose){ 
      std::cout << "\n\t\t\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      //Print trigger object collection and type
	    std::cout << "\t\t\tCollection: " << obj.collection() << std::endl;
	    std::cout << "\t\t\tFilter Label: " << filterLabel << std::endl;
		  std::cout << "\t\t\tFilter Name: " << filterName_ << std::endl;
    }

  }//trigger objects
  int obj_HLT_Mu7_IP4_ = 0, obj_HLT_Mu8_IP3_ = 0, obj_HLT_Mu8_IP5_ = 0, obj_HLT_Mu8_IP6_ = 0, obj_HLT_Mu8p5_IP3p5_ = 0;
  int obj_HLT_Mu9_IP0_ = 0, obj_HLT_Mu9_IP3_ = 0, obj_HLT_Mu9_IP4_ = 0, obj_HLT_Mu9_IP5_ = 0, obj_HLT_Mu9_IP6_ = 0, obj_HLT_Mu10p5_IP3p5_ = 0, obj_HLT_Mu12_IP6_ = 0;

  if ( triggerResults_handle.isValid()) {
    const edm::TriggerNames &TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
    regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9]_v[0-9]");

    for (unsigned int i = 0, n = triggerResults_handle->size(); i < n; ++i) {
      auto trgName = TheTriggerNames.triggerName(i);
      if(!regex_match(trgName, txt_regex_path)) continue;
      if (!(triggerPrescales->getPrescaleForIndex(i))) continue;
      if (verbose) {
        cout << "Trigger " << trgName << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << endl;
      }
      	//Trigger save non-sense 
				if (trgName.find("HLT_Mu7_IP4") != std::string::npos){
					obj_HLT_Mu7_IP4_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}
				if (trgName.find("HLT_Mu8_IP3") != std::string::npos){
					obj_HLT_Mu8_IP3_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}
				if (trgName.find("HLT_Mu8_IP5") != std::string::npos){
					obj_HLT_Mu8_IP5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}				
				if (trgName.find("HLT_Mu8_IP6") != std::string::npos){
					obj_HLT_Mu8_IP6_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}
				if (trgName.find("HLT_Mu8p5_IP3p5") != std::string::npos){
					obj_HLT_Mu8p5_IP3p5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}

				if (trgName.find("HLT_Mu9_IP0") != std::string::npos){
					obj_HLT_Mu9_IP0_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}					
				if (trgName.find("HLT_Mu9_IP3") != std::string::npos){
					obj_HLT_Mu9_IP3_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}				
				if (trgName.find("HLT_Mu9_IP4") != std::string::npos){
					obj_HLT_Mu9_IP4_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}	
				if (trgName.find("HLT_Mu9_IP5") != std::string::npos){
					obj_HLT_Mu9_IP5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}	
				if (trgName.find("HLT_Mu9_IP6") != std::string::npos){
					obj_HLT_Mu9_IP6_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}	
				if (trgName.find("HLT_Mu10p5_IP3p5") != std::string::npos){
					obj_HLT_Mu10p5_IP3p5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}	
				if (trgName.find("HLT_Mu12_IP6") != std::string::npos){
					obj_HLT_Mu12_IP6_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << trgName << std::endl;
				}	

    }
  }//end trigger results 
  
  obj_HLT_Mu7_IP4->push_back(obj_HLT_Mu7_IP4_); 
	obj_HLT_Mu8_IP3->push_back(obj_HLT_Mu8_IP3_); 
	obj_HLT_Mu8_IP5->push_back(obj_HLT_Mu8_IP5_);
	obj_HLT_Mu8_IP6->push_back(obj_HLT_Mu8_IP6_);
	obj_HLT_Mu8p5_IP3p5->push_back(obj_HLT_Mu8p5_IP3p5_);

  obj_HLT_Mu9_IP0->push_back(obj_HLT_Mu9_IP0_); 
	obj_HLT_Mu9_IP3->push_back(obj_HLT_Mu9_IP3_);
	obj_HLT_Mu9_IP4->push_back(obj_HLT_Mu9_IP4_); 
	obj_HLT_Mu9_IP5->push_back(obj_HLT_Mu9_IP5_); 
	obj_HLT_Mu9_IP6->push_back(obj_HLT_Mu9_IP6_); 
	obj_HLT_Mu10p5_IP3p5->push_back(obj_HLT_Mu10p5_IP3p5_); 
	obj_HLT_Mu12_IP6->push_back(obj_HLT_Mu12_IP6_);
  


  //Trigger muons . . .
  int mu_HLT_Mu7_IP4_ = 0, mu_HLT_Mu8_IP3_ = 0, mu_HLT_Mu8_IP5_ = 0, mu_HLT_Mu8_IP6_ = 0, mu_HLT_Mu8p5_IP3p5_ = 0;
  int mu_HLT_Mu9_IP0_ = 0, mu_HLT_Mu9_IP3_ = 0, mu_HLT_Mu9_IP4_ = 0, mu_HLT_Mu9_IP5_ = 0, mu_HLT_Mu9_IP6_ = 0, mu_HLT_Mu10p5_IP3p5_ = 0, mu_HLT_Mu12_IP6_ = 0;
  for(auto tag : triggerTags) {
    std::cout << "djfoeirwnfir: " << tag << std::endl;
  }
  for(edm::View<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); ++iMuon) {
    for(auto tag : triggerTags) {
      std::string triggerName = tag + "_part*_v*"; 
      std::cout << "checking trigger :" << triggerName << std::endl;
			if(iMuon->triggered(triggerName.c_str())){ 
        std::cout << "Muon Fired !! " << triggerName << std::endl;
        float mu_pt = iMuon->pt();
	      float mu_eta = iMuon->eta();
	      float mu_phi = iMuon->phi();
	      TriggerMu_pt->push_back(mu_pt);
	      TriggerMu_eta->push_back(mu_eta);
	      TriggerMu_phi->push_back(mu_phi);
	      float mu_ch = iMuon->charge();
	      TriggerMu_ch->push_back(mu_ch);
        float mu_ip = iMuon->muonBestTrack()->dxy(referencePos)/iMuon->muonBestTrack()->dxyError();
        TriggerMu_ip->push_back(mu_ip);
        
        if (triggerName.find("HLT_Mu7_IP4") != std::string::npos){
					mu_HLT_Mu7_IP4_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}
				if (triggerName.find("HLT_Mu8_IP3") != std::string::npos){
					mu_HLT_Mu8_IP3_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}
				if (triggerName.find("HLT_Mu8_IP5") != std::string::npos){
					mu_HLT_Mu8_IP5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}				
				if (triggerName.find("HLT_Mu8_IP6") != std::string::npos){
					mu_HLT_Mu8_IP6_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}
				if (triggerName.find("HLT_Mu8p5_IP3p5") != std::string::npos){
					mu_HLT_Mu8p5_IP3p5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}

				if (triggerName.find("HLT_Mu9_IP0") != std::string::npos){
					mu_HLT_Mu9_IP0_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}					
				if (triggerName.find("HLT_Mu9_IP3") != std::string::npos){
					mu_HLT_Mu9_IP3_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}				
				if (triggerName.find("HLT_Mu9_IP4") != std::string::npos){
					mu_HLT_Mu9_IP4_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}	
				if (triggerName.find("HLT_Mu9_IP5") != std::string::npos){
					mu_HLT_Mu9_IP5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}	
				if (triggerName.find("HLT_Mu9_IP6") != std::string::npos){
					mu_HLT_Mu9_IP6_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}	
				if (triggerName.find("HLT_Mu10p5_IP3p5") != std::string::npos){
					mu_HLT_Mu10p5_IP3p5_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}	
				if (triggerName.find("HLT_Mu12_IP6") != std::string::npos){
					mu_HLT_Mu12_IP6_ = 1;
					if (verbose > 1) std::cout<< "-Mu 1" << triggerName << std::endl;
				}

        mu_HLT_Mu7_IP4->push_back(mu_HLT_Mu7_IP4_); 
	      mu_HLT_Mu8_IP3->push_back(mu_HLT_Mu8_IP3_); 
	      mu_HLT_Mu8_IP5->push_back(mu_HLT_Mu8_IP5_);
	      mu_HLT_Mu8_IP6->push_back(mu_HLT_Mu8_IP6_);
	      mu_HLT_Mu8p5_IP3p5->push_back(mu_HLT_Mu8p5_IP3p5_);

        mu_HLT_Mu9_IP0->push_back(mu_HLT_Mu9_IP0_); 
	      mu_HLT_Mu9_IP3->push_back(mu_HLT_Mu9_IP3_);
	      mu_HLT_Mu9_IP4->push_back(mu_HLT_Mu9_IP4_); 
	      mu_HLT_Mu9_IP5->push_back(mu_HLT_Mu9_IP5_); 
	      mu_HLT_Mu9_IP6->push_back(mu_HLT_Mu9_IP6_); 
	      mu_HLT_Mu10p5_IP3p5->push_back(mu_HLT_Mu10p5_IP3p5_); 
	      mu_HLT_Mu12_IP6->push_back(mu_HLT_Mu12_IP6_);	
      }// end if triggered
    }
  }	
  

  tree->Fill();

  TriggerObj_pt->clear();  TriggerObj_eta->clear(); TriggerObj_phi->clear(); TriggerObj_ch->clear(); TriggerObj_ip->clear(); 
  obj_HLT_Mu7_IP4->clear(); obj_HLT_Mu8_IP3->clear(); obj_HLT_Mu8_IP5->clear(); obj_HLT_Mu8_IP6->clear(); obj_HLT_Mu8p5_IP3p5->clear(); 
  obj_HLT_Mu9_IP0->clear(); obj_HLT_Mu9_IP3->clear(); obj_HLT_Mu9_IP4->clear(); obj_HLT_Mu9_IP5->clear(); obj_HLT_Mu9_IP6->clear(); obj_HLT_Mu10p5_IP3p5->clear(); obj_HLT_Mu12_IP6->clear(); 

  TriggerMu_pt->clear();  TriggerMu_eta->clear(); TriggerMu_phi->clear(); TriggerMu_ch->clear(); TriggerMu_ip->clear(); 
  mu_HLT_Mu7_IP4->clear(); mu_HLT_Mu8_IP3->clear(); mu_HLT_Mu8_IP5->clear(); mu_HLT_Mu8_IP6->clear(); mu_HLT_Mu8p5_IP3p5->clear(); 
  mu_HLT_Mu9_IP0->clear(); mu_HLT_Mu9_IP3->clear(); mu_HLT_Mu9_IP4->clear(); mu_HLT_Mu9_IP5->clear(); mu_HLT_Mu9_IP6->clear(); mu_HLT_Mu10p5_IP3p5->clear(); mu_HLT_Mu12_IP6->clear(); 
  //if (verbose) {cout << "======================== " << endl;}
  return;

}

 
void ObjVsMuon_tupler::beginJob()
{

  std::cout << "Beginning analyzer " << std::endl;

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("T", "Trigger Objects and Trigger Muons TTree ");
  if(verbose) {cout << "\nCreating the branches in the output tree:\n";}

  tree->Branch("isRealData", &isRealData);
  tree->Branch("runNum", &runNum);
  tree->Branch("lumiNum", &lumiNum);
  tree->Branch("eventNum", &eventNum);

  //Trigger objects . . . 
  tree->Branch("TriggerObj_pt", &TriggerObj_pt);
  tree->Branch("TriggerObj_eta", &TriggerObj_eta);
  tree->Branch("TriggerObj_phi", &TriggerObj_phi);
  tree->Branch("TriggerObj_ch", &TriggerObj_ch);
  tree->Branch("TriggerObj_ip", &TriggerObj_ip);
  tree->Branch("TriggerObj_pt", &TriggerObj_pt);

  tree->Branch("obj_HLT_Mu7_IP4", &obj_HLT_Mu7_IP4);
  tree->Branch("obj_HLT_Mu8_IP3", &obj_HLT_Mu8_IP3);
  tree->Branch("obj_HLT_Mu8_IP5", &obj_HLT_Mu8_IP5);
  tree->Branch("obj_HLT_Mu8_IP6", &obj_HLT_Mu8_IP6);
  tree->Branch("obj_HLT_Mu8p5_IP3p5", &obj_HLT_Mu8p5_IP3p5);

  tree->Branch("obj_HLT_Mu9_IP0", &obj_HLT_Mu9_IP0);
  tree->Branch("obj_HLT_Mu9_IP3", &obj_HLT_Mu9_IP3);
  tree->Branch("obj_HLT_Mu9_IP4", &obj_HLT_Mu9_IP4);
  tree->Branch("obj_HLT_Mu9_IP5", &obj_HLT_Mu9_IP5);
  tree->Branch("obj_HLT_Mu9_IP6", &obj_HLT_Mu9_IP6);
  tree->Branch("obj_HLT_Mu10p5_IP3p5", &obj_HLT_Mu10p5_IP3p5);
  tree->Branch("obj_HLT_Mu12_IP6", &obj_HLT_Mu12_IP6);

  //Trigger muons  . . . 
  tree->Branch("TriggerMu_pt",  &TriggerMu_pt);
  tree->Branch("TriggerMu_eta", &TriggerMu_eta);
  tree->Branch("TriggerMu_phi", &TriggerMu_phi);
  tree->Branch("TriggerMu_ch",  &TriggerMu_ch);
  tree->Branch("TriggerMu_ip",  &TriggerMu_ip);
  tree->Branch("TriggerMu_pt",  &TriggerMu_pt);

  tree->Branch("mu_HLT_Mu7_IP4", &mu_HLT_Mu7_IP4);
  tree->Branch("mu_HLT_Mu8_IP3", &mu_HLT_Mu8_IP3);
  tree->Branch("mu_HLT_Mu8_IP5", &mu_HLT_Mu8_IP5);
  tree->Branch("mu_HLT_Mu8_IP6", &mu_HLT_Mu8_IP6);
  tree->Branch("mu_HLT_Mu8p5_IP3p5", &mu_HLT_Mu8p5_IP3p5);

  tree->Branch("mu_HLT_Mu9_IP0", &mu_HLT_Mu9_IP0);
  tree->Branch("mu_HLT_Mu9_IP3", &mu_HLT_Mu9_IP3);
  tree->Branch("mu_HLT_Mu9_IP4", &mu_HLT_Mu9_IP4);
  tree->Branch("mu_HLT_Mu9_IP5", &mu_HLT_Mu9_IP5);
  tree->Branch("mu_HLT_Mu9_IP6", &mu_HLT_Mu9_IP6);
  tree->Branch("mu_HLT_Mu10p5_IP3p5", &mu_HLT_Mu10p5_IP3p5);
  tree->Branch("mu_HLT_Mu12_IP6", &mu_HLT_Mu12_IP6);
}


// ------------ method called once each job just after ending the event loop  ------------
void ObjVsMuon_tupler::endJob() {
  std::cout << "Ending job" << std::endl;
  tree->GetDirectory()->cd();
  tree->Write();
  std::cout << "Ending job ok " << std::endl;
}


DEFINE_FWK_MODULE(ObjVsMuon_tupler);