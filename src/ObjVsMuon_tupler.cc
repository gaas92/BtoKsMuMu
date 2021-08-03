// Taken from https://github.com/ocerri/BPH_RDntuplizer/blob/master/plugins/ObjVsMuon_tupler.cc
// Original Author: Olmo Cerri 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
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

class ObjVsMuon_tupler : public edm::EDProducer {
   public:
      explicit ObjVsMuon_tupler(const edm::ParameterSet& iConfig);
      ~ObjVsMuon_tupler() {};

   private:
      void beginJob(const edm::EventSetup&) {};
      void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      void addToTree();

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;

      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

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

      int verbose = 0;
};


ObjVsMuon_tupler::ObjVsMuon_tupler(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( iConfig.getParameter<edm::InputTag>("triggerBits") ) ),
  triggerPrescalesSrc_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),

  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),

  //Trigger Muon Selecctor 
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),

  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  tree = fs->make<TTree>( "T", "Trigger Objects and Trigger Muons TTree ");
}

void ObjVsMuon_tupler::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  isRealData = iEvent.isRealData() ? 1 : 0 ;
  runNum     = iEvent.id().run();
  lumiNum    = iEvent.luminosityBlock();
  eventNum   = iEvent.id().event();

  //BPH trigger footprint
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9].*");
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
	  TriggerObj_eta->push_back(obj_phi);
	  TriggerObj_phi->push_back(obj_eta);
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


  //Trigger Muons ...

  
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    auto trgName = names.triggerName(i);
    if (!regex_match(trgName, txt_regex_path)) continue;
    if (verbose) {
      cout << "Trigger " << trgName << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << endl;
    }

    for (auto trgTag : triggerTags){
      bool match = trgName.substr(4, trgTag.size()) == trgTag.c_str();
      if (match && triggerPrescales->getPrescaleForIndex(i) > 0) trgActive[trgTag] = true;
      if (match && triggerBits->accept(i)){
        trgPassed[trgTag] = true;
        for(unsigned int k = 0; k < NTRIGGERS; i++){
		      if (trgTag.find(TriggersToTest[i]) != std::string::npos){
            new_trigger += (1<<k);
          }
        }
       //for (int part = 0; part <= 5; part++) {
       //  regex rule(Form("HLT_%s_part%d.*", trgTag.c_str(), part));
       //  if (regex_match(trgName, rule)) {
       //    outMap[Form("prescale_%s_part%d", trgTag.c_str(), part)] = triggerPrescales->getPrescaleForIndex(i);
       //    if (triggerPrescales->getPrescaleForIndex(i) > 0) outMap["prescale_" + trgTag]++;
       //  }
       //}

      }
    }
  }
  outMap["NewTrigger_Int"] = new_trigger;


  if(something_to_fill) addToTree();

  if (verbose) {cout << "======================== " << endl;}
  return;

}

void ObjVsMuon_tupler::addToTree() {
  if (!treeDeclared) {
    if(verbose) {cout << "\nCreating the branches in the output tree:\n";}
    tree->Branch("isRealData", &isRealData);
    tree->Branch("runNum", &runNum);
    tree->Branch("lumiNum", &lumiNum);
    tree->Branch("eventNum", &eventNum);

    for(auto& kv : outMap) {
      auto k = kv.first;
      if(verbose) {cout << "\t" << k;}
      tree->Branch(k.c_str(), &(outMap[k]));
    }
    treeDeclared = true;
    if(verbose) {cout << "\n\n";}
  }

  tree->Fill();
}
DEFINE_FWK_MODULE(ObjVsMuon_tupler);