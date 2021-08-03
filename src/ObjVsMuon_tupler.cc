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

      edm::Service<TFileService> fs;

      TTree* tree;
      map<string, float> outMap;
      bool treeDeclared = false;

      bool isRealData;
      unsigned int runNum;
      unsigned int lumiNum;
      unsigned long long eventNum;

      vector<string> triggerTag = {"Mu12_IP6", 
	                                 "Mu10p5_IP3p5", "Mu9_IP0", "Mu9_IP3", "Mu9_IP4", "Mu9_IP5", "Mu9_IP6", 
                                   "Mu8p5_IP3p5",  "Mu8_IP3", "Mu8_IP5", "Mu8_IP6", 
                                   "Mu7_IP4"};
      int verbose = 0;
};


ObjVsMuon_tupler::ObjVsMuon_tupler(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( iConfig.getParameter<edm::InputTag>("triggerBits") ) ),
  triggerPrescalesSrc_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),

  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  tree = fs->make<TTree>( "T", "Events Tree from N Vtx Produccer");
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

  isRealData = iEvent.isRealData() ? 1 : 0 ;
  runNum     = iEvent.id().run();
  lumiNum    = iEvent.luminosityBlock();
  eventNum   = iEvent.id().event();

  //BPH trigger footprint
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9].*");
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  if (verbose) {cout << "\n ==== TRIGGER PATHS ==== " << endl;}

  //for (auto trgTag : triggerTags) outMap["prescale_" + trgTag] = 0;

  //look and fill the BParked Trigger Objects ... 
  if ( triggerResults_handle.isValid()) {
  const edm::TriggerNames &TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);

    for (unsigned int i = 0, n = triggerResults_handle->size(); i < n; ++i) {
      auto trgName = TheTriggerNames.triggerName(i);
      //if (verbose) {
      // cout << "Trigger " << trgName << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << endl;
      //}

      for (unsigned int k = 0; k < NTRIGGERS; k++){
        std::string trgTag = TriggersToTest[k];	
	      bool found_ = false; 
        bool match = (trgName.find(trgTag) != std::string::npos && !found_);
        //bool match = trgName.substr(4, trgTag.size()) == trgTag.c_str();
        if (triggerPrescales->getPrescaleForIndex(i) < 1) continue;
        if (match && triggerResults_handle->accept(i)) {
		      trigger += (1<<k);
		      found_ = true;
	      }
      }
    }
  }

  
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