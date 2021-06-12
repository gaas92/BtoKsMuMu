// Taken from https://github.com/ocerri/BPH_RDntuplizer/blob/master/plugins/NumberOfVertexesProducer.cc
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

class NumberOfVertexesProducer : public edm::EDProducer {
   public:
      explicit NumberOfVertexesProducer(const edm::ParameterSet& iConfig);
      ~NumberOfVertexesProducer() {};

   private:
      void beginJob(const edm::EventSetup&) {};
      void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      void addToTree();

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;

      edm::Service<TFileService> fs;
      //map<string, TH1D*> hNvtx;
      //map<string, TH1D*> hNvtxPassed;
      //map<string, TH1D*> hZvtxPassed;

      TH1I* hAllNvts_Mu12_IP6;
      TH1I* hNvtxPassed_Mu12_IP6;
      TH1I* hAllVtxZ_Mu12_IP6;

      
      TH1I* hAllNvts_Mu9_IP6;
      TH1I* hNvtxPassed_Mu9_IP6;
      TH1I* hAllVtxZ_Mu9_IP6;

      TH1I* hAllNvts_Mu9_IP5;
      TH1I* hNvtxPassed_Mu9_IP5;
      TH1I* hAllVtxZ_Mu9_IP5;

      TH1I* hAllNvts_Mu9_IP4;
      TH1I* hNvtxPassed_Mu9_IP4;
      TH1I* hAllVtxZ_Mu9_IP4;

      TH1I* hAllNvts_Mu9_IP3;
      TH1I* hNvtxPassed_Mu9_IP3;
      TH1I* hAllVtxZ_Mu9_IP3;

      TH1I* hAllNvts_Mu9_IP0;
      TH1I* hNvtxPassed_Mu9_IP0;
      TH1I* hAllVtxZ_Mu9_IP0;


      TH1I* hAllNvts_Mu8_IP6;
      TH1I* hNvtxPassed_Mu8_IP6;
      TH1I* hAllVtxZ_Mu8_IP6;

      TH1I* hAllNvts_Mu8_IP5;
      TH1I* hNvtxPassed_Mu8_IP5;
      TH1I* hAllVtxZ_Mu8_IP5;

      TH1I* hAllNvts_Mu8_IP3;
      TH1I* hNvtxPassed_Mu8_IP3;
      TH1I* hAllVtxZ_Mu8_IP3;


      TH1I* hAllNvts_Mu7_IP4;
      TH1I* hNvtxPassed_Mu7_IP4;
      TH1I* hAllVtxZ_Mu7_IP4;

      TTree* tree;
      map<string, float> outMap;
      bool treeDeclared = false;

      bool isRealData;
      unsigned int runNum;
      unsigned int lumiNum;
      unsigned long long eventNum;

      vector<string> triggerTags = {"Mu12_IP6",
                                    "Mu9_IP6", "Mu9_IP5", "Mu9_IP4", "Mu9_IP3", "Mu9_IP0",
                                    "Mu8_IP6", "Mu8_IP5", "Mu8_IP3",
                                    "Mu7_IP4"};
      int verbose = 0;
};


NumberOfVertexesProducer::NumberOfVertexesProducer(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( iConfig.getParameter<edm::InputTag>("triggerBits") ) ),
  triggerPrescalesSrc_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),

  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  hAllNvts_Mu12_IP6 = fs->make<TH1I>("hAllNvts_Mu12_IP6", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu12_IP6 = fs->make<TH1I>("hNvtxPassed_Mu12_IP6", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu12_IP6 = fs->make<TH1I>("hAllVtxZ_Mu12_IP6", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);


  hAllNvts_Mu9_IP6 = fs->make<TH1I>("hAllNvts_Mu9_IP6", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu9_IP6 = fs->make<TH1I>("hNvtxPassed_Mu9_IP6", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP6 = fs->make<TH1I>("hAllVtxZ_Mu9_IP6", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP5 = fs->make<TH1I>("hAllNvts_Mu9_IP5", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu9_IP5 = fs->make<TH1I>("hNvtxPassed_Mu9_IP5", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP5 = fs->make<TH1I>("hAllVtxZ_Mu9_IP5", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP4 = fs->make<TH1I>("hAllNvts_Mu9_IP4", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu9_IP4 = fs->make<TH1I>("hNvtxPassed_Mu9_IP4", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP4 = fs->make<TH1I>("hAllVtxZ_Mu9_IP4", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP3 = fs->make<TH1I>("hAllNvts_Mu9_IP3", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu9_IP3 = fs->make<TH1I>("hNvtxPassed_Mu9_IP3", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP3 = fs->make<TH1I>("hAllVtxZ_Mu9_IP3", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP0 = fs->make<TH1I>("hAllNvts_Mu9_IP0", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu9_IP0 = fs->make<TH1I>("hNvtxPassed_Mu9_IP0", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP0 = fs->make<TH1I>("hAllVtxZ_Mu9_IP0", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);


  hAllNvts_Mu8_IP6 = fs->make<TH1I>("hAllNvts_Mu8_IP6", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu8_IP6 = fs->make<TH1I>("hNvtxPassed_Mu8_IP6", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu8_IP6 = fs->make<TH1I>("hAllVtxZ_Mu8_IP6", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu8_IP5 = fs->make<TH1I>("hAllNvts_Mu8_IP5", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu8_IP5 = fs->make<TH1I>("hNvtxPassed_Mu8_IP5", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu8_IP5 = fs->make<TH1I>("hAllVtxZ_Mu8_IP5", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu8_IP3 = fs->make<TH1I>("hAllNvts_Mu8_IP3", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu8_IP3 = fs->make<TH1I>("hNvtxPassed_Mu8_IP3", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu8_IP3 = fs->make<TH1I>("hAllVtxZ_Mu8_IP3", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);


  hAllNvts_Mu7_IP4 = fs->make<TH1I>("hAllNvts_Mu7_IP4", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hNvtxPassed_Mu7_IP4 = fs->make<TH1I>("hNvtxPassed_Mu7_IP4", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu7_IP4 = fs->make<TH1I>("hAllVtxZ_Mu7_IP4", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  tree = fs->make<TTree>( "T", "Events Tree from N Vtx Produccer");
}

void NumberOfVertexesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescalesSrc_, triggerPrescales);

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];

  isRealData = iEvent.isRealData() ? 1 : 0 ;
  runNum     = iEvent.id().run();
  lumiNum    = iEvent.luminosityBlock();
  eventNum   = iEvent.id().event();

  //BPH trigger footprint
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9].*");
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  if (verbose) {cout << "\n ==== TRIGGER PATHS ==== " << endl;}

  for (auto trgTag : triggerTags) outMap["prescale_" + trgTag] = 0;

  map<string, bool> trgActive;
  map<string, bool> trgPassed;
  for (auto trgTag : triggerTags) {
    trgActive[trgTag] = false;
    trgPassed[trgTag] = false;
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
      if (match && triggerBits->accept(i)) trgPassed[trgTag] = true;

       //for (int part = 0; part <= 5; part++) {
       //  regex rule(Form("HLT_%s_part%d.*", trgTag.c_str(), part));
       //  if (regex_match(trgName, rule)) {
       //    outMap[Form("prescale_%s_part%d", trgTag.c_str(), part)] = triggerPrescales->getPrescaleForIndex(i);
       //    if (triggerPrescales->getPrescaleForIndex(i) > 0) outMap["prescale_" + trgTag]++;
       //  }
       //}

    }

  }

  auto Nvtx = vtxHandle->size();
  bool something_to_fill = false;
  if(trgActive["Mu12_IP6"]) {
    if (verbose) {cout << "Filling active " << "Mu12_IP6" << endl;}
    hAllNvts_Mu12_IP6->Fill(Nvtx);
    outMap["Nvtx_Mu12_IP6_A"] = Nvtx; 
    something_to_fill = true
  }
  if(trgPassed["Mu12_IP6"]) {
    if (verbose) {cout << "Filling passed " << "Mu12_IP6" << endl;}
    hNvtxPassed_Mu12_IP6->Fill(Nvtx);
    for(auto vtx : (*vtxHandle)) hAllVtxZ_Mu12_IP6->Fill(vtx.position().z());
    outMap["Nvtx_Mu12_IP6_P"] = Nvtx; 
    something_to_fill = true
  }
  
  if(trgActive["Mu9_IP6"]) {
    if (verbose) {cout << "Filling active " << "Mu9_IP6" << endl;}
    hAllNvts_Mu9_IP6->Fill(Nvtx);
    outMap["Nvtx_Mu9_IP6_A"] = Nvtx; 
    something_to_fill = true
  }
  if(trgPassed["Mu9_IP6"]) {
    if (verbose) {cout << "Filling passed " << "Mu9_IP6" << endl;}
    hNvtxPassed_Mu9_IP6->Fill(Nvtx);
    for(auto vtx : (*vtxHandle)) hAllVtxZ_Mu9_IP6->Fill(vtx.position().z());
    outMap["Nvtx_Mu9_IP6_P"] = Nvtx; 
    something_to_fill = true
  }

  outMap["N_vertexes"] = Nvtx;

  if(something_to_fill) addToTree();

  if (verbose) {cout << "======================== " << endl;}
  return;

}

void NumberOfVertexesProducer::addToTree() {
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
DEFINE_FWK_MODULE(NumberOfVertexesProducer);