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

using namespace std;

class NumberOfVertexesProducer : public edm::EDProducer {
   public:
      explicit NumberOfVertexesProducer(const edm::ParameterSet& iConfig);
      ~NumberOfVertexesProducer() {};

   private:
      void beginJob(const edm::EventSetup&) {};
      void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;

      edm::Service<TFileService> fs;
      //map<string, TH1D*> hNvtx;
      //map<string, TH1D*> hNvtxPassed;
      //map<string, TH1D*> hZvtxPassed;

      TH1I* hAllNvts_Mu12_IP6;
      TH1I* hAllNTrueIntMC_Mu12_IP6;
      TH1I* hAllVtxZ_Mu12_IP6;

      
      TH1I* hAllNvts_Mu9_IP6;
      TH1I* hAllNTrueIntMC_Mu9_IP6;
      TH1I* hAllVtxZ_Mu9_IP6;

      TH1I* hAllNvts_Mu9_IP5;
      TH1I* hAllNTrueIntMC_Mu9_IP5;
      TH1I* hAllVtxZ_Mu9_IP5;

      TH1I* hAllNvts_Mu9_IP4;
      TH1I* hAllNTrueIntMC_Mu9_IP4;
      TH1I* hAllVtxZ_Mu9_IP4;

      TH1I* hAllNvts_Mu9_IP3;
      TH1I* hAllNTrueIntMC_Mu9_IP3;
      TH1I* hAllVtxZ_Mu9_IP3;

      TH1I* hAllNvts_Mu9_IP0;
      TH1I* hAllNTrueIntMC_Mu9_IP0;
      TH1I* hAllVtxZ_Mu9_IP0;


      TH1I* hAllNvts_Mu8_IP6;
      TH1I* hAllNTrueIntMC_Mu8_IP6;
      TH1I* hAllVtxZ_Mu8_IP6;

      TH1I* hAllNvts_Mu8_IP5;
      TH1I* hAllNTrueIntMC_Mu8_IP5;
      TH1I* hAllVtxZ_Mu8_IP5;

      TH1I* hAllNvts_Mu8_IP3;
      TH1I* hAllNTrueIntMC_Mu8_IP3;
      TH1I* hAllVtxZ_Mu8_IP3;


      TH1I* hAllNvts_Mu7_IP4;
      TH1I* hAllNTrueIntMC_Mu7_IP4;
      TH1I* hAllVtxZ_Mu7_IP4;

      TTree* tree;
      map<string, float> outMap;

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
  hAllNTrueIntMC_Mu12_IP6 = fs->make<TH1I>("hAllNTrueIntMC_Mu12_IP6", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu12_IP6 = fs->make<TH1I>("hAllVtxZ_Mu12_IP6", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);


  hAllNvts_Mu9_IP6 = fs->make<TH1I>("hAllNvts_Mu9_IP6", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu9_IP6 = fs->make<TH1I>("hAllNTrueIntMC_Mu9_IP6", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP6 = fs->make<TH1I>("hAllVtxZ_Mu9_IP6", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP5 = fs->make<TH1I>("hAllNvts_Mu9_IP5", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu9_IP5 = fs->make<TH1I>("hAllNTrueIntMC_Mu9_IP5", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP5 = fs->make<TH1I>("hAllVtxZ_Mu9_IP5", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP4 = fs->make<TH1I>("hAllNvts_Mu9_IP4", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu9_IP4 = fs->make<TH1I>("hAllNTrueIntMC_Mu9_IP4", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP4 = fs->make<TH1I>("hAllVtxZ_Mu9_IP4", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP3 = fs->make<TH1I>("hAllNvts_Mu9_IP3", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu9_IP3 = fs->make<TH1I>("hAllNTrueIntMC_Mu9_IP3", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP3 = fs->make<TH1I>("hAllVtxZ_Mu9_IP3", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu9_IP0 = fs->make<TH1I>("hAllNvts_Mu9_IP0", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu9_IP0 = fs->make<TH1I>("hAllNTrueIntMC_Mu9_IP0", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu9_IP0 = fs->make<TH1I>("hAllVtxZ_Mu9_IP0", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);


  hAllNvts_Mu8_IP6 = fs->make<TH1I>("hAllNvts_Mu8_IP6", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu8_IP6 = fs->make<TH1I>("hAllNTrueIntMC_Mu8_IP6", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu8_IP6 = fs->make<TH1I>("hAllVtxZ_Mu8_IP6", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu8_IP5 = fs->make<TH1I>("hAllNvts_Mu8_IP5", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu8_IP5 = fs->make<TH1I>("hAllNTrueIntMC_Mu8_IP5", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu8_IP5 = fs->make<TH1I>("hAllVtxZ_Mu8_IP5", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);

  hAllNvts_Mu8_IP3 = fs->make<TH1I>("hAllNvts_Mu8_IP3", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu8_IP3 = fs->make<TH1I>("hAllNTrueIntMC_Mu8_IP3", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ_Mu8_IP3 = fs->make<TH1I>("hAllVtxZ_Mu8_IP3", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);


  hAllNvts_Mu7_IP4 = fs->make<TH1I>("hAllNvts_Mu7_IP4", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC_Mu7_IP4 = fs->make<TH1I>("hAllNTrueIntMC_Mu7_IP4", "Number of true interactions generated in MC", 101, -0.5, 100.5);
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

  //BPH trigger footprint
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9].*");
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  if (verbose) {cout << "\n ==== TRIGGER PATHS ==== " << endl;}

  //for (auto trgTag : triggerTags){(*outputNtuplizer)["prescale_" + trgTag] = 0;}

}

DEFINE_FWK_MODULE(NumberOfVertexesProducer);