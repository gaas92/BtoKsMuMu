// system include files
#include <memory>
#include <iostream>
#include <string>
#include <regex>
#include <math.h>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2D.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Needed for Transient Tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
// Needed for the kinematic fit
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
// Needed for Chi2 computation
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// L1 trigger
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;

class Prescale_Producer : public edm::stream::EDFilter<> {
   public:

      explicit Prescale_Producer(const edm::ParameterSet&);
      ~Prescale_Producer() {
        cout << "Total events in output tree: " << tree->GetEntries() << endl;
      };

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      void addToTree();

      // ----------member data ---------------------------

      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescales_;

      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

      edm::Service<TFileService> fs;


      TTree* tree;
      map<string, float> outMap;
      bool treeDeclared = false;

      bool isRealData;
      unsigned int runNum;
      unsigned int lumiNum;
      unsigned long long eventNum;
      int verbose = 0;


};


Prescale_Producer::Prescale_Producer(const edm::ParameterSet& iConfig):
  triggerBitsSrc_(consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","HLT") )),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>( edm::InputTag("patTrigger") )),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("slimmedPatTrigger"))),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  tree = fs->make<TTree>( "T", "Events Tree Prescales");
}

bool Prescale_Producer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  isRealData = iEvent.isRealData() ? 1 : 0 ;
  runNum     = iEvent.id().run();
  lumiNum    = iEvent.luminosityBlock();
  eventNum   = iEvent.id().event();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);


  if (verbose) {cout << "\n\n =================  Event " << eventNum << " =================== " << endl;}

  vector<string> triggerTag = {"Mu12_IP6", 
	                             "Mu10p5_IP3p5", "Mu9_IP0", "Mu9_IP3", "Mu9_IP4", "Mu9_IP5", "Mu9_IP6", 
                               "Mu8p5_IP3p5",  "Mu8_IP3", "Mu8_IP5", "Mu8_IP6", 
                               "Mu7_IP4"};

  for(auto tag : triggerTag) outMap["prescale" + tag] = -1;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9]_v[0-9]");
  if (verbose) { cout << "\n == TRIGGER PATHS == \n";}
  for (unsigned int i = 0; i < triggerBits->size(); ++i) {
      auto n =  names.triggerName(i);
      if (verbose && triggerBits->accept(i) && false) {
        cout << n << ": PASS" << endl;
      }
      if(!regex_match(n, txt_regex_path)) continue;
      for(auto tag : triggerTag) {
        if(n.substr(4, tag.size()) == tag) {
          outMap["prescale" + tag] = triggerPrescales->getPrescaleForIndex(i);
          if(verbose) {cout << tag << "\t" << n << "\t" << triggerBits->wasrun(i) << "\t" << triggerPrescales->getPrescaleForIndex(i) << endl;}
        }
      }
  }
  addToTree();
  if (verbose) {cout << "======================== " << endl;}
  return true;
}


void Prescale_Producer::addToTree() {
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

DEFINE_FWK_MODULE(Prescale_Producer);