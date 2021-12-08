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
      double dR(double, double, double, double);
      RefCountedKinematicTree FitJpsi_mumu(const edm::EventSetup&, pat::Muon, pat::Muon, bool);
      struct kinFitResuts{
        bool isValid = false;
        double chi2 = -1;
        double dof = -1;
        double pval = -1;
        bool isGood = false;
      };
      kinFitResuts fitQuality(RefCountedKinematicTree, double = -1);

      explicit Prescale_Producer(const edm::ParameterSet&);
      ~Prescale_Producer() {
        cout << "Total events in output tree: " << tree->GetEntries() << endl;
      };

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      tuple<uint, float, float, float> matchL1Muon(pat::Muon muReco, BXVector<l1t::Muon> muonsL1, uint skipIdx=9999);
      void addToTree();

      // ----------member data ---------------------------
      edm::EDGetTokenT<BXVector<GlobalAlgBlk>> L1triggerBitsSrc_;
      edm::EDGetTokenT<edm::TriggerResults> L1triggerResultsSrc_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<BXVector<l1t::Muon>> l1MuonSrc_;
      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

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
          //if(verbose) {cout << tag << "\t" << n << "\t" << triggerBits->wasrun(i) << "\t" << triggerPrescales->getPrescaleForIndex(i) << endl;}
        }
      }
  }
  addToTree();
  if (verbose) {cout << "======================== " << endl;}
  return true;
}

tuple<uint, float, float, float> Prescale_Producer::matchL1Muon(pat::Muon muReco, BXVector<l1t::Muon> muonsL1, uint skipIdx) {
  uint idxMatch = 9999;
  float best_dR = 1e6;
  float best_dpt = 1e6;
  float best_pt = -1;

  // Approximately radius 2.5 m in a 3.8T magnetic field
  float dphi = 2 * TMath::ASin(2.5 * 0.3 * 3.8 / (2 * muReco.pt()) );
  float phiProp = muReco.phi() + TMath::Sign(1, muReco.pdgId()) * dphi;

  for (uint i=0; i < muonsL1.size(0); i++) {
    if (i == skipIdx) continue;
    auto m = muonsL1.at(0,i);
    if (m.hwQual() < 12) continue;
    float dR_ = dR(m.phi(), phiProp, m.eta(), muReco.eta());
    float dpt = fabs(muReco.pt() - m.pt())/muReco.pt();
    if ((dR_ < best_dR && dpt < best_dpt) || (dpt + dR_ < best_dpt + best_dR)) {
      best_dR = dR_;
      best_dpt = dpt;
      idxMatch = i;
      best_pt = m.pt();
    }
  }

  tuple<uint, float, float, float> out(idxMatch, best_dR, best_dpt, best_pt);
  return out;
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
double Prescale_Producer::dR(double p1, double p2, double e1, double e2){
    double dp = std::abs(p1 - p2);
    if (dp > (M_PI)){
      dp -= (2 * M_PI);
    }  
    return std::sqrt( (e1 - e2) * (e1 - e2) + dp * dp );

}
RefCountedKinematicTree Prescale_Producer::FitJpsi_mumu(const edm::EventSetup& iSetup, pat::Muon m1, pat::Muon m2, bool mass_constraint) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack m1_tk = TTBuilder->build(m1.muonBestTrack());
  reco::TransientTrack m2_tk = TTBuilder->build(m2.muonBestTrack());

  std::vector<RefCountedKinematicParticle> parts;
  KinematicParticleFactoryFromTransientTrack pFactory;
  double chi = 0, ndf = 0;
  float mMu = 0.1056583745, dmMu = 0.0000000024;
  parts.push_back(pFactory.particle(m1_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(m2_tk, mMu, chi, ndf, dmMu));

  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
    return KinTree;
  }
  else {
    ParticleMass mass = 3.096916;
    MultiTrackKinematicConstraint * mass_c = new TwoTrackMassKinematicConstraint(mass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree KinTree = kcVtxFitter.fit(parts, mass_c);
    return KinTree;
  }
}

Prescale_Producer::kinFitResuts Prescale_Producer::fitQuality(RefCountedKinematicTree t, double pval_thr){
  kinFitResuts out;
  if(t->isValid()) {
    out.isValid = true;
    t->movePointerToTheTop();
    out.chi2 = t->currentDecayVertex()->chiSquared();
    out.dof = t->currentDecayVertex()->degreesOfFreedom();
    out.pval = ChiSquaredProbability(out.chi2, out.dof);
    if (pval_thr > 0) {
      out.isGood = out.pval > pval_thr;
    }
  }
  return out;
}
DEFINE_FWK_MODULE(Prescale_Producer);