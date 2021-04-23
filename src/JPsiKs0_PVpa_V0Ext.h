#ifndef _JPsiKs0_PVpa_V0Ext_h
#define _JPsiKs0_PVpa_V0Ext_h

// system include files
#include <memory>

// user include files
//#include "myAnalyzers/BtoKsMuMu/interface/JPsif0PAT.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

//
// class decleration
//

class JPsiKs0_PVpa_V0Ext : public edm::EDAnalyzer {
public:
  explicit JPsiKs0_PVpa_V0Ext(const edm::ParameterSet&);
  ~JPsiKs0_PVpa_V0Ext();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const reco::Track& tk, const pat::Muon& mu);
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  bool IsTheSame(const pat::GenericParticle& tk1, const pat::GenericParticle& tk2);
  bool isAncestor(const reco::Candidate*, const reco::Candidate*);
  bool isAncestor(int, const reco::Candidate*);
  double GetLifetime(TLorentzVector, TVector3, TVector3);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  void printMCtree(const reco::Candidate *, int);
  void printMCtreeUP(const reco::Candidate *, int);
  std::string printName(int);
 
  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> v0PtrCollection_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> Lost_track_label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pPFC_track_label;

  //Trigger Muon Selector
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;

  bool OnlyBest_;
  bool isMC_;
  bool isRes_;
  bool OnlyGen_;
  bool doMC_;

  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;
  
  std::vector<float>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;
  std::vector<float>       *trg_dzm1, *trg_dzm2;

  std::vector<int>         *tri_Dim25, *tri_JpsiTk, *tri_JpsiTkTk;
 
  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  

  //Triger Selector
  std::vector<float>       *drTrg_m1, *drTrg_m2; 
 
  int                      muAcc, muTrig, weight;
 
  // vertice primario CON mejor pointing angle
  unsigned int             nVtx;
  std::vector<unsigned int> *nTks;
  std::vector<int>         *TrkIndex;
  std::vector<float>       *PVTriggDz;
  std::vector<float>       *priVtxX, *priVtxY, *priVtxZ, *priVtxXE, *priVtxYE, *priVtxZE, *priVtxCL;
  std::vector<float>       *priVtxXYE, *priVtxXZE, *priVtxYZE;
 
  // Track Container used for B candidate ...
  std::vector<int>         *trackContainer; // 1 = V0 Container, 2 = lost Tracks , 3 = Packed PF Candidates

  // ********************************** ************************************************************************

  std::vector<float>       *bDecayVtxX, *bDecayVtxY, *bDecayVtxZ;
  std::vector<double>      *bDecayVtxXE, *bDecayVtxYE, *bDecayVtxZE;
  std::vector<double>      *bDecayVtxXYE, *bDecayVtxXZE, *bDecayVtxYZE;

  std::vector<float>       *VDecayVtxX, *VDecayVtxY, *VDecayVtxZ;
  std::vector<float>       *VDecayVtxXE, *VDecayVtxYE, *VDecayVtxZE;
  std::vector<float>       *VDecayVtxXYE, *VDecayVtxXZE, *VDecayVtxYZE;

  // *************************************

  unsigned int             nB;
  unsigned int             nMu;

  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz;
 
  std::vector<float>       *B_Ks0_mass, *B_Ks0_px, *B_Ks0_py, *B_Ks0_pz;
  std::vector<float>       *B_Ks0_pt1, *B_Ks0_px1, *B_Ks0_py1, *B_Ks0_pz1;
  std::vector<float>       *B_Ks0_pt2, *B_Ks0_px2, *B_Ks0_py2, *B_Ks0_pz2;
  
  std::vector<float>       *B_Ks0_px1_track, *B_Ks0_py1_track, *B_Ks0_pz1_track;
  std::vector<float>       *B_Ks0_px2_track, *B_Ks0_py2_track, *B_Ks0_pz2_track;

  std::vector<float>       *pi1dxy, *pi2dxy, *pi1dz, *pi2dz;
  std::vector<float>       *pi1dxy_e, *pi2dxy_e, *pi1dz_e, *pi2dz_e;
  std::vector<int>         *B_Ks0_charge1, *B_Ks0_charge2;
  
  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;
  std::vector<float>       *B_J_pt1, *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_pt2, *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;
  
  std::vector<float>       *B_Ks0_chi2, *B_J_chi2, *B_chi2;
  std::vector<float>       *B_Prob, *B_J_Prob, *B_ks0_Prob;

  int  run, event;
  int  lumiblock;
  UInt_t trigger;
  
  TLorentzVector gen_b_p4,gen_jpsi_p4,gen_pion1_p4,gen_pion2_p4,gen_ks0_p4,gen_muon1_p4,gen_muon2_p4;
  TVector3       gen_b_vtx,gen_jpsi_vtx,gen_ks0_vtx;
  float          gen_b_ct, gen_ks0_ct;
  int ngen;

};

#endif
