// -*- C++ -*-
//
// Package:    JPsiLam0_PVpa
// Class:      JPsiLam0_PVpa
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  November 2020                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>

// user include files
#include "myAnalyzers/BtoKsMuMu/src/JPsiLam0_PVpa.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <vector>
#include <utility>
#include <string>
//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//
JPsiLam0_PVpa::JPsiLam0_PVpa(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secundaryVerticesPtr"))),	       
  //Trigger Muon Selecctor 
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),

  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))), 
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))), 
  
  
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  isRes_(iConfig.getParameter<bool>("isRes")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  trg_dzm1(0), trg_dzm2(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),
  //Trigger Selector
  drTrg_m1(0), drTrg_m2(0),

  nVtx(0), nTks(0), TrkIndex(0), PVTriggDz(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
 
  // ************************ ****************************************************

  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),
  
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0), VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 

  // *******************************************************
  nB(0), nMu(0),
  B_mass(0),
  LB_mass(0), B_px(0), B_py(0), B_pz(0),
  
  B_Ks0_mass(0),
  LB_L0_mass(0), B_Ks0_px(0), B_Ks0_py(0), B_Ks0_pz(0),
  B_Ks0_pt1(0), B_Ks0_px1(0), B_Ks0_py1(0), B_Ks0_pz1(0), 
  B_Ks0_pt2(0), B_Ks0_px2(0), B_Ks0_py2(0), B_Ks0_pz2(0), 

  B_Ks0_px1_track(0), B_Ks0_py1_track(0), B_Ks0_pz1_track(0), 
  B_Ks0_px2_track(0), B_Ks0_py2_track(0), B_Ks0_pz2_track(0), 

  pi1dxy(0), pi2dxy(0), pi1dz(0), pi2dz(0),
  pi1dxy_e(0), pi2dxy_e(0), pi1dz_e(0), pi2dz_e(0),
  B_Ks0_charge1(0), B_Ks0_charge2(0),

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0), 
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),

  B_Ks0_chi2(0), B_J_chi2(0), B_chi2(0),
  B_Prob(0), B_J_Prob(0), B_ks0_Prob(0),

  run(0), event(0),
  lumiblock(0),
  trigger(0)

{
   //now do what ever initialization is needed
}

JPsiLam0_PVpa::~JPsiLam0_PVpa()
{

}

// ------------ method called to for each event  ------------
void JPsiLam0_PVpa::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  //std::cout << "analyse ok" <<std::endl;
  //*********************************
  // Get event content information
  //*********************************  

  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> theV0PtrHandle;
  iEvent.getByToken(v0PtrCollection_,  theV0PtrHandle);

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  //Trigger info
  edm::Handle<edm::TriggerResults> triggerResults_handle;
  iEvent.getByToken(triggerResults_Label, triggerResults_handle);

  //Trigger Muon Selector
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  // Not so clear
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMPathsAndTriggerBits
  //edm::Handle<edm::TriggerResults> triggerBits;
  //iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults_handle);

  //*********************************
  // Get gen level information
  //*********************************
  edm::Handle<reco::GenParticleCollection> pruned;
  //edm::Handle<pat::PackedGenParticle> pruned; 
  iEvent.getByToken(genCands_, pruned);
  
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);
  //For simulated events, only a selected set of particles is stored because the simulated particle
  //format, called GenParticle, takes a lot of space. First, a set called pruned GenParticles that
  //includes initial partons, heavy flavor particles, electroweak bosons, and leptons is stored in full.
  //Second, a set called packed GenParticles that have only the four-momentum and particle type
  //is saved for particles representing the final state particles in the event. Generated jets and some
  //reference information is also saved.
  //With the packed GenParticles, an analyst can re-make the generated jets with various
  //algorithms. The pruned GenParticles enable event classification, flavor definition, and matching
  //to reconstructed physics objects. Links from each packed GenParticle to its last surviving
  //ancestor pruned GenParticle allow the decay chain of the event to be reconstructed.

  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_ks0_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_ks0_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -9999.;
  gen_ks0_ct = -9999.;
  // for resonant in jpsi 
  //MONTE CARLO NEEDS TO BE CHANGED FOR LAMBDA RESONANCE 
  //Solo encontre MC resoante /LambdaBToJpsiLambda_JpsiToMuMu_TuneCP5_13TeV-pythia8-evtgen/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM 
  if ( (isMC_ || OnlyGen_) && pruned.isValid() && isRes_) {   
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 5122) ) { //&& (dau->status() == 2) ) { //Lambda b 
	    foundit++;
		std::cout<< "Lambda B Found" << std::endl;
	    gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	    gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
	    for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	      const reco::Candidate *gdau = dau->daughter(k);
	      if (gdau->pdgId()==443 ) { // if Jpsi
	        foundit++;
	        gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	        gen_b_ct = GetLifetime(gen_b_p4,gen_b_vtx,gen_jpsi_vtx);     
	        int nm=0;
	        for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
	          const reco::Candidate *mm = gdau->daughter(l);
	          if (mm->pdgId()==13) { foundit++;
	    		if (mm->status()!=1) {
	    		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
	    		    const reco::Candidate *mu = mm->daughter(m);
	    		    if (mu->pdgId()==13 ) { //&& mu->status()==1) {
	    		      nm++;
	    		      gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
	    		      break;
	    		    }
	    		  }
	    		} 
			    else {
	    	    gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
	    	    nm++;
	    	    }
	          }
	          if (mm->pdgId()==-13) { foundit++;
	    		if (mm->status()!=1) {
	    	  		for (size_t m=0; m<mm->numberOfDaughters(); m++) {
	    	  		  const reco::Candidate *mu = mm->daughter(m);
	    	  		  if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
	    	  		    nm++;
	    	  		    gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
	    	  		    break;
	    	  		  }
	    	  		}
	    		}
			  	else {
	    	  		gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
	    	  		nm++;
	    	  	}
	          }
	        }// end for daugters of Jpsi
	        if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	        else foundit-=nm;
	      }//end if Jpsi
	    } // end for B daughters for jpsi
        for (size_t k=0; k<dau->numberOfDaughters(); k++){
			const reco::Candidate *gdau = dau->daughter(k);
			if (gdau->pdgId()==3122){ //is Lambda0
			  foundit++;
			  gen_ks0_vtx.SetXYZ(gdau->vx(), gdau->vy(), gdau->vz());
			  gen_ks0_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
			  gen_ks0_ct = GetLifetime(gen_ks0_p4, gen_jpsi_vtx, gen_ks0_vtx); 
			  // TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx 
		    }// end if K0s
		}// end of B daughters for Ks0
      } // end in B0
      if (foundit>=5) break; //1-B0, 2-JPsi, 3-mu1, 4-mu2, 5-Ks0,// NO PIONS 6-pi1, 7-pi2
    } // for i
    //if (foundit!=5) {
    //  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_ks0_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_b_vtx.SetXYZ(0.,0.,0.);
  	//  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  	//  gen_ks0_vtx.SetXYZ(0.,0.,0.);
  	//  gen_b_ct = -9999.;
  	//  gen_ks0_ct = -9999.;
    //  std::cout << "Does not found the given decay (res) " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    //}
  }
  
  
  // for the non-resonant channel 
  if ( (isMC_ || OnlyGen_) && pruned.isValid() && !isRes_) {
	std::cout<< "only gen test ok "<< std::endl;  
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 511) ) { //&& (dau->status() == 2) ) { //B0 
	    foundit++;
		//std::cout<< "|--XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX----XXX--|"<< std::endl;
		//std::cout<<"Found B0 printing Decay Tree ..."<< std::endl;
		//printMCtree(dau, 0);
	    gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	    gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
		int nm=0;
	    for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	      const reco::Candidate *mm = dau->daughter(k);     
	      if (mm->pdgId()==13 && !isAncestor(443,mm)) { foundit++;  // cames from B but not J/p
		    gen_jpsi_vtx.SetXYZ(mm->vx(),mm->vy(),mm->vz());
	        if (mm->status()!=1) {
	          for (size_t m=0; m<mm->numberOfDaughters(); m++) {
	            const reco::Candidate *mu = mm->daughter(m);
	            if (mu->pdgId()==13 ) { //&& mu->status()==1) {
	              nm++;
	              gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
	              break;
	            }
	          }
	        } 
		    else {
	          gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
	          nm++;
	        }
	      }
	      if (mm->pdgId()==-13 && !isAncestor(443,mm)) { foundit++;  // cames from B but not J/p
	        if (mm->status()!=1) {
	        	for (size_t m=0; m<mm->numberOfDaughters(); m++) {
	        	  const reco::Candidate *mu = mm->daughter(m);
	        	  if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
	        	    nm++;
	        	    gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
	        	    break;
	        	  }
	        	}
	        }
		    else {
	        	gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
	        	nm++;
	        }
	      }
	    }// end for daugters (for muons)
	    if (nm==2) {
			gen_jpsi_p4 = gen_muon1_p4 + gen_muon2_p4;
	        gen_b_ct = GetLifetime(gen_b_p4,gen_b_vtx,gen_jpsi_vtx);
			//std::cout << "muons ok"<< std::endl;
		} 	  
	    else foundit-=nm;
	    // end for B daughters for dimuon
        for (size_t k=0; k<dau->numberOfDaughters(); k++){
			const reco::Candidate *gdau = dau->daughter(k); 
			if (gdau->pdgId()==310){ //is K0s
			  foundit++;
			  gen_ks0_vtx.SetXYZ(gdau->vx(), gdau->vy(), gdau->vz());
			  gen_ks0_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
			  gen_ks0_ct = GetLifetime(gen_ks0_p4, gen_jpsi_vtx, gen_ks0_vtx); 
			   
		    }// end if K0s
	    }// end for B daughters for Ks0
      } // end if B0
      if (foundit>=5) break; //1-B0, 2-JPsi, 3-mu1, 4-mu2, 5-Ks0 //NOPIONS , 6-pi1, 7-pi2
    } // for gen particlea
    //if (foundit!=5) {
    //  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_ks0_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  	//  gen_b_vtx.SetXYZ(0.,0.,0.);
  	//  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  	//  gen_ks0_vtx.SetXYZ(0.,0.,0.);
  	//  gen_b_ct = -9999.;
  	//  gen_ks0_ct = -9999.;
    //  std::cout << "Does not found the given decay (non-res) " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    //}
  }
  //std::cout<< "Decay info: B0: "  << gen_b_p4.M() << std::endl;
  //std::cout<< "Decay info: Mu1: " << gen_muon1_p4.M() << std::endl;
  //std::cout<< "Decay info: Mu2: " << gen_muon2_p4.M() << std::endl;
  //std::cout<< "Decay info: Ks0: " << gen_ks0_p4.M()   << std::endl;
  //std::cout<< "Decay info: Jspi: "<< gen_jpsi_p4.M()  << std::endl;
  
  //Trigger info 
  //origen /afs/cern.ch/user/j/jmejiagu/work/public/data2015_RunII/HI/CMSSW_8_0_31/src/Ponia/OniaPhoton/src/Bu_JpsiK_PAT.cc
  trigger = 0;
  if ( triggerResults_handle.isValid()) {
   //std::cout << "Triggers ok ..." <<std::endl;	  
   const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
   //std::vector<std::string> const& names = TheTriggerNames.triggerNames(); //new 
   //for (unsigned i = 0; i < TheTriggerNames.size(); ++i) {
   //   std::cout << names[i] << "  " << TheTriggerNames.triggerName(i) << std::endl;
   //}
   unsigned int NTRIGGERS = 20;
   // de acuerdo con https://indico.cern.ch/event/988495/contributions/4161361/attachments/2166530/3656840/Slides_210104_ERDUpdate.pdf
   // los Trigger de BParking son : HLT_Mu7_IP4, HLT_Mu9_IP5, HLT_Mu9_IP6, HLT_Mu12_IP6
   //#V2.0
   //"HLT_Mu9_IP6_part","HLT_Mu8p5_IP3p5","HLT_Mu10p5_IP3p5","HLT_Mu8_IP3",
   //#V2.2
   //"HLT_Mu12_IP6","HLT_Mu9_IP5","HLT_Mu7_IP4","HLT_Mu9_IP4","HLT_Mu8_IP5","HLT_Mu8_IP6",
   //#V3.5
   //"HLT_Mu9_IP3","HLT_Mu9_IP0")
   std::string TriggersToTest[NTRIGGERS] = {
     "HLT_Mu12_IP6", //0
	 "HLT_Mu9_IP0","HLT_Mu9_IP3", "HLT_Mu9_IP4", "HLT_Mu9_IP5", "HLT_Mu9_IP6", //1-5
     "HLT_Mu8_IP3","HLT_Mu8_IP5", "HLT_Mu8_IP6", //6-8
     "HLT_Mu7_IP4", //9
  	 "L1_SingleMu22", "L1_SingleMu25", "L1_SingleMu18er1p5", "L1_SingleMu14er1p5", "L1_SingleMu12er1p5", "L1_SingleMu10er1p5", //10-15
	 "L1_SingleMu9er1p5", "L1_SingleMu8er1p5", "L1_SingleMu7er1p5", "L1_SingleMu6er1p5"}; //16-19

   for (unsigned int i = 0; i < NTRIGGERS; i++) {
	 for (unsigned int h = 0; h < TheTriggerNames.size(); ++h){
		std::string triggerName = TheTriggerNames.triggerName(h); 
		if (triggerName.find(TriggersToTest[i]) != std::string::npos ){
			trigger += (1<<i);
			//std::cout<< TriggersToTest[i] <<" found in " << triggerName << " STR !!" << std::endl; 
		} 
		//std::cout<< "Trigger to check: " << TheTriggerNames.triggerName(h)<< " | looking for: "<< TriggersToTest[i] << std::endl;
	 }
     //for (int version = 1; version < 9; version++) {
     //  std::stringstream ss;
     //  ss << TriggersToTest[i] << "_v" << version;
	 //  //std::cout << "Looking for: " << ss.rdbuf() << std::endl; 
     //  unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
	 //  //std::cout << "Resulting bit: "<< bit << std::endl;
	 //  //std::cout << "trigger handle size: "<< triggerResults_handle->size() << std::endl;
     //  if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
     //    trigger += (1<<i);
     //    break;
     //  }
     //}
   }
  } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
  
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  reco::Vertex bestVtxBS;

  // get primary vertex Change to get PV by best Pointing Angle
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);
  //bestVtx = *(recVtxs->begin());
  //
  //priVtxX = bestVtx.x();
  //priVtxY = bestVtx.y();
  //priVtxZ = bestVtx.z();
  //priVtxXE = bestVtx.covariance(0, 0);
  //priVtxYE = bestVtx.covariance(1, 1);
  //priVtxZE = bestVtx.covariance(2, 2);
  //priVtxXYE = bestVtx.covariance(0, 1);
  //priVtxXZE = bestVtx.covariance(0, 2);
  //priVtxYZE = bestVtx.covariance(1, 2);
  //priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 

  nVtx = recVtxs->size();

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event(); 


  //emulate BParking MuonTriggerSelector 
  std::vector<pat::TriggerObjectStandAlone> triggeringMuons;
  int int_obj = 0;
  bool debug = false;
  const double maxdR_ = 0.1;
  // std::cout << "\n\n\n------------>>>>>>NEW RECORD NEW RECORD NEW RECORD NEW RECORD"<<"\n";
 
  //Itera sobre los objetos 
  //corrobora que al menos haya un filterID = 83
  //corrobora que al menos uno tenga 'hltL3' y 'Park'
  // el objeto que satisfaga estas dos condiciones se agrega
  // al vector triggeringMuons, definido previamente
  // std::cout << "\n\n ------------>>>>>>>>TriggerObjectStandAlone TriggerObjectStandAlone TriggerObjectStandAlone TriggerObjectStandAlone"<<"\n";
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames   
    int_obj++; 
    // std::cout << "---->>obj number = " <<int_obj << "\n";
    obj.unpackFilterLabels(iEvent, *triggerResults_handle);
    obj.unpackPathNames(names);
    bool isTriggerMuon = false;

    // checa que al menos un elemento sea un muon (ID=83)
    // que pasa con los demas?
    // std::cout << "\tfilterIds size:   " << obj.filterIds().size()<< "\n";
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
    // std::cout << "\tfilterLabels size:  " << obj.filterLabels().size()<< "\n";
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){   
        std::string filterName = obj.filterLabels()[h];
        //std::cout << "\t\tfilterlabes:  " << h << filterName << "\n";
        if(filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
            isTriggerMuon = true;
            //std::cout << "\t\tVVVVVVVV  Filter:   " << filterName<<"\n"; 
        }
	    //else{ isTriggerMuon = false; }
    }
    //std::cout << "\n\n\n";
    if(!isTriggerMuon) continue;
    triggeringMuons.push_back(obj);


    if(debug){ 
        std::cout << "\t\t\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      	//Print trigger object collection and type
	    std::cout << "\t\t\tCollection: " << obj.collection() << std::endl;
    }
  }//trigger objects

  if(debug){    
    // std::cout << "\n Total n of triggering muons = " << triggeringMuons.size() << std::endl;
    for(auto ij : triggeringMuons){
	    std::cout << " \t>>> components (pt, eta, phi) = (" << ij.pt() << ", " << ij.eta() << ", " << ij.phi() << ")\n";
    }
  }
  
  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  unsigned int nMu_tmp = thePATMuonHandle->size();
  if (!OnlyGen_){
   for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  //std::cout << "muons working "<< std::endl;		
	  if(iMuon1==iMuon2) continue;

	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------
	  //std::cout<< "XXXXXXXXXXXXXXXXXXXXX MUON COMPARISON XXXXXXXXXXXXXXXXXXXXX"<< std::endl;
	  //if(iMuon1->charge() == -1){
	  //	std::cout<< " MUON - REC( pT: "<< iMuon1->pt() << ", Eta: "<<iMuon1->eta() << ", Phi: "<< iMuon1->phi() << std::endl;
	  //	std::cout<< " MUON - GEN( pT: "<< gen_muon1_p4.Pt() << ", Eta: "<<gen_muon1_p4.Eta() << ", Phi: "<< gen_muon1_p4.Phi() << std::endl;
	  //	std::cout<< "\n" <<std::endl;
	  //	std::cout<< " MUON + REC( pT: "<< iMuon2->pt() << ", Eta: "<<iMuon2->eta() << ", Phi: "<< iMuon2->phi() << std::endl;
	  //	std::cout<< " MUON + GEN( pT: "<< gen_muon2_p4.Pt() << ", Eta: "<<gen_muon2_p4.Eta() << ", Phi: "<< gen_muon2_p4.Phi() << std::endl;  
	  //}
	  //else{
	  //	std::cout<< " MUON - REC( pT: "<< iMuon2->pt() << ", Eta: "<<iMuon2->eta() << ", Phi: "<< iMuon2->phi() << std::endl;
	  //	std::cout<< " MUON - GEN( pT: "<< gen_muon1_p4.Pt() << ", Eta: "<<gen_muon1_p4.Eta() << ", Phi: "<< gen_muon1_p4.Phi() << std::endl;
	  //	std::cout<<"\n"<<std::endl;  
	  //	std::cout<< " MUON + REC( pT: "<< iMuon1->pt() << ", Eta: "<<iMuon1->eta() << ", Phi: "<< iMuon1->phi() << std::endl;
	  //	std::cout<< " MUON + GEN( pT: "<< gen_muon2_p4.Pt() << ", Eta: "<<gen_muon2_p4.Eta() << ", Phi: "<< gen_muon2_p4.Phi() << std::endl;
	  //} 
	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<1.0) continue;
	  if(iMuon2->track()->pt()<1.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
	  
	  //Let's check the vertex and mass
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	  // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  //Match with TriggerMuons, BParking Nano MuonTriggerSelector emulation 
	  float dRMuonMatching1 = -1.; 	
	  float dRMuonMatching2 = -1.;
	  float dzm1_trg = 10.0;
	  float dzm2_trg = 10.0;
	   	
      for(unsigned int iTrg=0; iTrg<triggeringMuons.size(); ++iTrg){
        float dR1 = reco::deltaR(triggeringMuons[iTrg], *iMuon1);
        float dR2 = reco::deltaR(triggeringMuons[iTrg], *iMuon2);
        // std::cout << "\n\t\tdR = " << dR << "\n";
		if (abs(triggeringMuons[iTrg].vz() - iMuon1->vz()) < dzm1_trg){
		  dzm1_trg = abs(triggeringMuons[iTrg].vz() - iMuon1->vz());
		}
		if (abs(triggeringMuons[iTrg].vz() - iMuon2->vz())) {
		  dzm2_trg = abs(triggeringMuons[iTrg].vz() - iMuon2->vz());
		}
	    if((dR1 < dRMuonMatching1 || dRMuonMatching1 == -1) && dR1 < maxdR_){
        	dRMuonMatching1 = dR1;
        	// float eta = muon.eta() - triggeringMuons[iTrg].eta();
        	// float phi = muon.phi() - triggeringMuons[iTrg].phi();
        	// dR_H = std::sqrt(eta*eta+phi*phi);
        	// std::cout << "\n\t\t dR_H"<< iTrg <<" = " << dR_H
        	//   << "\n\t\treco (pt, eta, phi) = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << " " 
        	//   << "\n\t\tHLT (pt, eta, phi)  = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " " << triggeringMuons[iTrg].phi()
        	//   << std::endl;
	    }
	    if((dR2 < dRMuonMatching2 || dRMuonMatching2 == -1) && dR2 < maxdR_){
        	dRMuonMatching2 = dR2;	 
		}
      }
      //if ((dRMuonMatching1 != -1) | (dRMuonMatching2 != -1)) std::cout << "matching ok ..." << std::endl;
	  //std::cout << "distaNCE TO TRIGGER :" << dzm1_trg << " & " << dzm2_trg << std::endl;
	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************

	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    //std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }

	  KinematicParticleVertexFitter fitter;   

	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    //std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }
      //std::cout << "pass fit continues ... "<< std::endl;
	  if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }
      try{
	  	psiVertexFitTree->movePointerToTheTop();
	  }
	  catch(...){
		if(debug) std::cout<< "psi vertex empry or not valid" <<std::endl;
		continue;  
	  }
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }
	  
	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	   //if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;
	   if(psi_vFit_noMC->currentState().mass()<1.0 || psi_vFit_noMC->currentState().mass()>22.0) continue;

	   //  ***************  

	   if ( theV0PtrHandle->size()>0 && thePATMuonHandle->size()>=2 ) 
	     {
	       
	       for ( vector<VertexCompositePtrCandidate>::const_iterator iVee = theV0PtrHandle->begin();   iVee != theV0PtrHandle->end(); ++iVee )
		 {
			 
		     //get Lam tracks from V0 candidate fot the B --> Ks0 Decay, omit all continues 
		     vector<pat::PackedCandidate> v0daughters;
		     vector<Track> theDaughterTracks;
			 
			 const pat::PackedCandidate* track1 = dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(0));
			 const pat::PackedCandidate* track2 = dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(1));

			 if (track1->pt() < 0.55 || track2->pt() < 0.55) continue;


		     v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(0))) );
		     v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(1))) );
		     		     
		     for(unsigned int j = 0; j < v0daughters.size(); ++j)
		       {
			 theDaughterTracks.push_back(v0daughters[j].pseudoTrack());
		       }
		     			
		     // it does not have sences here. 
		     //if ( IsTheSame(*theDaughterTracks[0],*iMuon1) || IsTheSame(*theDaughterTracks[0],*iMuon2) ) continue;
		     //if ( IsTheSame(*theDaughterTracks[1],*iMuon1) || IsTheSame(*theDaughterTracks[1],*iMuon2) ) continue;
		     
		     //Now let's see if these two tracks make a vertex
		     reco::TransientTrack pion1TT((*theB).build(theDaughterTracks[0]));
		     reco::TransientTrack pion2TT((*theB).build(theDaughterTracks[1]));		     
		     
		     ParticleMass pion_mass = 0.13957018;
		     ParticleMass Ks0_mass = 0.497614;
		     float pion_sigma = pion_mass*1.e-6;
		     float Ks0_sigma = Ks0_mass*1.e-6;
		     ///*
		     //initial chi2 and ndf before kinematic fits.
		     float chi = 0.;
		     float ndf = 0.;
		     vector<RefCountedKinematicParticle> pionParticles;
		     // vector<RefCountedKinematicParticle> muonParticles;
		     try {
		       pionParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
		       pionParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
		     }
		     catch(...) {
		       //std::cout<<" Exception caught ... continuing 3 "<<std::endl;
		       continue;
		     }
		     
		     RefCountedKinematicTree Ks0VertexFitTree;
		     try{
		       Ks0VertexFitTree = fitter.fit(pionParticles); 
		     }
		     catch(...) {
		       //std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
		       continue;
		     }
		     if (!Ks0VertexFitTree->isValid()) {
			 //std::cout << "invalid vertex from the Ks0 vertex fit" << std::endl;
			 //continue; 
		     }
			 bool goodKs0VtxFit = true;  
			 float b_mass_tmp = -1.;
			 float Ks0_mass_tmp = -1.;
			 try{  
		     	Ks0VertexFitTree->movePointerToTheTop();
			 }
			 catch(...){
				if (debug) std::cout<<"Ks0 Vertex fit empty or not valid" << std::endl;
				goodKs0VtxFit = false;
				//continue; we are not interested in the Ks0 Fit
			 }
			 if (goodKs0VtxFit){
		     	RefCountedKinematicParticle Ks0_vFit_noMC = Ks0VertexFitTree->currentParticle();
		     	RefCountedKinematicVertex Ks0_vFit_vertex_noMC = Ks0VertexFitTree->currentDecayVertex();
		     	//aca chinga a su madre todo

             	try {
		     	   if( Ks0_vFit_vertex_noMC->chiSquared() < 0 )
		     	     { 
			 	   //std::cout << "negative chisq from ks fit" << endl;
			 	   //continue;
		     	     }
		     	   //std::cout << "pass Fit continues ... "<< std::endl;
		     	   //some loose cuts go here

		     	   //if(Ks0_vFit_vertex_noMC->chiSquared()>50) continue;
		     	   //if(Ks0_vFit_noMC->currentState().mass()<0.45 || Ks0_vFit_noMC->currentState().mass()>0.55) continue; 

		     	   Ks0VertexFitTree->movePointerToTheFirstChild();
		     	   RefCountedKinematicParticle T1CandMC = Ks0VertexFitTree->currentParticle();

		     	   Ks0VertexFitTree->movePointerToTheNextChild();
		     	   RefCountedKinematicParticle T2CandMC = Ks0VertexFitTree->currentParticle();

		     	   //  Ks0  mass constrain
		     	   // do mass constrained vertex fit
		     	   // creating the constraint with a small sigma to put in the resulting covariance 
		     	   // matrix in order to avoid singularities
		     	   // JPsi mass constraint is applied in the final B fit

		     	   KinematicParticleFitter csFitterKs;
		     	   KinematicConstraint * ks_c = new MassKinematicConstraint(Ks0_mass,Ks0_sigma);
		     	   // add mass constraint to the ks0 fit to do a constrained fit:  

		     	   Ks0VertexFitTree = csFitterKs.fit(ks_c,Ks0VertexFitTree);
		     	   if (!Ks0VertexFitTree->isValid()){
		     	     //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
		     	     //continue; 
		     	   }
		     	   //std::cout << "pass 424 continues ... "<< std::endl;
			 	   try {
		     	   	Ks0VertexFitTree->movePointerToTheTop();
					}
					catch(...){
						if(debug) std::cout << "Ks0 mass constrained vertex fit empty or not valid" << std::endl;
						//continue; we are not interested in the ks0 fit 
					}
		     	   RefCountedKinematicParticle ks0_vFit_withMC = Ks0VertexFitTree->currentParticle();
		     	   //Now we are ready to combine!
		     	   // JPsi mass constraint is applied in the final Bd fit,
		     	   vector<RefCountedKinematicParticle> vFitMCParticles;
		     	   vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		     	   vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		     	   vFitMCParticles.push_back(ks0_vFit_withMC);
		     	   MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
		     	   //KinematicConstrainedVertexFitter kcvFitter;
		     	   //RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
			 	   //no mass constrain 
			 	   KinematicParticleVertexFitter kcvFitter;
			 	   RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles);
		     	   if (!vertexFitTree->isValid()) {
		     	     //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		     	     //continue;
		     	   }
					try {
		     	   	vertexFitTree->movePointerToTheTop();		     
					}
					catch (...){
						if(debug) std::cout<< "B Vertex empty or not valid" << std::endl;
						//continue ; // we are not interested in the B decay ...
					}
		     	   RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		     	   RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		     	   if (!bDecayVertexMC->vertexIsValid()){
		     	     //std::cout << "B MC fit vertex is not valid" << endl;
		     	     //continue;
		     	   }
		     	   //if(bCandMC->currentState().mass()<4.5 || bCandMC->currentState().mass()>6.0) continue;
		     	   if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ) 
		     	     {
			 	       //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
			 	       //continue;
		     	     }
		     	   //std::cout << "pass 461 continues ... "<< std::endl;
		     	   double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		     	   if(B_Prob_tmp<0.01)
		     	     {
			 	      //continue;
		     	     }		     
		     	   //std::cout << "pass 467" <<std::endl;
					b_mass_tmp = bCandMC->currentState().mass();
					Ks0_mass_tmp =  Ks0_vFit_noMC->currentState().mass();

			 }
			 catch (...){
				 std::cout << "B0 Ks0 Vertex fit chinga a su madre " << std::endl;
			 }
			 }
             //*/
			 //MAKE LAMBDA FIT 
             //Now let's see if these two tracks make a vertex
			 //if (debug) std::cout << "Ks0 - B0 Vtx ok" << std::endl;
			 int pion_I = 0;
			 int proton_I = 0;
			 if (v0daughters[0].p() > v0daughters[1].p() ){
				pion_I = 1;
				proton_I = 0; 
			 }
			 else {
				pion_I = 0;
				proton_I = 1;
			 }
		     reco::TransientTrack pionTT((*theB).build(theDaughterTracks[pion_I]));
		     reco::TransientTrack protonTT((*theB).build(theDaughterTracks[proton_I]));		     
		     
		     ParticleMass proton_mass = 0.938272;
			 ParticleMass lambda0_mass = 1.115683;

			 float proton_sigma = proton_mass*1.e-6;
		     float lambda0_sigma = lambda0_mass*1.e-6;
		     
		     //initial chi2 and ndf before kinematic fits.
		     //float chi = 0.;
		     //float ndf = 0.;
		     vector<RefCountedKinematicParticle> trackParticles;
		     // vector<RefCountedKinematicParticle> muonParticles;
		     try {
		       trackParticles.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));
		       trackParticles.push_back(pFactory.particle(protonTT,proton_mass,chi,ndf,proton_sigma));
		     }
		     catch(...) {
		       //std::cout<<" Exception caught ... continuing 3 "<<std::endl;
		       continue;
		     }
		     
		     RefCountedKinematicTree Lam0VertexFitTree;
		     try{
		       Lam0VertexFitTree = fitter.fit(trackParticles); 
		     }
		     catch(...) {
		       //std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
		       continue;
		     }
		     if (!Lam0VertexFitTree->isValid()) 
		       {
			 //std::cout << "invalid vertex from the Ks0 vertex fit" << std::endl;
			   continue; 
		       }
			 try {  
		     	Lam0VertexFitTree->movePointerToTheTop();
			 }
			 catch (...){
				if (debug) std::cout<< "Lamb0 Vertex Fit Tree empty or not valid"<< std::endl;
				continue;
			 }
		     RefCountedKinematicParticle Lam0_vFit_noMC = Lam0VertexFitTree->currentParticle();
		     RefCountedKinematicVertex Lam0_vFit_vertex_noMC = Lam0VertexFitTree->currentDecayVertex();
		     
		     if( Lam0_vFit_vertex_noMC->chiSquared() < 0 )
		       { 
			 //std::cout << "negative chisq from ks fit" << endl;
			 //continue;
		       }
		     //std::cout << "pass Fit continues ... "<< std::endl;
		     //some loose cuts go here
		     
		     //if(Ks0_vFit_vertex_noMC->chiSquared()>50) continue;
		     //if(Ks0_vFit_noMC->currentState().mass()<0.45 || Ks0_vFit_noMC->currentState().mass()>0.55) continue; 
		     
		     Lam0VertexFitTree->movePointerToTheFirstChild();
		     RefCountedKinematicParticle TpiCandMC = Lam0VertexFitTree->currentParticle();
		     
		     Lam0VertexFitTree->movePointerToTheNextChild();
		     RefCountedKinematicParticle TprCandMC = Lam0VertexFitTree->currentParticle();
		     
		     //  Lambda0  mass constrain
		     // do mass constrained vertex fit
		     // creating the constraint with a small sigma to put in the resulting covariance 
		     // matrix in order to avoid singularities
		     // JPsi mass constraint is applied in the final LambdaB fit
		     
		     KinematicParticleFitter csFitterL0;
		     KinematicConstraint * l0_c = new MassKinematicConstraint(lambda0_mass,lambda0_sigma);
		     // add mass constraint to the ks0 fit to do a constrained fit:  
		     
		     Lam0VertexFitTree = csFitterL0.fit(l0_c,Lam0VertexFitTree);
		     if (!Lam0VertexFitTree->isValid()){
		       //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
		       continue; 
		     }
		     //std::cout << "pass 424 continues ... "<< std::endl;
		     //aca chinga a su madre todo
			 try{
		     	Lam0VertexFitTree->movePointerToTheTop();
			 }
			 catch(...){
				if (debug) std::cout<< "Lamb0 mass constrained Vertex Fit Tree empty or not valid" <<std::endl; 
				continue;
			 }

			 RefCountedKinematicParticle Lam0_vFit_withMC = Lam0VertexFitTree->currentParticle();
		     
		     //Now we are ready to combine!
		     
		     vector<RefCountedKinematicParticle> LvFitMCParticles;
		     LvFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		     LvFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		     LvFitMCParticles.push_back(Lam0_vFit_withMC);
		     
		     //KinematicConstrainedVertexFitter kcvFitter;
		     //RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);

			 //no mass constrain 
			 KinematicParticleVertexFitter LbcvFitter;
			 RefCountedKinematicTree LbvertexFitTree = LbcvFitter.fit(LvFitMCParticles);
		     if (!LbvertexFitTree->isValid()) {
		       //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		       //continue;
		     }
		     
			 try {
		     	LbvertexFitTree->movePointerToTheTop();		     
			 }
			 catch (...){
				 if (debug) std::cout<<"Lambda B Vertex empty or not valid" << std::endl;
				 continue;
			 }
		     RefCountedKinematicParticle LbCandMC = LbvertexFitTree->currentParticle();
		     RefCountedKinematicVertex LbDecayVertexMC = LbvertexFitTree->currentDecayVertex();
		     if (!LbDecayVertexMC->vertexIsValid()){
		       //std::cout << "B MC fit vertex is not valid" << endl;
		       continue;
		     }
		     
		     if(LbCandMC->currentState().mass()<4.7 || LbCandMC->currentState().mass()>6.7) continue;
		     
		     if(LbDecayVertexMC->chiSquared()<0 || LbDecayVertexMC->chiSquared()>50 ) 
		       {
			     //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
			     continue;
		       }
		     //std::cout << "pass 461 continues ... "<< std::endl;
		     double LB_Prob_tmp       = TMath::Prob(LbDecayVertexMC->chiSquared(),(int)LbDecayVertexMC->degreesOfFreedom());
		     if(LB_Prob_tmp<0.01)
		       {
			    continue;
		       }		     



			 //END LAMBDA FIT 
		   // get children from final B fit
		   LbvertexFitTree->movePointerToTheFirstChild();
		   RefCountedKinematicParticle mu1CandMC = LbvertexFitTree->currentParticle();
		   LbvertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle mu2CandMC = LbvertexFitTree->currentParticle();
		   
		   LbvertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle L0CandMC = LbvertexFitTree->currentParticle();
		   
		   KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMupKP;
		   KinematicParameters psiMumKP;
	       
		   if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
		   if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
		   if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
		   if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

 		   GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
				       mu1CandMC->currentState().globalMomentum().y(),
 				       mu1CandMC->currentState().globalMomentum().z());

 		   GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
				       mu2CandMC->currentState().globalMomentum().y(),
 				       mu2CandMC->currentState().globalMomentum().z());

 		   GlobalVector L0p1vec(TprCandMC->currentState().globalMomentum().x(),
				        TprCandMC->currentState().globalMomentum().y(),
 				        TprCandMC->currentState().globalMomentum().z());

 		   GlobalVector L0p2vec(TpiCandMC->currentState().globalMomentum().x(),
					TpiCandMC->currentState().globalMomentum().y(),
					TpiCandMC->currentState().globalMomentum().z());

		   KinematicParameters L0PioP = TpiCandMC->currentState().kinematicParameters();
		   KinematicParameters L0ProP = TprCandMC->currentState().kinematicParameters();
		   KinematicParameters Ks0PipKP;
		   KinematicParameters Ks0PimKP;
	       
		   if ( TpiCandMC->currentState().particleCharge() > 0 ) Ks0PipKP = L0PioP;
		   if ( TpiCandMC->currentState().particleCharge() < 0 ) Ks0PimKP = L0PioP;
		   if ( TprCandMC->currentState().particleCharge() > 0 ) Ks0PipKP = L0ProP;
		   if ( TprCandMC->currentState().particleCharge() < 0 ) Ks0PimKP = L0ProP;	 
		   //Get bestPVtx by best ointing angle 	
		   // ********************* loop over all the primary vertices and we choose the one with the best pointing angle ****************

           Double_t priVtxX_t = -10000.0;
           Double_t priVtxY_t = -10000.0;
           Double_t priVtxZ_t = -10000.0;
           Double_t priVtxXE_t = -10000.0;
           Double_t priVtxYE_t = -10000.0;
           Double_t priVtxZE_t = -10000.0;
           Double_t priVtxXZE_t = -10000.0;
           Double_t priVtxXYE_t = -10000.0;
           Double_t priVtxYZE_t = -10000.0;
           Double_t priVtxCL_t = -10000.0;
           Double_t lip1 = -1000000.0;
		   Double_t PVTriggDz_t = 100000.0;
		   unsigned int TrkIndex_t = 0;
		   int nTks_t = 0;
           for(size_t i = 0; i < recVtxs->size(); ++i) {
                const reco::Vertex &vtx = (*recVtxs)[i];

                Double_t dx1 = (*LbDecayVertexMC).position().x() - vtx.x();
                Double_t dy1 = (*LbDecayVertexMC).position().y() - vtx.y();
                Double_t dz1 = (*LbDecayVertexMC).position().z() - vtx.z();
                float cosAlphaXYb1 = ( LbCandMC->currentState().globalMomentum().x() * dx1 + LbCandMC->currentState().globalMomentum().y()*dy1 + LbCandMC->currentState().globalMomentum().z()*dz1  )/( sqrt(dx1*dx1+dy1*dy1+dz1*dz1)* LbCandMC->currentState().globalMomentum().mag() );

                if(cosAlphaXYb1>lip1){
                    lip1 = cosAlphaXYb1 ;
                    priVtxX_t = vtx.x();
                    priVtxY_t = vtx.y();
                    priVtxZ_t = vtx.z();
                    priVtxXE_t = vtx.covariance(0, 0);
                    priVtxYE_t = vtx.covariance(1, 1);
                    priVtxZE_t = vtx.covariance(2, 2);
                    priVtxXYE_t = vtx.covariance(0, 1);
                    priVtxXZE_t = vtx.covariance(0, 2);
                    priVtxYZE_t = vtx.covariance(1, 2);
                    priVtxCL_t = ChiSquaredProbability((double)(vtx.chi2()),(double)(vtx.ndof())); 
					TrkIndex_t = i;
					nTks_t    = vtx.tracksSize();
					for(unsigned int iTrg=0; iTrg<triggeringMuons.size(); ++iTrg){
						Double_t PVTriggDz_tt = abs(triggeringMuons[iTrg].vz() - vtx.z());
						if (PVTriggDz_tt < PVTriggDz_t){
							PVTriggDz_t = PVTriggDz_tt;
						}
					}
                    bestVtx = vtx;
                }
            }

		   // fill candidate variables now
		   
		   if(nB==0){		    
		     nMu  = nMu_tmp;
		     // cout<< "*Number of Muons : " << nMu_tmp << endl;
		   } // end nB==0		     
		
		   B_mass->push_back(b_mass_tmp);
		   LB_mass->push_back(LbCandMC->currentState().mass());
		   B_px->push_back(LbCandMC->currentState().globalMomentum().x());
		   B_py->push_back(LbCandMC->currentState().globalMomentum().y());
		   B_pz->push_back(LbCandMC->currentState().globalMomentum().z());
		
		   B_Ks0_mass->push_back( Ks0_mass_tmp );
		   LB_L0_mass->push_back( Lam0_vFit_noMC->currentState().mass() );
		   B_Ks0_px->push_back( Lam0_vFit_noMC->currentState().globalMomentum().x() );
		   B_Ks0_py->push_back( Lam0_vFit_noMC->currentState().globalMomentum().y() );
		   B_Ks0_pz->push_back( Lam0_vFit_noMC->currentState().globalMomentum().z() );

		   B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		   B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		   B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		   B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );

		   B_Ks0_pt1->push_back(L0p1vec.perp());
		   B_Ks0_px1->push_back(L0ProP.momentum().x());
		   B_Ks0_py1->push_back(L0ProP.momentum().y());
		   B_Ks0_pz1->push_back(L0ProP.momentum().z());
		   B_Ks0_px1_track->push_back(v0daughters[proton_I].px());
		   B_Ks0_py1_track->push_back(v0daughters[proton_I].py());
		   B_Ks0_pz1_track->push_back(v0daughters[proton_I].pz());
		   B_Ks0_charge1->push_back(TprCandMC->currentState().particleCharge());

		   B_Ks0_pt2->push_back(L0p2vec.perp());
		   B_Ks0_px2->push_back(L0PioP.momentum().x());
		   B_Ks0_py2->push_back(L0PioP.momentum().y());
		   B_Ks0_pz2->push_back(L0PioP.momentum().z());
		   B_Ks0_px2_track->push_back(v0daughters[pion_I].px());
		   B_Ks0_py2_track->push_back(v0daughters[pion_I].py());
		   B_Ks0_pz2_track->push_back(v0daughters[pion_I].pz());
		   B_Ks0_charge2->push_back(TpiCandMC->currentState().particleCharge());

		   B_J_pt1->push_back(Jp1vec.perp());
		   B_J_px1->push_back(psiMu1KP.momentum().x());
		   B_J_py1->push_back(psiMu1KP.momentum().y());
		   B_J_pz1->push_back(psiMu1KP.momentum().z());
		   B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());

		   B_J_pt2->push_back(Jp2vec.perp());
		   B_J_px2->push_back(psiMu2KP.momentum().x());
		   B_J_py2->push_back(psiMu2KP.momentum().y());
		   B_J_pz2->push_back(psiMu2KP.momentum().z());
		   B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());

		   B_Ks0_chi2->push_back(Lam0_vFit_vertex_noMC->chiSquared());
		   B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		   B_chi2->push_back(LbDecayVertexMC->chiSquared());

		   //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		   double L0_Prob_tmp  = TMath::Prob(Lam0_vFit_vertex_noMC->chiSquared(),(int)Lam0_vFit_vertex_noMC->degreesOfFreedom());
		   B_Prob    ->push_back(LB_Prob_tmp);
		   B_J_Prob  ->push_back(J_Prob_tmp);
		   B_ks0_Prob ->push_back(L0_Prob_tmp);

           //Best Pointing angle PV
           priVtxX->push_back(priVtxX_t);
		   priVtxY->push_back(priVtxY_t);
		   priVtxZ->push_back(priVtxZ_t);
		   priVtxXE->push_back(priVtxXE_t); 
		   priVtxYE->push_back(priVtxYE_t); 
		   priVtxZE->push_back(priVtxZE_t); 
		   priVtxXZE->push_back(priVtxXZE_t);
		   priVtxXYE->push_back(priVtxXYE_t);
		   priVtxYZE->push_back(priVtxYZE_t);
		   priVtxCL->push_back(priVtxCL_t); 
		   TrkIndex->push_back(TrkIndex_t);
           nTks->push_back(nTks_t);
		   PVTriggDz->push_back(PVTriggDz_t);

	       // ************
		   bDecayVtxX->push_back((*LbDecayVertexMC).position().x());
		   bDecayVtxY->push_back((*LbDecayVertexMC).position().y());
		   bDecayVtxZ->push_back((*LbDecayVertexMC).position().z());
		   bDecayVtxXE->push_back(LbDecayVertexMC->error().cxx());
		   bDecayVtxYE->push_back(LbDecayVertexMC->error().cyy());
		   bDecayVtxZE->push_back(LbDecayVertexMC->error().czz());
		   bDecayVtxXYE->push_back(LbDecayVertexMC->error().cyx());
		   bDecayVtxXZE->push_back(LbDecayVertexMC->error().czx());
		   bDecayVtxYZE->push_back(LbDecayVertexMC->error().czy());

		   VDecayVtxX->push_back( Lam0_vFit_vertex_noMC->position().x() );
		   VDecayVtxY->push_back( Lam0_vFit_vertex_noMC->position().y() );
		   VDecayVtxZ->push_back( Lam0_vFit_vertex_noMC->position().z() );
		   VDecayVtxXE->push_back( Lam0_vFit_vertex_noMC->error().cxx() );
		   VDecayVtxYE->push_back( Lam0_vFit_vertex_noMC->error().cyy() );
		   VDecayVtxZE->push_back( Lam0_vFit_vertex_noMC->error().czz() );
		   VDecayVtxXYE->push_back( Lam0_vFit_vertex_noMC->error().cyx() );
		   VDecayVtxXZE->push_back( Lam0_vFit_vertex_noMC->error().czx() );
		   VDecayVtxYZE->push_back( Lam0_vFit_vertex_noMC->error().czy() );

           // ********************* muon-trigger-machint**************** 
		   
		   const pat::Muon* muon1 = &(*iMuon1);
		   const pat::Muon* muon2 = &(*iMuon2);

		   int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		   
		   if (muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr) tri_JpsiTkTk_tmp = 1;
		   
		   tri_Dim25->push_back( tri_Dim25_tmp );	       
		   tri_JpsiTk->push_back( tri_JpsiTk_tmp );
           tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );

 	       // ************
		  
		   mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
		   mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
		   mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
		   mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
		   mu1PF->push_back(iMuon1->isPFMuon());
		   mu2PF->push_back(iMuon2->isPFMuon());
		   mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		   mu2loose->push_back(muon::isLooseMuon(*iMuon2));

           //Trigger Selector
           drTrg_m1->push_back(dRMuonMatching1);
	       drTrg_m2->push_back(dRMuonMatching2);
		   std::cout << "pushing " << dRMuonMatching1 << " & " << dRMuonMatching2 << std::endl;

		   mumC2->push_back( glbTrackM->normalizedChi2() );
		   mumNHits->push_back( glbTrackM->numberOfValidHits() );
		   mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		   mupC2->push_back( glbTrackP->normalizedChi2() );
		   mupNHits->push_back( glbTrackP->numberOfValidHits() );
		   mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
           mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
		   mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
		   mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
		   mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
		   muon_dca->push_back(dca);
		   trg_dzm1->push_back(dzm1_trg);
		   trg_dzm2->push_back(dzm2_trg);

		   pi1dxy->push_back(v0daughters[proton_I].dxy());
		   pi2dxy->push_back(v0daughters[pion_I].dxy());
		   pi1dz->push_back(v0daughters[proton_I].dz());
		   pi2dz->push_back(v0daughters[pion_I].dz());

		   pi1dxy_e->push_back(v0daughters[proton_I].dxyError());
		   pi2dxy_e->push_back(v0daughters[pion_I].dxyError());
		   pi1dz_e->push_back(v0daughters[proton_I].dzError());
		   pi2dz_e->push_back(v0daughters[pion_I].dzError());

		   // try refitting the primary without the tracks in the B reco candidate		   
		   //std::cout<< "pass all" << std::endl;
		   nB++;	       
		   
		   /////////////////////////////////////////////////
		   //pionParticles.clear();
		   muonParticles.clear();
		   //vFitMCParticles.clear();
		  

		   }
	     }
	}
    }
  }
   
   //fill the tree and clear the vectors
   if (nB > 0 || OnlyGen_) 
     {
       std::cout << "filling tree" << std::endl;
       tree_->Fill();
	   std::cout << "fill ok" << std::endl;
     }
   // *********

   nB = 0; nMu = 0;
   trigger = 0;
   if (!OnlyGen_){	
   		B_mass->clear(); 
		LB_mass->clear();   B_px->clear();    B_py->clear();    B_pz->clear();
   		B_Ks0_mass->clear();
		LB_L0_mass->clear(); B_Ks0_px->clear(); B_Ks0_py->clear(); B_Ks0_pz->clear();

   		B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   		B_Ks0_pt1->clear(); B_Ks0_px1->clear(); B_Ks0_py1->clear(); B_Ks0_pz1->clear(); B_Ks0_charge1->clear(); 
   		B_Ks0_pt2->clear(); B_Ks0_px2->clear(); B_Ks0_py2->clear(); B_Ks0_pz2->clear(); B_Ks0_charge2->clear(); 

   		B_Ks0_px1_track->clear(); B_Ks0_py1_track->clear(); B_Ks0_pz1_track->clear(); 
   		B_Ks0_px2_track->clear(); B_Ks0_py2_track->clear(); B_Ks0_pz2_track->clear(); 

   		B_J_pt1->clear();  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   		B_J_pt2->clear();  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   		B_Ks0_chi2->clear(); B_J_chi2->clear(); B_chi2->clear();
   		B_Prob->clear(); B_J_Prob->clear(); B_ks0_Prob->clear();

		// *********
   		nTks->clear(); priVtxCL->clear(); TrkIndex->clear();
		PVTriggDz->clear();   
   		priVtxX->clear(); priVtxY->clear(); priVtxZ->clear(); 
   		priVtxXE->clear(); priVtxYE->clear(); priVtxZE->clear(); priVtxCL->clear();
   		priVtxXYE->clear(); priVtxXZE->clear(); priVtxYZE->clear();   
   } 
   nVtx = 0;
   if (!OnlyGen_){
   		bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear(); 
   		bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear(); 
   		bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();  

   		VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
   		VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
   		VDecayVtxXYE->clear(); VDecayVtxXZE->clear(); VDecayVtxYZE->clear();

   		pi1dxy->clear(); pi2dxy->clear(); pi1dz->clear(); pi2dz->clear();
   		pi1dxy_e->clear(); pi2dxy_e->clear(); pi1dz_e->clear(); pi2dz_e->clear();

   		mumC2->clear();
   		mumNHits->clear(); mumNPHits->clear();
   		mupC2->clear();
   		mupNHits->clear(); mupNPHits->clear();
   		mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();
		trg_dzm1->clear(); trg_dzm2->clear();   

   		tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();

   		mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   		mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 
   		//Trigger Selector
   		drTrg_m1->clear();
   		drTrg_m2->clear();
   }
}

bool JPsiLam0_PVpa::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool JPsiLam0_PVpa::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}
bool JPsiLam0_PVpa::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}
double JPsiLam0_PVpa::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}
//Try to print mc Tree
std::string JPsiLam0_PVpa::printName(int pdgid){
    std::unordered_map<int, std::string> umap;
    umap[11]  = "e-";
    umap[-11] = "e+";
    umap[13]  = "mu-";
    umap[-13] = "mu+";
    umap[22]  = "gamma";
    umap[23]  = "Z";
    umap[12]  = "Ve";
    umap[14]  = "Vmu";
    umap[-14] = "-Vmu";
    umap[15]  = "tau-";
    umap[-15] = "tau+";
    umap[16]  = "Vtau";
    umap[-16] = "-Vtau";
    umap[111] = "pi0";
    umap[211] = "pi+";
    umap[-211]= "pi-";
	umap[310] = "Ks0";
	umap[443] = "Jpsi";
	umap[511] = "B0";
	umap[-511]= "B0-";
	umap[521] = "B+";
	umap[-521]= "B-";
	umap[311] = "K0";
	umap[130] = "Kl0";
	umap[321] = "K+";
	umap[-321]= "K-";
	umap[413] = "D*(2010)+";
	umap[-413]= "D*(2010)-";
	umap[421] = "D0";
	umap[-421]= "D0-";
	umap[411] = "D+";
	umap[-411]= "D-";
	umap[415] = "D*2(2460)+";
	umap[-415]= "D*2(2460)-";
	umap[423] = "D*(2007)0";
	umap[-423]= "D*(2007)0-"; 
	umap[10411] = "D*0(2400)+";
	umap[-10411] = "D*0(2400)-";
	umap[2212]  = "p";
	umap[-2212] = "p-";
	umap[2112]  = "n";
	umap[-2112] = "n-";
	umap[21] = "g";
	umap[20213] = "a1(1260)+";
	umap[-20213]= "a1(1260)-";


    std::string retstr;
    if (umap.find(pdgid) == umap.end()){
        retstr = std::to_string(pdgid);
    }
    else{
        retstr = umap.at(pdgid);
    }

    return retstr;
}
void JPsiLam0_PVpa::printMCtree(const reco::Candidate* mother, int indent=0){
    if (mother == NULL){
         std::cout << "end tree" << std::endl;
         return;
    }
    if (mother->numberOfDaughters() > 0){
        if(indent){
                std::cout << std::setw(indent) << " ";
        }
        std::cout << printName(mother->pdgId()) <<" has "<< mother->numberOfDaughters() <<" daughters " <<std::endl;
    }
    int extraIndent = 0;
    for (size_t i=0; i< mother->numberOfDaughters(); i++){
        const reco::Candidate * daughter = mother->daughter(i);
        if (mother->numberOfDaughters() > 0){
            if(indent){
                std::cout << std::setw(indent) << " ";
            }
            std::cout << " daugter "<< i+1 <<": "<<  printName(daughter->pdgId()) << " with Pt: ";
            std::cout << daughter->pt() << " | Eta: "<< daughter->eta() << " | Phi: "<< daughter->phi() << " | mass: "<< daughter->mass() << std::endl;
            extraIndent+=4;
        }
        if (daughter->numberOfDaughters()) printMCtree(daughter, indent+extraIndent);
    }
}

void JPsiLam0_PVpa::printMCtreeUP(const reco::Candidate* daughter, int indent = 0){
    if (daughter == NULL){
        std::cout << "end tree" << std::endl;
        return;
    }
    int extraIndent = 0;
    for(size_t i = 0; i < daughter->numberOfMothers(); i++){
        const reco::Candidate* mother = daughter->mother(i);
		if (mother->numberOfMothers() > 0){
        	if(indent){
        	    std::cout << std::setw(indent) << " ";
        	}
        	std::cout<< "mother "<< i+1 << ": "<<printName(mother->pdgId()) << std::endl;
        	extraIndent+=1;
	    }
		if (mother->pdgId() == daughter->pdgId()) return;
        //if (mother->numberOfMothers() > 1) printMCtreeUP(daughter, indent+extraIndent);
		return;
    }
}
// ------------ method called once each job just before starting event loop  ------------

void 
JPsiLam0_PVpa::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->mu mu Ks0 ntuple");
  if (!OnlyGen_){
     tree_->Branch("nB",&nB,"nB/i");
     tree_->Branch("nMu",&nMu,"nMu/i");
   
     tree_->Branch("B_mass", &B_mass);
     tree_->Branch("LB_mass", &LB_mass);
     tree_->Branch("B_px", &B_px);
     tree_->Branch("B_py", &B_py);
     tree_->Branch("B_pz", &B_pz);
   
     tree_->Branch("B_Ks0_mass", &B_Ks0_mass);
     tree_->Branch("LB_L0_mass", &LB_L0_mass);
     tree_->Branch("B_Ks0_px", &B_Ks0_px);
     tree_->Branch("B_Ks0_py", &B_Ks0_py);
     tree_->Branch("B_Ks0_pz", &B_Ks0_pz);
    
     tree_->Branch("B_J_mass", &B_J_mass);
     tree_->Branch("B_J_px", &B_J_px);
     tree_->Branch("B_J_py", &B_J_py);
     tree_->Branch("B_J_pz", &B_J_pz);
   
     tree_->Branch("B_Ks0_pt1", &B_Ks0_pt1);
     tree_->Branch("B_Ks0_px1", &B_Ks0_px1);
     tree_->Branch("B_Ks0_py1", &B_Ks0_py1);
     tree_->Branch("B_Ks0_pz1", &B_Ks0_pz1);
     tree_->Branch("B_Ks0_px1_track", &B_Ks0_px1_track);
     tree_->Branch("B_Ks0_py1_track", &B_Ks0_py1_track);
     tree_->Branch("B_Ks0_pz1_track", &B_Ks0_pz1_track);
     tree_->Branch("B_Ks0_charge1", &B_Ks0_charge1); 
    
     tree_->Branch("B_Ks0_pt2", &B_Ks0_pt2);
     tree_->Branch("B_Ks0_px2", &B_Ks0_px2);
     tree_->Branch("B_Ks0_py2", &B_Ks0_py2);
     tree_->Branch("B_Ks0_pz2", &B_Ks0_pz2);
     tree_->Branch("B_Ks0_px2_track", &B_Ks0_px2_track);
     tree_->Branch("B_Ks0_py2_track", &B_Ks0_py2_track);
     tree_->Branch("B_Ks0_pz2_track", &B_Ks0_pz2_track);
     tree_->Branch("B_Ks0_charge2", &B_Ks0_charge2);
   
     tree_->Branch("B_J_pt1", &B_J_pt1);
     tree_->Branch("B_J_px1", &B_J_px1);
     tree_->Branch("B_J_py1", &B_J_py1);
     tree_->Branch("B_J_pz1", &B_J_pz1);
     tree_->Branch("B_J_charge1", &B_J_charge1);
   
     tree_->Branch("B_J_pt2", &B_J_pt2);
     tree_->Branch("B_J_px2", &B_J_px2);
     tree_->Branch("B_J_py2", &B_J_py2);
     tree_->Branch("B_J_pz2", &B_J_pz2);
     tree_->Branch("B_J_charge2", &B_J_charge2);
   
     tree_->Branch("B_chi2", &B_chi2);
     tree_->Branch("B_Ks0_chi2", &B_Ks0_chi2);
     tree_->Branch("B_J_chi2", &B_J_chi2);
   
     tree_->Branch("B_Prob",    &B_Prob);
     tree_->Branch("B_ks0_Prob", &B_ks0_Prob);
     tree_->Branch("B_J_Prob",  &B_J_Prob);
          
     // *************************
   
     tree_->Branch("priVtxX",&priVtxX);
     tree_->Branch("priVtxY",&priVtxY);
     tree_->Branch("priVtxZ",&priVtxZ);
     tree_->Branch("priVtxXE",&priVtxXE);
     tree_->Branch("priVtxYE",&priVtxYE);
     tree_->Branch("priVtxZE",&priVtxZE);
     tree_->Branch("priVtxXYE",&priVtxXYE);
     tree_->Branch("priVtxXZE",&priVtxXZE);
     tree_->Branch("priVtxYZE",&priVtxYZE);
     tree_->Branch("priVtxCL",&priVtxCL);
     tree_->Branch("nTks",&nTks);
     tree_->Branch("TrkIndex",&TrkIndex);
     tree_->Branch("PVTriggDz",&PVTriggDz);
	 
   
     tree_->Branch("nVtx",       &nVtx);
     tree_->Branch("run",        &run,       "run/I");
     tree_->Branch("event",        &event,     "event/I");
     tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
	 tree_->Branch("trigger",  &trigger, "trigger/i");
   
     tree_->Branch("bDecayVtxX",&bDecayVtxX);
     tree_->Branch("bDecayVtxY",&bDecayVtxY);
     tree_->Branch("bDecayVtxZ",&bDecayVtxZ);
     tree_->Branch("bDecayVtxXE",&bDecayVtxXE);
     tree_->Branch("bDecayVtxYE",&bDecayVtxYE);
     tree_->Branch("bDecayVtxZE",&bDecayVtxZE);
     tree_->Branch("bDecayVtxXYE",&bDecayVtxXYE);
     tree_->Branch("bDecayVtxXZE",&bDecayVtxXZE);
     tree_->Branch("bDecayVtxYZE",&bDecayVtxYZE);
   
     tree_->Branch("VDecayVtxX",&VDecayVtxX);
     tree_->Branch("VDecayVtxY",&VDecayVtxY);
     tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
     tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
     tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
     tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
     tree_->Branch("VDecayVtxXYE",&VDecayVtxXYE);
     tree_->Branch("VDecayVtxXZE",&VDecayVtxXZE);
     tree_->Branch("VDecayVtxYZE",&VDecayVtxYZE);
   
     tree_->Branch("pi1dxy",&pi1dxy);
     tree_->Branch("pi2dxy",&pi2dxy);
     tree_->Branch("pi1dz",&pi1dz);
     tree_->Branch("pi2dz",&pi2dz);
   
     tree_->Branch("pi1dxy_e",&pi1dxy_e);
     tree_->Branch("pi2dxy_e",&pi2dxy_e);
     tree_->Branch("pi1dz_e",&pi1dz_e);
     tree_->Branch("pi2dz_e",&pi2dz_e);
   
     tree_->Branch("mumC2",&mumC2);  
     tree_->Branch("mumNHits",&mumNHits);
     tree_->Branch("mumNPHits",&mumNPHits);
     tree_->Branch("mupC2",&mupC2);  
     tree_->Branch("mupNHits",&mupNHits);
     tree_->Branch("mupNPHits",&mupNPHits);
     tree_->Branch("mumdxy",&mumdxy);
     tree_->Branch("mupdxy",&mupdxy);
     tree_->Branch("muon_dca",&muon_dca);
     tree_->Branch("trg_dzm1",&trg_dzm1);
     tree_->Branch("trg_dzm2",&trg_dzm2);
   
     tree_->Branch("tri_Dim25",&tri_Dim25);
     tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
     tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk); 
    
     tree_->Branch("mumdz",&mumdz);
     tree_->Branch("mupdz",&mupdz);
     tree_->Branch("mu1soft",&mu1soft);
     tree_->Branch("mu2soft",&mu2soft);
     tree_->Branch("mu1tight",&mu1tight);
     tree_->Branch("mu2tight",&mu2tight);
     tree_->Branch("mu1PF",&mu1PF);
     tree_->Branch("mu2PF",&mu2PF);
     tree_->Branch("mu1loose",&mu1loose);
     tree_->Branch("mu2loose",&mu2loose);
	 //Trigger Selector
	 tree_->Branch("drTrg_m1", &drTrg_m1);
	 tree_->Branch("drTrg_m2", &drTrg_m2);
  }
    // gen
  if (isMC_ || OnlyGen_) {

     tree_->Branch("gen_b_p4",      "TLorentzVector",  &gen_b_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_ks0_p4" ,   "TLorentzVector",  &gen_ks0_p4);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_pion1_p4",  "TLorentzVector",  &gen_pion1_p4);
     tree_->Branch("gen_pion2_p4",  "TLorentzVector",  &gen_pion2_p4);
     tree_->Branch("gen_b_vtx",     "TVector3",        &gen_b_vtx);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_ks0_vtx",   "TVector3",        &gen_ks0_vtx);
     tree_->Branch("gen_b_ct",      &gen_b_ct,        "gen_b_ct/F");
     tree_->Branch("gen_ks0_ct",    &gen_ks0_ct,      "gen_ks0_ct/F");
  }

}


// ------------ method called once each job just after ending the event loop  ------------
void JPsiLam0_PVpa::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiLam0_PVpa);

