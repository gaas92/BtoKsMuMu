// -*- C++ -*-
//
// Package:    Analyze/McGenAnalyzer
// Class:      McGenAnalyzer
//
/**\class McGenAnalyzer McGenAnalyzer.cc Analyze/McGenAnalyzer/plugins/McGenAnalyzer.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Horacio Crotte Ledesma
//         Created:  Wed, 02 Sep 2020 19:00:43 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

/*
HCL
 I will try to follow all information given here:
 https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface
*/
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//WE NEED THESE ONES FOR MAKING THE NTUPLES
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <utility>
#include <string>
#include "Math/GenVector/Boost.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include <Math/VectorUtil.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//
// FROM JHOVANNYS CODE
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

using reco::TrackCollection;

//class McGenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
class McGenAnalyzer : public edm::EDAnalyzer {
   public:
      explicit McGenAnalyzer(const edm::ParameterSet&);
      ~McGenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      bool    isAncestor(int, const reco::Candidate*);
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genCands_;
      
      TLorentzVector gen_b_p4,gen_jpsi_p4,gen_pion1_p4,gen_pion2_p4,gen_ks0_p4,gen_muon1_p4,gen_muon2_p4;
      TVector3       gen_b_vtx,gen_jpsi_vtx,gen_ks0_vtx;

      TTree *tree_;
      UInt_t    run;
      ULong64_t event;


      float Cos_T_L1;
      float Cos_T_L2;
      float Cos_T_LL; 
};


//
//genCands_(consumes<rec::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))), constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
McGenAnalyzer::McGenAnalyzer(const edm::ParameterSet& iConfig)
 //:number_daughters(0), number_daughtersJ(0), bplus(0), costhetaL(0.0), costhetaKL(0.0), costhetaLJ(0.0), costhetaKLJ(0.0)
{
  //std::cout << "INITIALIZER?" << std::endl;
  genCands_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  
}


McGenAnalyzer::~McGenAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//recursively check is a given particle has ancestor with given pdg_id
bool McGenAnalyzer::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}
//
// member functions
//

// ------------ method called for each event  ------------
void McGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //std::cout<< "frnklfnrsdf"<< std::endl;
  bool debug = true;

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

  Cos_T_L1  = 0;
  Cos_T_L2  = 0;
  Cos_T_LL  = 0; 


  if (debug) std::cout << "HELLO FROM ANALYZER! " << std::endl;
 
  //edm::Handle<reco::GenParticleCollection> pruned;
  edm::Handle<std::vector<reco::GenParticle>> pruned; 
  iEvent.getByToken(genCands_, pruned);

  //JHOVANNYS
  if (debug) std::cout << "PRUNED? " << pruned.isValid() std::endl;
  if ( pruned.isValid() ) {
    std::cout << "SIZE = " << pruned->size() << std::endl;
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
			   
		    }// end if K0s
	    }// end for B daughters for Ks0
      } // end if B0
      if (foundit>=5) {
        // Aqui creemos el boost al CM del dilepton
        math::XYZTLorentzVector gen_dilep(gen_jpsi_p4.Px(), gen_jpsi_p4.Py(), gen_jpsi_p4.Pz(), gen_jpsi_p4.M());
        ROOT::Math::Boost gen_cmboost(gen_dilep.BoostToCM());

        math::XYZTLorentzVector gen_kaon(gen_ks0_p4.Px(), gen_ks0_p4.Py(), gen_ks0_p4.Pz(), gen_ks0_p4.M());
        math::XYZTLorentzVector  gen_kaonCM(  gen_cmboost( gen_kaon )  );
        math::XYZTLorentzVector  gen_muonCMP, gen_muonCMN;

        math::XYZTLorentzVector gen_muon1(gen_muon1_p4.Px(), gen_muon1_p4.Py(), gen_muon1_p4.Pz(), gen_muon1_p4.M());
         math::XYZTLorentzVector gen_muon2(gen_muon2_p4.Px(), gen_muon2_p4.Py(), gen_muon2_p4.Pz(), gen_muon2_p4.M());

        gen_muonCMP = gen_cmboost(gen_muon2);
        gen_muonCMN = gen_cmboost(gen_muon1);
        
        Cos_T_LL = ( gen_muonCMP.x()*gen_muonCMN.x() 
                   + gen_muonCMP.y()*gen_muonCMN.y() 
                   + gen_muonCMP.z()*gen_muonCMN.z() ) / (gen_muonCMP.P()*gen_muonCMN.P() );

        Cos_T_L1 = ( gen_muonCMP.x()*gen_kaonCM.x() 
                   + gen_muonCMP.y()*gen_kaonCM.y() 
                   + gen_muonCMP.z()*gen_kaonCM.z() ) / (gen_muonCMP.P()*gen_kaonCM.P() );
        Cos_T_L2 = ( gen_muonCMN.x()*gen_kaonCM.x() 
                   + gen_muonCMN.y()*gen_kaonCM.y() 
                   + gen_muonCMN.z()*gen_kaonCM.z() ) / (gen_muonCMN.P()*gen_kaonCM.P() );

        tree_->Fill();
        break; //1-B0, 2-JPsi, 3-mu1, 4-mu2, 5-Ks0 //NOPIONS , 6-pi1, 7-pi2
      }
    }// end for pruned  
  }// end if pruned
}// end analize


// ------------ method called once each job just before starting event loop  ------------
void
McGenAnalyzer::beginJob()
{
  //std::cout << "Beginning analyzer job" << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make < TTree > ("tree", "Tree of gen B0");

  tree_->Branch("run",      &run,      "run/i");
  tree_->Branch("event",    &event,    "event/i");

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
  tree_->Branch("Cos_T_L1",    &Cos_T_L1,      "Cos_T_L1/F");
  tree_->Branch("Cos_T_L2",    &Cos_T_L2,      "Cos_T_L2/F");
  tree_->Branch("Cos_T_LL",    &Cos_T_LL,      "Cos_T_LL/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void
McGenAnalyzer::endJob()
{
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
McGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(McGenAnalyzer);