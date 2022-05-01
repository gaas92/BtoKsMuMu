// -*- C++ -*-
//
// Package:   JPsiKs0_onlyGen
// Class:     JPsiKs0_onlyGen
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  November 2020                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>

// user include files
#include "myAnalyzers/BtoKsMuMu/src/JPsiKs0_onlyGen.h"

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
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

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
JPsiKs0_onlyGen::JPsiKs0_onlyGen(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secundaryVerticesPtr"))),	
  Lost_track_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Lost_Tracks"))), 
  pPFC_track_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pPFC_Tracks"))), 
  triggerPrescalesSrc_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),

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

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

JPsiKs0_onlyGen::~JPsiKs0_onlyGen()
{

}

// ------------ method called to for each event  ------------
void JPsiKs0_onlyGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using std::vector;
  	using namespace edm;
  	using namespace reco;
  	using namespace std;
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
  	
	// for resonant in jpsi 
  	if (pruned.isValid() && isRes_) {   
		std::cout << "Res Pruned valid? " << pruned.isValid() << std::endl;  
    	int foundit = 0;
    	for (size_t i=0; i<pruned->size(); i++) {
      		foundit = 0;
      		const reco::Candidate *dau = &(*pruned)[i];
      		if ( (abs(dau->pdgId()) == 511) ) { //&& (dau->status() == 2) ) { //B0 
	    		foundit++;
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
					if (gdau->pdgId()==310){ //is K0s
			  			foundit++;
			  			gen_ks0_vtx.SetXYZ(gdau->vx(), gdau->vy(), gdau->vz());
			  			gen_ks0_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
			  			gen_ks0_ct = GetLifetime(gen_ks0_p4, gen_jpsi_vtx, gen_ks0_vtx); 
			  			// TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx 
		    		}// end if K0s
				}// end of B daughters for Ks0
      		} // end in B0
      		if (foundit>=5){
				ngen++;
	  		}  //1-B0, 2-JPsi, 3-mu1, 4-mu2, 5-Ks0,// NO PIONS 6-pi1, 7-pi2
    	} // for i
	}//end valid pruned 
	std::cout << "only gen, test ok"  << std::endl;

	tree_->Fill();

}

bool JPsiKs0_onlyGen::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool JPsiKs0_onlyGen::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool JPsiKs0_onlyGen::IsTheSame(const pat::GenericParticle& tk1, const pat::GenericParticle& tk2){
  double DeltaEta = fabs(tk2.eta()-tk1.eta());
  double DeltaP   = fabs(tk2.p()-tk1.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool JPsiKs0_onlyGen::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}
bool JPsiKs0_onlyGen::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}
double JPsiKs0_onlyGen::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}
//Try to print mc Tree
std::string JPsiKs0_onlyGen::printName(int pdgid){
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
void JPsiKs0_onlyGen::printMCtree(const reco::Candidate* mother, int indent=0){
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

void JPsiKs0_onlyGen::printMCtreeUP(const reco::Candidate* daughter, int indent = 0){
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
JPsiKs0_onlyGen::beginJob()
{

  std::cout << "Beginning analyzer " << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->mu mu Ks0 ntuple");
    // gen
  //if (isMC_ || OnlyGen_) {
  if (1) {
	 std::cout << "Definition of only gen TTree" << std::endl; 
	 tree_->Branch("ngen", &ngen, "ngen/I");
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
	 std::cout << "Definition of only gen TTree ok ..." << std::endl; 

  }

}


// ------------ method called once each job just after ending the event loop  ------------
void JPsiKs0_onlyGen::endJob() {
  std::cout << "Ending job" << std::endl;
  tree_->GetDirectory()->cd();
  tree_->Write();
  std::cout << "Ending job ok " << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKs0_onlyGen);

