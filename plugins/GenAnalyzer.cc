// -*- C++ -*-
//
// Package:    Gen/GenAnalyzer
// Class:      GenAnalyzer
//
/**\class GenAnalyzer GenAnalyzer.cc Gen/GenAnalyzer/plugins/GenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bhim Bam
//         Created:  Sun, 16 Jul 2023 18:47:17 GMT
//
//


// system include files
#include <memory>
// #include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//for the GenParticleCollection and GenParticles
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TTree.h"
//
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
//


// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



// int pdgid_ = 25; // A/H to tautau
int pdgid_ = 553; // Upsilon 1S to tauatau


using std::vector;
using reco::GenParticle;

int ntotal_event ;
int npassed_event ;

unsigned int runId_;
unsigned int lumiId_;
unsigned long long eventId_;

TH1D *H_tau_att_genA1_M_inv;
TH1D *H_tau_att_genA1_M;
TH1D *H_tau_att_dR_A1_Tau1;
TH1D *H_tau_att_dR_A1_Tau2;
TH1D *H_tau_att_dR_Tau1_Tau2;

TH1D *H_tau_att_A1_pt;
TH1D *H_tau_att_Tau1_pt;
TH1D *H_tau_att_Tau2_pt;
TH1D *H_tau_att_A1_eta;

TH1D *H_tau_att_Tau1_eta;
TH1D *H_tau_att_Tau2_eta;
TH1D *H_tau_att_A1_phi;
TH1D *H_tau_att_Tau1_phi;
TH1D *H_tau_att_Tau2_phi;

TH1D *H_tau_att_Tau1_Tau2_deta;
TH1D *H_tau_att_Tau1_Tau2_dphi;

TH2D *H_tau_att_Tau1_Tau2_dphi_deta;

float V_att_genA1_M_inv;
float V_att_genA1_M;
float V_att_dR_A1_Tau1;
float V_att_dR_A1_Tau2;
float V_att_dR_Tau1_Tau2;

float V_att_A1_pt;
float V_att_Tau1_pt;
float V_att_Tau2_pt;
float V_att_A1_eta;
float V_att_Tau1_eta;
float V_att_Tau2_eta;
float V_att_A1_phi;
float V_att_Tau1_phi;
float V_att_Tau2_phi;

float V_att_Tau1_Tau2_deta;
float V_att_Tau1_Tau2_dphi;

// vector<int> vAIdxs;
vector<float> V_att_genA1_M_inv_;
vector<float> V_att_genA1_M_;
vector<float> V_att_dR_A1_Tau1_;
vector<float> V_att_dR_A1_Tau2_;
vector<float> V_att_dR_Tau1_Tau2_;

vector<float> V_att_A1_pt_;
vector<float> V_att_Tau1_pt_;
vector<float> V_att_Tau2_pt_;
vector<float> V_att_A1_eta_;
vector<float> V_att_Tau1_eta_;
vector<float> V_att_Tau2_eta_;
vector<float> V_att_A1_phi_;
vector<float> V_att_Tau1_phi_;
vector<float> V_att_Tau2_phi_;

vector<float> V_att_Tau1_Tau2_deta_;
vector<float> V_att_Tau1_Tau2_dphi_;


TLorentzVector SetTaus(Float_t tau_pt, Float_t tau_eta, Float_t tau_phi, Float_t tau_mass){
  TLorentzVector Tau_Candidate;
  Tau_Candidate.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_mass);
  return Tau_Candidate;
}

// TLorentzVector SetMETs(Float_t met, Float_t metphi){
//   TLorentzVector Met;
//   Met.SetPtEtaPhiM(met, 0, metphi, 0.);
//   return Met;
// }

class GenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenAnalyzer(const edm::ParameterSet&);
      ~GenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      // edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::Service<TFileService> fs;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_;
      edm::InputTag genParticles_;

      TTree *RHTree;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
 // :
  // tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
   //now do what ever initialization is needed
   RHTree = fs->make<TTree>("RHTree","Gen info Tree");

   H_tau_att_genA1_M_inv     = fs->make<TH1D>("h_genA1_M_inv"   , "m^{gen_inv A1};m^{gen_inv A1};Events"                     ,  30,  3, 15);
   H_tau_att_genA1_M     = fs->make<TH1D>("h_genA1_M"   , "m^{gen A1};m^{gen A1};Events"                     ,  30,  3, 15);
   H_tau_att_dR_A1_Tau1     = fs->make<TH1D>("h_dR_A1_Tau1"   , "dR^{gen A1_Tau1};dR^{gen A1_Tau1};Events"                     ,  10,  0, 5);
   H_tau_att_dR_A1_Tau2     = fs->make<TH1D>("h_dR_A1_Tau2"   , "dR^{gen A1_Tau2};dR^{gen A1_Tau2};Events"                     ,  10,  0, 5);
   H_tau_att_dR_Tau1_Tau2     = fs->make<TH1D>("h_dR_Tau1_Tau2"   , "dR^{gen Tau1_Tau2};dR^{gen Tau1_Tau2};Events"                     ,  10,  0, .6);
   H_tau_att_A1_pt     = fs->make<TH1D>("h_A1_pt"   , "pt^{gen A1};pt^{gen A1};Events"                     ,  100,  0, 200);
   H_tau_att_Tau1_pt     = fs->make<TH1D>("h_Tau1_pt"   , "pt^{gen Tau1};pt^{gen Tau1};Events"                     ,  70,  10, 150);
   H_tau_att_Tau2_pt     = fs->make<TH1D>("h_Tau2_pt"   , "pt^{gen Tau2};pt^{gen Tau2};Events"                     ,  70,  10, 150);
   H_tau_att_A1_eta     = fs->make<TH1D>("h_A1_eta"   , "eta^{gen A1};eta^{gen A1};Events"                     ,  20,  -5, 5);
   H_tau_att_Tau1_eta     = fs->make<TH1D>("h_Tau1_eta"   , "eta^{gen Tau1};eta^{gen Tau1};Events"                     ,  20,  -5, 5);
   H_tau_att_Tau2_eta     = fs->make<TH1D>("h_Tau2_eta"   , "eta^{gen Tau2};eta^{gen Tau2};Events"                     ,  20,  -5, 5);
   H_tau_att_A1_phi     = fs->make<TH1D>("h_A1_phi"   , "phi^{gen A1};phi^{gen A1};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau1_phi     = fs->make<TH1D>("h_Tau1_phi"   , "phi^{gen Tau1};phi^{gen Tau1};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau2_phi     = fs->make<TH1D>("h_Tau2_phi"   , "phi^{gen Tau2};phi^{gen Tau2};Events"                     ,  20,  -3.2, 3.2);

   H_tau_att_Tau1_Tau2_deta     = fs->make<TH1D>("h_Tau1_Tau2_deta"   , "deta^{gen Tau1_Tau2};deta^{gen Tau1_Tau2};Events"                     ,  20,  0, .8);
   H_tau_att_Tau1_Tau2_dphi     = fs->make<TH1D>("h_Tau1_Tau2_dphi"   , "dphi^{gen Tau1_Tau2};dphi^{gen Tau1_Tau2};Events"                     ,  20,  0, .8);

   H_tau_att_Tau1_Tau2_dphi_deta     = fs->make<TH2D>("h_Tau1_Tau2_dphi_deta"   , "dphi vs deta^{gen Tau1_Tau2};dphi;deta" , 20,  0, .5 , 20,  0, .5);

   RHTree->Branch("Event",  &eventId_);
   RHTree->Branch("Run",  &runId_);
   RHTree->Branch("LumiSection",  &lumiId_);

   RHTree->Branch("GenA1_inv",  &V_att_genA1_M_inv);
   RHTree->Branch("GenA1",  &V_att_genA1_M);
   RHTree->Branch("dR_A1_Tau1",  &V_att_dR_A1_Tau1);
   RHTree->Branch("dR_A1_Tau2",  &V_att_dR_A1_Tau2);
   RHTree->Branch("dR_Tau1_Tau2",  &V_att_dR_Tau1_Tau2);

   RHTree->Branch("A1_pt",  &V_att_A1_pt);
   RHTree->Branch("Tau1_pt",  &V_att_Tau1_pt);
   RHTree->Branch("Tau2_pt",  &V_att_Tau2_pt);
   RHTree->Branch("A1_eta",  &V_att_A1_eta);
   RHTree->Branch("Tau1_eta",  &V_att_Tau1_eta);
   RHTree->Branch("Tau2_eta",  &V_att_Tau2_eta);
   RHTree->Branch("A1_phi",  &V_att_A1_phi);
   RHTree->Branch("Tau1_phi",  &V_att_Tau1_phi);
   RHTree->Branch("Tau2_phi",  &V_att_Tau2_phi);

   RHTree->Branch("Tau1_Tau2_deta",  &V_att_Tau1_Tau2_deta);
   RHTree->Branch("Tau1_Tau2_dphi",  &V_att_Tau1_Tau2_dphi);

   genParticlesToken_   = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));

}


GenAnalyzer::~GenAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  eventId_ = iEvent.id().event();
  runId_ = iEvent.id().run();
  lumiId_ = iEvent.id().luminosityBlock();
  // std::cout<<"Event info"<<eventId_<<":  "<<runId_<<":  "<<lumiId_<<std::endl;

   V_att_genA1_M_inv    = -1111.1111;
   V_att_genA1_M    = -1111.1111;
   V_att_dR_A1_Tau1    = -1111.1111;
   V_att_dR_A1_Tau2    = -1111.1111;
   V_att_dR_Tau1_Tau2    = -1111.1111;
   V_att_A1_pt    = -1111.1111;
   V_att_Tau1_pt    = -1111.1111;
   V_att_Tau2_pt    = -1111.1111;
   V_att_A1_eta    = -1111.1111;
   V_att_Tau1_eta    = -1111.1111;
   V_att_Tau2_eta    = -1111.1111;
   V_att_A1_phi    = -1111.1111;
   V_att_Tau1_phi    = -1111.1111;
   V_att_Tau2_phi    = -1111.1111;

   V_att_Tau1_Tau2_deta    = -1111.1111;
   V_att_Tau1_Tau2_dphi    = -1111.1111;

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken(genParticlesToken_,   genParticles);


  // unsigned int NAs = 0;
  // // unsigned int NTau_fromA = 0;
  // // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  // //   if ( abs(iGen->pdgId()) != 15 || abs(iGen->mother()->pdgId()) != 25) continue;
  // //   NTau_fromA++;
  // // }
  // // std::cout << "  >>>>>> Number Tau from  A <<<<<"<<"    "  <<  NTau_fromA << std::endl;
  // // vAIdxs.clear();
  // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  //   if ( iGen->pdgId() != -25 || abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15) continue;
  //   NAs++;
  //
  //   // vAIdxs.push_back(NAs-1);
  //   // std::cout<<"Size------------"<<vAIdxs.size()<<std::endl;
  // }
  // std::cout<<" Anti particle "<<NAs<<std::endl;
  // std::cout << "  >>>>>> Number of A giving Tau <<<<<"<<"    "  <<  NAs << std::endl;

float genA1_mass_inv = -1111.1111;
float genA1_mass = -1111.1111;
float A1_Tau1_dR = -1111.1111;
float A1_Tau2_dR = -1111.1111;
float Tau1_Tau2_dR = -1111.1111;
float A1_pt = -1111.1111;
float Tau1_pt = -1111.1111;
float Tau2_pt = -1111.1111;
float A1_eta = -1111.1111;
float Tau1_eta = -1111.1111;
float Tau2_eta = -1111.1111;
float A1_phi = -1111.1111;
float Tau1_phi = -1111.1111;
float Tau2_phi = -1111.1111;

float Tau1_Tau2_deta = -1111.1111;
float Tau1_Tau2_dphi = -1111.1111;



for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) { //Gen loop
  bool pass = false;
  if (abs(iGen->pdgId()) != pdgid_ || iGen->numberOfDaughters() != 2 || abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;
  ntotal_event++;
  if ( abs(iGen->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->status()) != 2 || iGen->daughter(0)->numberOfMothers() < 1 || iGen->daughter(1)->numberOfMothers() < 1) continue;


  float dR_A1_Tau1 = reco::deltaR( iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->eta(), iGen->phi());
  float dR_A1_Tau2 = reco::deltaR( iGen->daughter(1)->eta(), iGen->daughter(1)->phi(), iGen->eta(), iGen->phi());
  float dR_Tau1_Tau2 = reco::deltaR( iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());

  if(dR_Tau1_Tau2 > 0.4) continue;
  pass = true;


  TLorentzVector GenTau1  = SetTaus(iGen->daughter(0)->pt(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(0)->mass());
  TLorentzVector GenTau2  = SetTaus(iGen->daughter(1)->pt(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi(), iGen->daughter(1)->mass());
  TLorentzVector GenA1 = GenTau1 + GenTau2;


  genA1_mass_inv = GenA1.M();

  genA1_mass = iGen->mass();




  A1_Tau1_dR = dR_A1_Tau1;
  A1_Tau2_dR = dR_A1_Tau2;
  Tau1_Tau2_dR = dR_Tau1_Tau2;

  A1_pt = iGen->pt();
  Tau1_pt = iGen->daughter(0)->pt();
  Tau2_pt = iGen->daughter(1)->pt();
  A1_eta = iGen->eta();
  Tau1_eta = iGen->daughter(0)->eta();
  Tau2_eta = iGen->daughter(1)->eta();
  A1_phi = iGen->phi();
  Tau1_phi = iGen->daughter(0)->phi();
  Tau2_phi = iGen->daughter(1)->phi();

  Tau1_Tau2_deta = abs(Tau1_eta-Tau2_eta);
  Tau1_Tau2_dphi = abs(Tau1_phi-Tau2_phi);


  std::cout << "  >>>>>> A (25)/ U (553) gen  <<<<<"<<"<<< status: "<<iGen->status()<<"<<<pt:  "<<iGen->pt()<<"<<<eta:  "<<iGen->eta()<<"<<<phi:  "<<iGen->phi()<<"<<<mass:  "<<iGen->mass() <<"<<< pdgid:  "<<iGen->pdgId()<< std::endl;
  std::cout << "  >>>>>> Gen A (25)/ U (553) LorentzVector  <<<<<"<<"<<<mass:  "<< GenA1.M() << std::endl;
  std::cout << "  >>>>>> Tau1 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(0)->mass() << std::endl;
  std::cout << "  >>>>>> Tau2 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(1)->mass() << std::endl;
  std::cout << "  >>>>>> dR_A1_Tau1:  "<<dR_A1_Tau1<< std::endl;
  std::cout << "  >>>>>> dR_A1_Tau2:  "<<dR_A1_Tau2<< std::endl;
  std::cout << "  >>>>>> dR_Tau1_Tau2:  "<<dR_Tau1_Tau2<< std::endl;

  // }

if (pass){
npassed_event++;


V_att_genA1_M_inv    = genA1_mass_inv;
V_att_genA1_M    = genA1_mass;
V_att_dR_A1_Tau1    = A1_Tau1_dR;
V_att_dR_A1_Tau2    = A1_Tau2_dR;
V_att_dR_Tau1_Tau2    = Tau1_Tau2_dR;

V_att_A1_pt    = A1_pt;
V_att_Tau1_pt    = Tau1_pt;
V_att_Tau2_pt    = Tau2_pt;
V_att_A1_eta    = A1_eta;
V_att_Tau1_eta    = Tau1_eta;
V_att_Tau2_eta    = Tau2_eta;
V_att_A1_phi    = A1_phi;
V_att_Tau1_phi    = Tau1_phi;
V_att_Tau2_phi    = Tau2_phi;

V_att_Tau1_Tau2_deta    = Tau1_Tau2_deta;
V_att_Tau1_Tau2_dphi    = Tau1_Tau2_dphi;

V_att_genA1_M_inv_.push_back( V_att_genA1_M_inv );
V_att_genA1_M_.push_back( V_att_genA1_M );
V_att_dR_A1_Tau1_.push_back( V_att_dR_A1_Tau1 );
V_att_dR_A1_Tau2_.push_back( V_att_dR_A1_Tau2 );
V_att_dR_Tau1_Tau2_.push_back( V_att_dR_Tau1_Tau2 );

V_att_A1_pt_.push_back( V_att_A1_pt );
V_att_Tau1_pt_.push_back( V_att_Tau1_pt );
V_att_Tau2_pt_.push_back( V_att_Tau2_pt );
V_att_A1_eta_.push_back( V_att_A1_eta );
V_att_Tau1_eta_.push_back( V_att_Tau1_eta );
V_att_Tau2_eta_.push_back( V_att_Tau2_eta );
V_att_A1_phi_.push_back( V_att_A1_phi );
V_att_Tau1_phi_.push_back( V_att_Tau1_phi );
V_att_Tau2_phi_.push_back( V_att_Tau2_phi );

V_att_Tau1_Tau2_deta_.push_back( V_att_Tau1_Tau2_deta );
V_att_Tau1_Tau2_dphi_.push_back( V_att_Tau1_Tau2_dphi );




// // for ( unsigned iG(0); iG != vAIdxs.size(); ++iG ) {
//
H_tau_att_genA1_M_inv->Fill( V_att_genA1_M_inv);
H_tau_att_genA1_M->Fill( V_att_genA1_M);
H_tau_att_dR_A1_Tau1->Fill( V_att_dR_A1_Tau1);
H_tau_att_dR_A1_Tau2->Fill( V_att_dR_A1_Tau2);
H_tau_att_dR_Tau1_Tau2->Fill( V_att_dR_Tau1_Tau2);

H_tau_att_A1_pt->Fill( V_att_A1_pt);
H_tau_att_Tau1_pt->Fill( V_att_Tau1_pt);
H_tau_att_Tau2_pt->Fill( V_att_Tau2_pt);
H_tau_att_A1_eta->Fill( V_att_A1_eta);
H_tau_att_Tau1_eta->Fill( V_att_Tau1_eta);
H_tau_att_Tau2_eta->Fill( V_att_Tau2_eta);
H_tau_att_A1_phi->Fill( V_att_A1_phi);
H_tau_att_Tau1_phi->Fill( V_att_Tau1_phi);
H_tau_att_Tau2_phi->Fill( V_att_Tau2_phi);

H_tau_att_Tau1_Tau2_deta->Fill( V_att_Tau1_Tau2_deta);
H_tau_att_Tau1_Tau2_dphi->Fill( V_att_Tau1_Tau2_dphi);

H_tau_att_Tau1_Tau2_dphi_deta->Fill(V_att_Tau1_Tau2_dphi,V_att_Tau1_Tau2_deta);

RHTree->Fill();

}//pass loop
}//Gen loop


V_att_genA1_M_inv_.clear();
V_att_genA1_M_.clear();
V_att_dR_A1_Tau1_.clear();
V_att_dR_A1_Tau2_.clear();
V_att_dR_Tau1_Tau2_.clear();

V_att_A1_pt_.clear();
V_att_Tau1_pt_.clear();
V_att_Tau2_pt_.clear();
V_att_A1_eta_.clear();
V_att_Tau1_eta_.clear();
V_att_Tau2_eta_.clear();
V_att_A1_phi_.clear();
V_att_Tau1_phi_.clear();
V_att_Tau2_phi_.clear();

V_att_Tau1_Tau2_deta_.clear();
V_att_Tau1_Tau2_dphi_.clear();


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
GenAnalyzer::beginJob()
{
  ntotal_event = 0;
  npassed_event = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenAnalyzer::endJob()

{


  std::cout << "  >>>>>> Total events selected events <<<<<  "<<npassed_event<<"/"<<ntotal_event<<std::endl;


}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(GenAnalyzer);
