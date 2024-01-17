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
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
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
// This will improve performance in multithreaded jobs
bool debug = true;
using std::vector;
using reco::GenParticle;

int ntotal_event ;
int npassed_event ;

unsigned int runId_;
unsigned int lumiId_;
unsigned long long eventId_;


TH1D *H_tau_att_genHiggs_M_inv;
TH1D *H_tau_att_genA1_M_inv;
TH1D *H_tau_att_genA2_M_inv;
TH1D *H_tau_att_genHiggs_M;
TH1D *H_tau_att_genA1_M;
TH1D *H_tau_att_genA2_M;
TH1D *H_tau_att_dR_A1_A2;
TH1D *H_tau_att_dR_H_A1;
TH1D *H_tau_att_dR_H_A2;
TH1D *H_tau_att_dR_A1_Tau1;
TH1D *H_tau_att_dR_A1_Tau2;
TH1D *H_tau_att_dR_A2_Tau3;
TH1D *H_tau_att_dR_A2_Tau4;
TH1D *H_tau_att_dR_Tau1_Tau2;
TH1D *H_tau_att_dR_Tau3_Tau4;

TH1D *H_tau_att_H_pt;
TH1D *H_tau_att_A1_pt;
TH1D *H_tau_att_A2_pt;
TH1D *H_tau_att_Tau1_pt;
TH1D *H_tau_att_Tau2_pt;
TH1D *H_tau_att_Tau3_pt;
TH1D *H_tau_att_Tau4_pt;
TH1D *H_tau_att_H_eta;
TH1D *H_tau_att_A1_eta;
TH1D *H_tau_att_A2_eta;
TH1D *H_tau_att_Tau1_eta;
TH1D *H_tau_att_Tau2_eta;
TH1D *H_tau_att_Tau3_eta;
TH1D *H_tau_att_Tau4_eta;
TH1D *H_tau_att_H_phi;
TH1D *H_tau_att_A1_phi;
TH1D *H_tau_att_A2_phi;
TH1D *H_tau_att_Tau1_phi;
TH1D *H_tau_att_Tau2_phi;
TH1D *H_tau_att_Tau3_phi;
TH1D *H_tau_att_Tau4_phi;

TH1D *H_tau_att_Tau1_Tau2_deta;
TH1D *H_tau_att_Tau1_Tau2_dphi;
TH1D *H_tau_att_Tau3_Tau4_deta;
TH1D *H_tau_att_Tau3_Tau4_dphi;

TH2D *H_tau_att_Tau1_Tau2_dphi_deta;
TH2D *H_tau_att_Tau3_Tau4_dphi_deta;

float V_att_genHiggs_M_inv;
float V_att_genA1_M_inv;
float V_att_genA2_M_inv;
float V_att_genHiggs_M;
float V_att_genA1_M;
float V_att_genA2_M;
float V_att_dR_A1_A2;
float V_att_dR_H_A1;
float V_att_dR_H_A2;
float V_att_dR_A1_Tau1;
float V_att_dR_A1_Tau2;
float V_att_dR_A2_Tau3;
float V_att_dR_A2_Tau4;
float V_att_dR_Tau1_Tau2;
float V_att_dR_Tau3_Tau4;

float V_att_H_pt;
float V_att_A1_pt;
float V_att_A2_pt;
float V_att_Tau1_pt;
float V_att_Tau2_pt;
float V_att_Tau3_pt;
float V_att_Tau4_pt;
float V_att_H_eta;
float V_att_A1_eta;
float V_att_A2_eta;
float V_att_Tau1_eta;
float V_att_Tau2_eta;
float V_att_Tau3_eta;
float V_att_Tau4_eta;
float V_att_H_phi;
float V_att_A1_phi;
float V_att_A2_phi;
float V_att_Tau1_phi;
float V_att_Tau2_phi;
float V_att_Tau3_phi;
float V_att_Tau4_phi;

float V_att_Tau1_Tau2_deta;
float V_att_Tau1_Tau2_dphi;
float V_att_Tau3_Tau4_deta;
float V_att_Tau3_Tau4_dphi;


vector<float> V_att_genHiggs_M_inv_;
vector<float> V_att_genA1_M_inv_;
vector<float> V_att_genA2_M_inv_;
vector<float> V_att_genHiggs_M_;
vector<float> V_att_genA1_M_;
vector<float> V_att_genA2_M_;
vector<float> V_att_dR_A1_A2_;
vector<float> V_att_dR_H_A1_;
vector<float> V_att_dR_H_A2_;
vector<float> V_att_dR_A1_Tau1_;
vector<float> V_att_dR_A1_Tau2_;
vector<float> V_att_dR_A2_Tau3_;
vector<float> V_att_dR_A2_Tau4_;
vector<float> V_att_dR_Tau1_Tau2_;
vector<float> V_att_dR_Tau3_Tau4_;

vector<float> V_att_H_pt_;
vector<float> V_att_A1_pt_;
vector<float> V_att_A2_pt_;
vector<float> V_att_Tau1_pt_;
vector<float> V_att_Tau2_pt_;
vector<float> V_att_Tau3_pt_;
vector<float> V_att_Tau4_pt_;
vector<float> V_att_H_eta_;
vector<float> V_att_A1_eta_;
vector<float> V_att_A2_eta_;
vector<float> V_att_Tau1_eta_;
vector<float> V_att_Tau2_eta_;
vector<float> V_att_Tau3_eta_;
vector<float> V_att_Tau4_eta_;
vector<float> V_att_H_phi_;
vector<float> V_att_A1_phi_;
vector<float> V_att_A2_phi_;
vector<float> V_att_Tau1_phi_;
vector<float> V_att_Tau2_phi_;
vector<float> V_att_Tau3_phi_;
vector<float> V_att_Tau4_phi_;

vector<float> V_att_Tau1_Tau2_deta_;
vector<float> V_att_Tau1_Tau2_dphi_;
vector<float> V_att_Tau3_Tau4_deta_;
vector<float> V_att_Tau3_Tau4_dphi_;


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

      edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
      edm::InputTag genJets_;

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

   H_tau_att_genHiggs_M_inv     = fs->make<TH1D>("h_genHiggs_M_inv"   , "m^{gen_inv H};m^{gen_inv H};Events"                 ,  10,  120, 130);
   H_tau_att_genA1_M_inv     = fs->make<TH1D>("h_genA1_M_inv"   , "m^{gen_inv A1};m^{gen_inv A1};Events"                     ,  30,  3, 15);
   H_tau_att_genA2_M_inv     = fs->make<TH1D>("h_genA2_M_inv"   , "m^{gen_inv A2};m^{gen_inv A2};Events"                     ,  30,  3, 15);
   H_tau_att_genHiggs_M     = fs->make<TH1D>("h_genHiggs_M"   , "m^{gen H};m^{gen H};Events"                 ,  10,  120, 130);
   H_tau_att_genA1_M     = fs->make<TH1D>("h_genA1_M"   , "m^{gen A1};m^{gen A1};Events"                     ,  30,  3, 15);
   H_tau_att_genA2_M     = fs->make<TH1D>("h_genA2_M"   , "m^{gen A2};m^{gen A2};Events"                     ,  30,  3, 15);
   H_tau_att_dR_A1_A2     = fs->make<TH1D>("h_dR_A1_A2"   , "dR^{gen A1_A2};dR^{gen A1_A2};Events"                     ,  10,  0, 5);
   H_tau_att_dR_H_A1     = fs->make<TH1D>("h_dR_H_A1"   , "dR^{ gen H_A1};dR^{ gen H_A1};Events"                     ,  10,  0, 5);
   H_tau_att_dR_H_A2     = fs->make<TH1D>("h_dR_H_A2"   , "dR^{gen H_A2};dR^{gen H_A2};Events"                     ,  10,  0, 5);
   H_tau_att_dR_A1_Tau1     = fs->make<TH1D>("h_dR_A1_Tau1"   , "dR^{gen A1_Tau1};dR^{gen A1_Tau1};Events"                     ,  10,  0, 5);
   H_tau_att_dR_A1_Tau2     = fs->make<TH1D>("h_dR_A1_Tau2"   , "dR^{gen A1_Tau2};dR^{gen A1_Tau2};Events"                     ,  10,  0, 5);
   H_tau_att_dR_A2_Tau3     = fs->make<TH1D>("h_dR_A2_Tau3"   , "dR^{gen A2_Tau3};dR^{gen A2_Tau3};Events"                     ,  10,  0, 5);
   H_tau_att_dR_A2_Tau4     = fs->make<TH1D>("h_dR_A2_Tau4"   , "dR^{gen A2_Tau4};dR^{gen A2_Tau4};Events"                     ,  10,  0, 5);
   H_tau_att_dR_Tau1_Tau2     = fs->make<TH1D>("h_dR_Tau1_Tau2"   , "dR^{gen Tau1_Tau2};dR^{gen Tau1_Tau2};Events"                     ,  10,  0, .6);
   H_tau_att_dR_Tau3_Tau4     = fs->make<TH1D>("h_dR_Tau3_Tau4"   , "dR^{gen Tau43_Tau4};dR^{gen Tau43_Tau4};Events"                     ,  10,  0, .6);

   H_tau_att_H_pt     = fs->make<TH1D>("h_H_pt"   , "pt^{gen H};pt^{gen H};Events"                     ,  200,  0, 400);
   H_tau_att_A1_pt     = fs->make<TH1D>("h_A1_pt"   , "pt^{gen A1};pt^{gen A1};Events"                     ,  100,  0, 200);
   H_tau_att_A2_pt     = fs->make<TH1D>("h_A2_pt"   , "pt^{gen A2};pt^{gen A2};Events"                     ,  100,  0, 200);
   H_tau_att_Tau1_pt     = fs->make<TH1D>("h_Tau1_pt"   , "pt^{gen Tau1};pt^{gen Tau1};Events"                     ,  70,  10, 150);
   H_tau_att_Tau2_pt     = fs->make<TH1D>("h_Tau2_pt"   , "pt^{gen Tau2};pt^{gen Tau2};Events"                     ,  70,  10, 150);
   H_tau_att_Tau3_pt     = fs->make<TH1D>("h_Tau3_pt"   , "pt^{gen Tau3};pt^{gen Tau3};Events"                     ,  70,  10, 150);
   H_tau_att_Tau4_pt     = fs->make<TH1D>("h_Tau4_pt"   , "pt^{gen Tau4};pt^{gen Tau4};Events"                     ,  70,  10, 150);
   H_tau_att_H_eta     = fs->make<TH1D>("h_H_eta"   , "eta^{gen H};eta^{gen H};Events"                     ,  20,  -5, 5);
   H_tau_att_A1_eta     = fs->make<TH1D>("h_A1_eta"   , "eta^{gen A1};eta^{gen A1};Events"                     ,  20,  -5, 5);
   H_tau_att_A2_eta     = fs->make<TH1D>("h_A2_eta"   , "eta^{gen A2};eta^{gen A2};Events"                     ,  20,  -5, 5);
   H_tau_att_Tau1_eta     = fs->make<TH1D>("h_Tau1_eta"   , "eta^{gen Tau1};eta^{gen Tau1};Events"                     ,  20,  -5, 5);
   H_tau_att_Tau2_eta     = fs->make<TH1D>("h_Tau2_eta"   , "eta^{gen Tau2};eta^{gen Tau2};Events"                     ,  20,  -5, 5);
   H_tau_att_Tau3_eta     = fs->make<TH1D>("h_Tau3_eta"   , "eta^{gen Tau3};eta^{gen Tau3};Events"                     ,  20,  -5, 5);
   H_tau_att_Tau4_eta     = fs->make<TH1D>("h_Tau4_eta"   , "eta^{gen Tau4};eta^{gen Tau4};Events"                     ,  20,  -5, 5);
   H_tau_att_H_phi     = fs->make<TH1D>("h_H_phi"   , "phi^{gen H};phi^{gen H};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_A1_phi     = fs->make<TH1D>("h_A1_phi"   , "phi^{gen A1};phi^{gen A1};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_A2_phi     = fs->make<TH1D>("h_A2_phi"   , "phi^{gen A2};phi^{gen A2};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau1_phi     = fs->make<TH1D>("h_Tau1_phi"   , "phi^{gen Tau1};phi^{gen Tau1};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau2_phi     = fs->make<TH1D>("h_Tau2_phi"   , "phi^{gen Tau2};phi^{gen Tau2};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau3_phi     = fs->make<TH1D>("h_Tau3_phi"   , "phi^{gen Tau3};phi^{gen Tau3};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau4_phi     = fs->make<TH1D>("h_Tau4_phi"   , "phi^{gen Tau4};phi^{gen Tau4};Events"                     ,  20,  -3.2, 3.2);

   H_tau_att_Tau1_Tau2_deta     = fs->make<TH1D>("h_Tau1_Tau2_deta"   , "deta^{gen Tau1_Tau2};deta^{gen Tau1_Tau2};Events"                     ,  20,  0, .8);
   H_tau_att_Tau1_Tau2_dphi     = fs->make<TH1D>("h_Tau1_Tau2_dphi"   , "dphi^{gen Tau1_Tau2};dphi^{gen Tau1_Tau2};Events"                     ,  20,  0, .8);
   H_tau_att_Tau3_Tau4_deta     = fs->make<TH1D>("h_Tau3_Tau4_deta"   , "deta^{gen Tau3_Tau4};deta^{gen Tau3_Tau4};Events"                     ,  20,  0, .8);
   H_tau_att_Tau3_Tau4_dphi     = fs->make<TH1D>("h_Tau3_Tau4_dphi"   , "dphi^{gen Tau3_Tau4};dphi^{gen Tau3_Tau4};Events"                     ,  20,  0, .8);

   H_tau_att_Tau1_Tau2_dphi_deta     = fs->make<TH2D>("h_Tau1_Tau2_dphi_deta"   , "dphi vs deta^{gen Tau1_Tau2};dphi;deta" , 20,  0, .5 , 20,  0, .5);
   H_tau_att_Tau3_Tau4_dphi_deta     = fs->make<TH2D>("h_Tau3_Tau4_dphi_deta"   , "dphi vs deta^{gen Tau3_Tau4};dphi;deta" , 20,  0, .5 , 20,  0, .5);

   RHTree->Branch("Event",  &eventId_);
   RHTree->Branch("Run",  &runId_);
   RHTree->Branch("LumiSection",  &lumiId_);

   RHTree->Branch("GenHiggs_inv",  &V_att_genHiggs_M_inv_);
   RHTree->Branch("GenA1_inv",  &V_att_genA1_M_inv_);
   RHTree->Branch("GenA2_inv",  &V_att_genA2_M_inv_);
   RHTree->Branch("GenHiggs",  &V_att_genHiggs_M_);
   RHTree->Branch("GenA1",  &V_att_genA1_M_);
   RHTree->Branch("GenA2",  &V_att_genA2_M_);
   RHTree->Branch("dR_A1_A2",  &V_att_dR_A1_A2_);
   RHTree->Branch("dR_H_A1",  &V_att_dR_H_A1_);
   RHTree->Branch("dR_H_A2",  &V_att_dR_H_A2_);
   RHTree->Branch("dR_A1_Tau1",  &V_att_dR_A1_Tau1_);
   RHTree->Branch("dR_A1_Tau2",  &V_att_dR_A1_Tau2_);
   RHTree->Branch("dR_A2_Tau3",  &V_att_dR_A2_Tau3_);
   RHTree->Branch("dR_A2_Tau4",  &V_att_dR_A2_Tau4_);
   RHTree->Branch("dR_Tau1_Tau2",  &V_att_dR_Tau1_Tau2_);
   RHTree->Branch("dR_Tau3_Tau4",  &V_att_dR_Tau3_Tau4_);

   RHTree->Branch("H_pt",  &V_att_H_pt_);
   RHTree->Branch("A1_pt",  &V_att_A1_pt_);
   RHTree->Branch("A2_pt",  &V_att_A2_pt_);
   RHTree->Branch("Tau1_pt",  &V_att_Tau1_pt_);
   RHTree->Branch("Tau2_pt",  &V_att_Tau2_pt_);
   RHTree->Branch("Tau3_pt",  &V_att_Tau3_pt_);
   RHTree->Branch("Tau4_pt",  &V_att_Tau4_pt_);
   RHTree->Branch("H_eta",  &V_att_H_eta_);
   RHTree->Branch("A1_eta",  &V_att_A1_eta_);
   RHTree->Branch("A2_eta",  &V_att_A2_eta_);
   RHTree->Branch("Tau1_eta",  &V_att_Tau1_eta_);
   RHTree->Branch("Tau2_eta",  &V_att_Tau2_eta_);
   RHTree->Branch("Tau3_eta",  &V_att_Tau3_eta_);
   RHTree->Branch("Tau4_eta",  &V_att_Tau4_eta_);
   RHTree->Branch("H_phi",  &V_att_H_phi_);
   RHTree->Branch("A1_phi",  &V_att_A1_phi_);
   RHTree->Branch("A2_phi",  &V_att_A2_phi_);
   RHTree->Branch("Tau1_phi",  &V_att_Tau1_phi_);
   RHTree->Branch("Tau2_phi",  &V_att_Tau2_phi_);
   RHTree->Branch("Tau3_phi",  &V_att_Tau3_phi_);
   RHTree->Branch("Tau4_phi",  &V_att_Tau4_phi_);

   RHTree->Branch("Tau1_Tau2_deta",  &V_att_Tau1_Tau2_deta_);
   RHTree->Branch("Tau1_Tau2_dphi",  &V_att_Tau1_Tau2_dphi_);
   RHTree->Branch("Tau3_Tau4_deta",  &V_att_Tau3_Tau4_deta_);
   RHTree->Branch("Tau3_Tau4_dphi",  &V_att_Tau3_Tau4_dphi_);




   genParticlesToken_   = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
   genJetsToken_   = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"));

}


GenAnalyzer::~GenAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
// Define struct to handle mapping for gen pho<->matched reco photons<->matched presel photons
struct jet_tau_map {
  unsigned int idx;
  std::vector<unsigned int> matchedgenJetIdxs;
  std::vector<unsigned int> matchedRecoTauIdxs;
};
std::vector<jet_tau_map> vTauJets;
// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   eventId_ = iEvent.id().event();
   runId_ = iEvent.id().run();
   lumiId_ = iEvent.id().luminosityBlock();

   V_att_genHiggs_M_inv = -1111.1111;
   V_att_genA1_M_inv    = -1111.1111;
   V_att_genA2_M_inv    = -1111.1111;
   V_att_genHiggs_M = -1111.1111;
   V_att_genA1_M    = -1111.1111;
   V_att_genA2_M    = -1111.1111;
   V_att_dR_A1_A2    = -1111.1111;
   V_att_dR_H_A1    = -1111.1111;
   V_att_dR_H_A2    = -1111.1111;
   V_att_dR_A1_Tau1    = -1111.1111;
   V_att_dR_A1_Tau2    = -1111.1111;
   V_att_dR_A2_Tau3    = -1111.1111;
   V_att_dR_A2_Tau4    = -1111.1111;
   V_att_dR_Tau1_Tau2    = -1111.1111;
   V_att_dR_Tau3_Tau4    = -1111.1111;

   V_att_H_pt    = -1111.1111;
   V_att_A1_pt    = -1111.1111;
   V_att_A2_pt    = -1111.1111;
   V_att_Tau1_pt    = -1111.1111;
   V_att_Tau2_pt    = -1111.1111;
   V_att_Tau3_pt    = -1111.1111;
   V_att_Tau4_pt    = -1111.1111;
   V_att_H_eta    = -1111.1111;
   V_att_A1_eta    = -1111.1111;
   V_att_A2_eta    = -1111.1111;
   V_att_Tau1_eta    = -1111.1111;
   V_att_Tau2_eta    = -1111.1111;
   V_att_Tau3_eta    = -1111.1111;
   V_att_Tau4_eta    = -1111.1111;
   V_att_H_phi    = -1111.1111;
   V_att_A1_phi    = -1111.1111;
   V_att_A2_phi    = -1111.1111;
   V_att_Tau1_phi    = -1111.1111;
   V_att_Tau2_phi    = -1111.1111;
   V_att_Tau3_phi    = -1111.1111;
   V_att_Tau4_phi    = -1111.1111;

   V_att_Tau1_Tau2_deta    = -1111.1111;
   V_att_Tau1_Tau2_dphi    = -1111.1111;
   V_att_Tau3_Tau4_deta    = -1111.1111;
   V_att_Tau3_Tau4_dphi    = -1111.1111;

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken(genParticlesToken_,   genParticles);

   edm::Handle<std::vector<reco::GenJet> > genJets;
   iEvent.getByToken(genJetsToken_,   genJets);


  // unsigned int NAs = 0;
  // unsigned int NTau_fromA = 0;
  // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  //   if ( abs(iGen->pdgId()) != 15 || abs(iGen->mother()->pdgId()) != 25) continue;
  //   NTau_fromA++;
  // }
  // if (debug) std::cout << "  >>>>>> Number Tau from  A <<<<<"<<"    "  <<  NTau_fromA << std::endl;
  // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  //   if ( abs(iGen->pdgId()) != 25 || abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15) continue;
  //   NAs++;
  // }
  // if (debug) std::cout << "  >>>>>> Number of A giving Tau <<<<<"<<"    "  <<  NAs << std::endl;

float genHiggs_mass_inv = -1111.1111;
float genA1_mass_inv = -1111.1111;
float genA2_mass_inv = -1111.1111;
float genHiggs_mass = -1111.1111;
float genA1_mass = -1111.1111;
float genA2_mass = -1111.1111;
float A1_A2_dR = -1111.1111;
float H_A1_dR = -1111.1111;
float H_A2_dR = -1111.1111;
float A1_Tau1_dR = -1111.1111;
float A1_Tau2_dR = -1111.1111;
float A2_Tau3_dR = -1111.1111;
float A2_Tau4_dR = -1111.1111;
float Tau1_Tau2_dR = -1111.1111;
float Tau3_Tau4_dR = -1111.1111;

float H_pt = -1111.1111;
float A1_pt = -1111.1111;
float A2_pt = -1111.1111;
float Tau1_pt = -1111.1111;
float Tau2_pt = -1111.1111;
float Tau3_pt = -1111.1111;
float Tau4_pt = -1111.1111;
float H_eta = -1111.1111;
float A1_eta = -1111.1111;
float A2_eta = -1111.1111;
float Tau1_eta = -1111.1111;
float Tau2_eta = -1111.1111;
float Tau3_eta = -1111.1111;
float Tau4_eta = -1111.1111;
float H_phi = -1111.1111;
float A1_phi = -1111.1111;
float A2_phi = -1111.1111;
float Tau1_phi = -1111.1111;
float Tau2_phi = -1111.1111;
float Tau3_phi = -1111.1111;
float Tau4_phi = -1111.1111;

float Tau1_Tau2_deta = -1111.1111;
float Tau1_Tau2_dphi = -1111.1111;
float Tau3_Tau4_deta = -1111.1111;
float Tau3_Tau4_dphi = -1111.1111;


//-----------------------------------------------------------------------------------------------------------------------------------------------
float dR;
std::vector<unsigned int> vGenTauIdxs;
std::vector<unsigned int> vJetIdxs;
bool passedGenSel;
// unsigned int iGenParticle = 0;
passedGenSel = false;
// Find matched gen Tau index
for (unsigned int iG=0; iG < genParticles->size(); iG++) {

  reco::GenParticleRef iGen( genParticles, iG );
  if ((std::abs(iGen->pdgId()) != 15 || iGen->status() != 2) || std::abs(iGen->mother()->pdgId()) != 25 || std::abs(iGen->mother()->mother()->pdgId()) != 35 || iGen->mother()->mother()->numberOfDaughters() != 2) continue;
  float dR_Tau1_Tau2 = reco::deltaR( iGen->mother()->mother()->daughter(0)->daughter(0)->eta(), iGen->mother()->mother()->daughter(0)->daughter(0)->phi(), iGen->mother()->mother()->daughter(0)->daughter(1)->eta(), iGen->mother()->mother()->daughter(0)->daughter(1)->phi());
  float dR_Tau3_Tau4 = reco::deltaR( iGen->mother()->mother()->daughter(1)->daughter(0)->eta(), iGen->mother()->mother()->daughter(1)->daughter(0)->phi(), iGen->mother()->mother()->daughter(1)->daughter(1)->eta(), iGen->mother()->mother()->daughter(1)->daughter(1)->phi());
  if (dR_Tau1_Tau2 < 0.4 || dR_Tau3_Tau4 < 0.4) continue;

  vGenTauIdxs.push_back(iG);
  // if (debug) std::cout << "  >>>>>> Genparticle idx:  "<<iG<< std::endl;
} //genparticles  vGenTauIdxs
if (debug) std::cout << "  >>>>>> Total gen tau matched :  "<<vGenTauIdxs.size()<< std::endl;
if (vGenTauIdxs.size() > 3){
passedGenSel = true;
}
if (debug) std::cout << "  >>>>>> passedGenSel>>> :  "<<passedGenSel<< std::endl;


float minDR = 100.;
int minDR_idx = -1;
vTauJets.clear();

// / Find matched gen Jet index matched to gen Tau
////////// Build gen Tau-jet mapping //////////

// Create mapping between gen tau<->matched jets
// For each gen tau, find "reco" jets matched to it,
// Loop over valid gen Tau idxs
  for ( auto& iG : vGenTauIdxs ) {
    reco::GenParticleRef iGenTau( genParticles, iG );
    std::vector<unsigned int> vMatchedgenJetIdxs;
    // Do dR match to closest reco jets
    minDR = 100.;
    minDR_idx = -1;

    for ( unsigned int iJ = 0; iJ < genJets->size(); iJ++ ) {
       reco::GenJetRef iJet( genJets, iJ );
       dR = reco::deltaR( iJet->eta(),iJet->phi(), iGenTau->eta(),iGenTau->phi() );
       if ( dR > minDR ) continue;
       minDR = dR;
       minDR_idx = iJ;

      } // gen jets
      // Require minimum dR to declare match
      // Protects against matching to PU
      // minDR only needs to be generous enough so that one of the gen taus match to a reco jets for analysis
    if ( minDR > 0.4 ) continue;
    // Declare gen jet matching to gen tau: only store unique reco idxs
    if ( std::find(vMatchedgenJetIdxs.begin(), vMatchedgenJetIdxs.end(), minDR_idx) != vMatchedgenJetIdxs.end() ) continue;
      vMatchedgenJetIdxs.push_back( minDR_idx );

    // store the matched jets ID to a vector
    for ( auto& iJ : vMatchedgenJetIdxs ) {
      vJetIdxs.push_back( iJ );
      // if (debug) std::cout << " ----> matched jet [" << iJ << "] " <<  std::endl;
      }
    // Store this mappin
    jet_tau_map iTau_obj = { iG, vMatchedgenJetIdxs };
    vTauJets.push_back( iTau_obj );
    // if (debug) std::cout << " ----> jet_tau_map " << iTau_obj <<  std::endl;
   } // selected genTua
if (debug) std::cout << " ----> Total matched gen matched jet " << vTauJets.size() <<  std::endl;


// float inv_jet_mass=0;
// float inv_tau_mass=0;
for ( auto const& ii: vTauJets ) {

    // Skip electrons which fails HE edge cut
    if(ii.matchedgenJetIdxs.empty())continue;
    if ( std::find(vJetIdxs.begin(), vJetIdxs.end(), ii.matchedgenJetIdxs[0]) == vJetIdxs.end()) continue;

    reco::GenParticleRef iGen( genParticles, ii.idx );
    reco::GenJetRef iJet( genJets, ii.matchedgenJetIdxs[0] );

    // if (debug) std::cout << " ii.idx: "<< ii.idx  << std::endl;
    // if (debug) std::cout << " ii.matchedgenJetIdxs[0]: "<< ii.matchedgenJetIdxs[0] << std::endl;
    // if (debug) std::cout << " Gen pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;
    // if (debug) std::cout << " Jet pt: "<< iJet->pt() << " eta: " <<iJet->eta() << " phi: " <<iJet->phi() << std::endl;
    if (debug) std::cout << "  >>>>>> Jet  <<<<<"<<"<<< status: "<<iJet->status()<<"<<<pt:  "<<iJet->pt()<<"<<<eta:  "<<iJet->eta()<<"<<<phi:  "<<iJet->phi()<<"<<<mass:  "<<iJet->mass() << std::endl;
    if (debug) std::cout << "  >>>>>> gen Tau  <<<<<"<<"<<< status: "<<iGen->status()<<"<<<pt:  "<<iGen->pt()<<"<<<eta:  "<<iGen->eta()<<"<<<phi:  "<<iGen->phi()<<"<<<mass:  "<<iGen->mass() << std::endl;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------

// Match gen to gen Jet
if (vTauJets.size()==4 && passedGenSel) {
if (debug) std::cout<<">>>>>genJet size"<<genJets->size()<<std::endl;
for ( unsigned int iJ1(0); iJ1 != genJets->size(); ++iJ1 ) {
   reco::GenJetRef iJet1( genJets, iJ1 );
    bool skip_tau = true;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      if (!(iGen->isLastCopy()) || iGen->pdgId() != 35 || iGen->daughter(0)->pdgId() != 25 || abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->numberOfDaughters() != 2 || iGen->daughter(0)->daughter(0)->status() != 2) continue;
      float jettaudR1 = reco::deltaR( iJet1->eta(),iJet1->phi(), iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi() );
      float gendRtautau12 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi());
      if ( jettaudR1 < 0.4 && gendRtautau12 > 0.4 ){

        if (debug) std::cout << "  >>>>>> gen Tau1 <<<<<"<<"<<< status: "<<iGen->daughter(0)->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(0)->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(0)->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(0)->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(0)->mass() << std::endl;
        skip_tau = false;
        break;
      }
    }
    if ( skip_tau ) {
    // if (debug) std::cout << " JET DO NOT MATCH A GEN TAU1" << std::endl;
      continue;
    }

  // unsigned int tau_combinations1 = 0;
  for ( unsigned int iJ2(0); iJ2 != genJets->size(); ++iJ2 ) {
    if ( iJ2 == iJ1 ) continue;
    reco::GenJetRef iJet2( genJets, iJ2 );

      bool skip_tau = true;
      for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {

        if (!(iGen->isLastCopy()) || iGen->pdgId() != 35 || iGen->daughter(0)->pdgId() != 25 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->numberOfDaughters() != 2 || iGen->daughter(0)->daughter(1)->status() != 2) continue;
        float jettaudR2 = reco::deltaR( iJet2->eta(),iJet2->phi(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi() );
        float gendRtautau12 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi() );
        if ( jettaudR2 < 0.4 && gendRtautau12 > 0.4 ){
          skip_tau = false;
          if (debug) std::cout << "  >>>>>> gen Tau2 <<<<<"<<"<<< status: "<<iGen->daughter(0)->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(0)->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(0)->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(0)->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(1)->mass() << std::endl;
          break;
        }
      }
      if ( skip_tau ) {
        // if (debug) std::cout << " JET DO NOT MATCH A GEN TAU2" << std::endl;
        continue;
      }

    float jetjetdR12 = reco::deltaR( iJet1->eta(),iJet1->phi(), iJet2->eta(),iJet2->phi() );
    if ( jetjetdR12 < 0.4 ) continue;
    if (debug) std::cout << " Jet dR = " << jetjetdR12 << std::endl;
    // ++tau_combinations;




  for ( unsigned int iJ3(0); iJ3 != genJets->size(); ++iJ3 ) {
    if ( iJ3 == iJ1 || iJ3==iJ2) continue;
     reco::GenJetRef iJet3( genJets, iJ3 );
      bool skip_tau = true;
      for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
        if (!(iGen->isLastCopy()) || iGen->pdgId() != 35 || iGen->daughter(1)->pdgId() != 25 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15 || iGen->numberOfDaughters() != 2 || iGen->daughter(1)->numberOfDaughters() != 2 || iGen->daughter(1)->daughter(0)->status() != 2) continue;
        float jettaudR3 = reco::deltaR( iJet3->eta(),iJet3->phi(), iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi() );
        float gendRtautau34 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi());
        if ( jettaudR3 < 0.4 && gendRtautau34 > 0.4 ){
          skip_tau = false;
          if (debug) std::cout << "  >>>>>> gen Tau3 <<<<<"<<"<<< status: "<<iGen->daughter(1)->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(1)->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(1)->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(1)->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(1)->daughter(0)->mass() << std::endl;
          break;
        }
      }
      if ( skip_tau ) {
      // if (debug) std::cout << " JET DO NOT MATCH A GEN TAU3" << std::endl;
        continue;
      }

    // unsigned int tau_combinations1 = 0;
    for ( unsigned int iJ4(0); iJ4 != genJets->size(); ++iJ4 ) {
      if ( iJ4 == iJ1 || iJ4==iJ2 || iJ4==iJ3  ) continue;
      reco::GenJetRef iJet4( genJets, iJ4 );
        bool skip_tau = true;
        for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {

          if (!(iGen->isLastCopy()) || iGen->pdgId() != 35 || iGen->daughter(1)->pdgId() != 25 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 15 || iGen->numberOfDaughters() != 2 || iGen->daughter(1)->numberOfDaughters() != 2 || iGen->daughter(1)->daughter(1)->status() != 2) continue;
          float jettaudR4 = reco::deltaR( iJet4->eta(),iJet4->phi(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi() );
          float gendRtautau34 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi() );

          if ( jettaudR4 < 0.4 && gendRtautau34 > 0.4 ){
            skip_tau = false;
            if (debug) std::cout << "  >>>>>> gen Tau4 <<<<<"<<"<<< status: "<<iGen->daughter(1)->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(1)->daughter(1)->mass() << std::endl;
            break;
          }
        }
        if ( skip_tau ) {
          // if (debug) std::cout << " JET DO NOT MATCH A GEN TAU4" << std::endl;
          continue;
        }

      float jetjetdR34 = reco::deltaR( iJet3->eta(),iJet3->phi(), iJet4->eta(),iJet4->phi() );
      if ( jetjetdR34 < 0.4 ) continue;
      if (debug) std::cout << " Jet dR = " << jetjetdR34 << std::endl;


      if (debug) std::cout << "  >>>>>> Jet1  <<<<<"<<"<<< status: "<<iJet1->status()<<"<<<pt:  "<<iJet1->pt()<<"<<<eta:  "<<iJet1->eta()<<"<<<phi:  "<<iJet1->phi()<<"<<<mass:  "<<iJet1->mass() << std::endl;
      if (debug) std::cout << "  >>>>>> Jet2  <<<<<"<<"<<< status: "<<iJet2->status()<<"<<<pt:  "<<iJet2->pt()<<"<<<eta:  "<<iJet2->eta()<<"<<<phi:  "<<iJet2->phi()<<"<<<mass:  "<<iJet2->mass() << std::endl;
      if (debug) std::cout << "  >>>>>> Jet3  <<<<<"<<"<<< status: "<<iJet3->status()<<"<<<pt:  "<<iJet3->pt()<<"<<<eta:  "<<iJet3->eta()<<"<<<phi:  "<<iJet3->phi()<<"<<<mass:  "<<iJet3->mass() << std::endl;
      if (debug) std::cout << "  >>>>>> Jet4  <<<<<"<<"<<< status: "<<iJet4->status()<<"<<<pt:  "<<iJet4->pt()<<"<<<eta:  "<<iJet4->eta()<<"<<<phi:  "<<iJet4->phi()<<"<<<mass:  "<<iJet4->mass() << std::endl;
      TLorentzVector l_gen_jet1  = SetTaus(iJet1->pt(), iJet1->eta(), iJet1->phi(), iJet1->mass());
      TLorentzVector l_gen_jet2  = SetTaus(iJet2->pt(), iJet2->eta(), iJet2->phi(), iJet2->mass());
      TLorentzVector l_gen_jet3  = SetTaus(iJet3->pt(), iJet3->eta(), iJet3->phi(), iJet3->mass());
      TLorentzVector l_gen_jet4  = SetTaus(iJet4->pt(), iJet4->eta(), iJet4->phi(), iJet4->mass());
      TLorentzVector  l_GenA1_jet = l_gen_jet1 + l_gen_jet2;
      TLorentzVector  l_GenA2_jet = l_gen_jet3 + l_gen_jet4;
      if (debug) std::cout << "  >>>>>> Invariant of A1 from gen Jet <<<<<"<<l_GenA1_jet.M()<<std::endl;
      if (debug) std::cout << "  >>>>>> Invariant of A2 from gen Jet <<<<<"<<l_GenA2_jet.M()<<std::endl;




}//iJ4
}//iJ3
} // iJ2

  // if ( tau_combinations1 == 0 ) continue;
} //end iJ1
}
// ------------------------------------------------------------------------------------------------------------------------------------------------------------
bool pass = false;
// if (vTauJets.size()==4) {
for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  if (pass) continue;
  if ( !(iGen->isLastCopy()) || iGen->pdgId() != 35 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 25 || iGen->daughter(1)->pdgId() != 25 ) continue;
  if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 15 ) continue;
  if ( abs(iGen->daughter(0)->daughter(0)->status()) != 2 || abs(iGen->daughter(0)->daughter(1)->status()) != 2 || abs(iGen->daughter(1)->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->daughter(1)->status()) != 2 ) continue;


  TLorentzVector GenTau1  = SetTaus(iGen->daughter(0)->daughter(0)->pt(), iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(0)->mass());
  TLorentzVector GenTau2  = SetTaus(iGen->daughter(0)->daughter(1)->pt(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi(), iGen->daughter(0)->daughter(1)->mass());
  TLorentzVector GenA1 = GenTau1 + GenTau2;
  TLorentzVector GenTau3  = SetTaus(iGen->daughter(1)->daughter(0)->pt(), iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(0)->mass());
  TLorentzVector GenTau4  = SetTaus(iGen->daughter(1)->daughter(1)->pt(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi(), iGen->daughter(1)->daughter(1)->mass());
  TLorentzVector GenA2 = GenTau3 + GenTau4;
  TLorentzVector GenHiggs = GenA1 + GenA2;

  genHiggs_mass_inv = GenHiggs.M();
  genA1_mass_inv = GenA1.M();
  genA2_mass_inv = GenA2.M();
  genHiggs_mass = iGen->mass();
  genA1_mass = iGen->daughter(0)->mass();
  genA2_mass = iGen->daughter(1)->mass();



  float dR_A1_A2 = reco::deltaR( iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
  float dR_H_A1 = reco::deltaR( iGen->eta(), iGen->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
  float dR_H_A2 = reco::deltaR( iGen->eta(), iGen->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
  float dR_A1_Tau1 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
  float dR_A1_Tau2 = reco::deltaR( iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi());
  float dR_A2_Tau3 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
  float dR_A2_Tau4 = reco::deltaR( iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());
  float dR_Tau1_Tau2 = reco::deltaR( iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi());
  float dR_Tau3_Tau4 = reco::deltaR( iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi());

  A1_A2_dR = dR_A1_A2;
  H_A1_dR = dR_H_A1;
  H_A2_dR = dR_H_A2;
  A1_Tau1_dR = dR_A1_Tau1;
  A1_Tau2_dR = dR_A1_Tau2;
  A2_Tau3_dR = dR_A2_Tau3;
  A2_Tau4_dR = dR_A2_Tau4;
  Tau1_Tau2_dR = dR_Tau1_Tau2;
  Tau3_Tau4_dR = dR_Tau3_Tau4;

  H_pt = iGen->pt();
  A1_pt = iGen->daughter(0)->pt();
  A2_pt = iGen->daughter(1)->pt();
  Tau1_pt = iGen->daughter(0)->daughter(0)->pt();
  Tau2_pt = iGen->daughter(0)->daughter(1)->pt();
  Tau3_pt = iGen->daughter(1)->daughter(0)->pt();
  Tau4_pt = iGen->daughter(1)->daughter(1)->pt();
  H_eta = iGen->eta();
  A1_eta = iGen->daughter(0)->eta();
  A2_eta = iGen->daughter(1)->eta();
  Tau1_eta = iGen->daughter(0)->daughter(0)->eta();
  Tau2_eta = iGen->daughter(0)->daughter(1)->eta();
  Tau3_eta = iGen->daughter(1)->daughter(0)->eta();
  Tau4_eta = iGen->daughter(1)->daughter(1)->eta();
  H_phi = iGen->phi();
  A1_phi = iGen->daughter(0)->phi();
  A2_phi = iGen->daughter(1)->phi();
  Tau1_phi = iGen->daughter(0)->daughter(0)->phi();
  Tau2_phi = iGen->daughter(0)->daughter(1)->phi();
  Tau3_phi = iGen->daughter(1)->daughter(0)->phi();
  Tau4_phi = iGen->daughter(1)->daughter(1)->phi();

  Tau1_Tau2_deta = abs(Tau1_eta-Tau2_eta);
  Tau1_Tau2_dphi = abs(Tau1_phi-Tau2_phi);
  Tau3_Tau4_deta = abs(Tau3_eta-Tau4_eta);
  Tau3_Tau4_dphi = abs(Tau3_phi-Tau4_phi);


  // if (debug) std::cout << "  >>>>>> Higgs gen (35) <<<<<"<<"<<< status: "<<iGen->status()<<"<<<pt:  "<<iGen->pt()<<"<<<eta:  "<<iGen->eta()<<"<<<phi:  "<<iGen->phi()<<"<<<mass:  "<<iGen->mass() << std::endl;
  // if (debug) std::cout << "  >>>>>> Higgs LorentzVector (35) <<<<<"<<"<<<mass:  "<<GenHiggs.M() << std::endl;
  // if (debug) std::cout << "  >>>>>> A1 gen (25) <<<<<"<<"<<< status: "<<iGen->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(0)->mass() << std::endl;
  // if (debug) std::cout << "  >>>>>> A1 LorentzVector (25) <<<<<"<<"<<<mass:  "<< GenA1.M() << std::endl;
  // if (debug) std::cout << "  >>>>>> A2 gen (25) <<<<<"<<"<<< status: "<<iGen->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(1)->mass() << std::endl;
  // if (debug) std::cout << "  >>>>>> A2 LorentzVector (25) <<<<<"<<"<<<mass:  "<< GenA2.M() << std::endl;
  // if (debug) std::cout << "  >>>>>> Tau1 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(0)->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(0)->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(0)->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(0)->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(0)->mass() << std::endl;
  // if (debug) std::cout << "  >>>>>> Tau2 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(0)->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(0)->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(0)->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(0)->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(1)->mass() << std::endl;
  // if (debug) std::cout << "  >>>>>> Tau3 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(1)->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(1)->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(1)->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(1)->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(1)->daughter(0)->mass() << std::endl;
  // if (debug) std::cout << "  >>>>>> Tau4 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(1)->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(1)->daughter(1)->mass() << std::endl;
  //
  // if (debug) std::cout << "  >>>>>> dR_A1_A2:  "<<dR_A1_A2<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_H_A1:  "<<dR_H_A1<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_H_A2:  "<<dR_H_A2<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_A1_Tau1:  "<<dR_A1_Tau1<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_A1_Tau2:  "<<dR_A1_Tau2<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_A2_Tau3:  "<<dR_A2_Tau3<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_A2_Tau4:  "<<dR_A2_Tau4<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_Tau1_Tau2:  "<<dR_Tau1_Tau2<< std::endl;
  // if (debug) std::cout << "  >>>>>> dR_Tau3_Tau4:  "<<dR_Tau3_Tau4<< std::endl;
  pass = true;
}
// }
ntotal_event++;
if (pass) {
npassed_event++;

V_att_genHiggs_M_inv    = genHiggs_mass_inv;
V_att_genA1_M_inv    = genA1_mass_inv;
V_att_genA2_M_inv    = genA2_mass_inv;
V_att_genHiggs_M    = genHiggs_mass;
V_att_genA1_M    = genA1_mass;
V_att_genA2_M    = genA2_mass;
V_att_dR_A1_A2    = A1_A2_dR;
V_att_dR_H_A1    = H_A1_dR;
V_att_dR_H_A2    = H_A2_dR;
V_att_dR_A1_Tau1    = A1_Tau1_dR;
V_att_dR_A1_Tau2    = A1_Tau2_dR;
V_att_dR_A2_Tau3    = A2_Tau3_dR;
V_att_dR_A2_Tau4    = A2_Tau4_dR;
V_att_dR_Tau1_Tau2    = Tau1_Tau2_dR;
V_att_dR_Tau3_Tau4    = Tau3_Tau4_dR;

V_att_H_pt    = H_pt;
V_att_A1_pt    = A1_pt;
V_att_A2_pt    = A2_pt;
V_att_Tau1_pt    = Tau1_pt;
V_att_Tau2_pt    = Tau2_pt;
V_att_Tau3_pt    = Tau3_pt;
V_att_Tau4_pt    = Tau4_pt;
V_att_H_eta    = H_eta;
V_att_A1_eta    = A1_eta;
V_att_A2_eta    = A2_eta;
V_att_Tau1_eta    = Tau1_eta;
V_att_Tau2_eta    = Tau2_eta;
V_att_Tau3_eta    = Tau3_eta;
V_att_Tau4_eta    = Tau4_eta;
V_att_H_phi    = H_phi;
V_att_A1_phi    = A1_phi;
V_att_A2_phi    = A2_phi;
V_att_Tau1_phi    = Tau1_phi;
V_att_Tau2_phi    = Tau2_phi;
V_att_Tau3_phi    = Tau3_phi;
V_att_Tau4_phi    = Tau4_phi;

V_att_Tau1_Tau2_deta    = Tau1_Tau2_deta;
V_att_Tau1_Tau2_dphi    = Tau1_Tau2_dphi;
V_att_Tau3_Tau4_deta    = Tau3_Tau4_deta;
V_att_Tau3_Tau4_dphi    = Tau3_Tau4_dphi;

V_att_genHiggs_M_inv_.push_back( V_att_genHiggs_M_inv );
V_att_genA1_M_inv_.push_back( V_att_genA1_M_inv );
V_att_genA2_M_inv_.push_back( V_att_genA2_M_inv );
V_att_genHiggs_M_.push_back( V_att_genHiggs_M );
V_att_genA1_M_.push_back( V_att_genA1_M );
V_att_genA2_M_.push_back( V_att_genA2_M );
V_att_dR_A1_A2_.push_back( V_att_dR_A1_A2 );
V_att_dR_H_A1_.push_back( V_att_dR_H_A1 );
V_att_dR_H_A2_.push_back( V_att_dR_H_A2 );
V_att_dR_A1_Tau1_.push_back( V_att_dR_A1_Tau1 );
V_att_dR_A1_Tau2_.push_back( V_att_dR_A1_Tau2 );
V_att_dR_A2_Tau3_.push_back( V_att_dR_A2_Tau3 );
V_att_dR_A2_Tau4_.push_back( V_att_dR_A2_Tau4 );
V_att_dR_Tau1_Tau2_.push_back( V_att_dR_Tau1_Tau2 );
V_att_dR_Tau3_Tau4_.push_back( V_att_dR_Tau3_Tau4 );

V_att_H_pt_.push_back( V_att_H_pt );
V_att_A1_pt_.push_back( V_att_A1_pt );
V_att_A2_pt_.push_back( V_att_A2_pt );
V_att_Tau1_pt_.push_back( V_att_Tau1_pt );
V_att_Tau2_pt_.push_back( V_att_Tau2_pt );
V_att_Tau3_pt_.push_back( V_att_Tau3_pt );
V_att_Tau4_pt_.push_back( V_att_Tau4_pt );
V_att_H_eta_.push_back( V_att_H_eta );
V_att_A1_eta_.push_back( V_att_A1_eta );
V_att_A2_eta_.push_back( V_att_A2_eta );
V_att_Tau1_eta_.push_back( V_att_Tau1_eta );
V_att_Tau2_eta_.push_back( V_att_Tau2_eta );
V_att_Tau3_eta_.push_back( V_att_Tau3_eta );
V_att_Tau4_eta_.push_back( V_att_Tau4_eta );
V_att_H_phi_.push_back( V_att_H_phi );
V_att_A1_phi_.push_back( V_att_A1_phi );
V_att_A2_phi_.push_back( V_att_A2_phi );
V_att_Tau1_phi_.push_back( V_att_Tau1_phi );
V_att_Tau2_phi_.push_back( V_att_Tau2_phi );
V_att_Tau3_phi_.push_back( V_att_Tau3_phi );
V_att_Tau4_phi_.push_back( V_att_Tau4_phi );

V_att_Tau1_Tau2_deta_.push_back( V_att_Tau1_Tau2_deta );
V_att_Tau1_Tau2_dphi_.push_back( V_att_Tau1_Tau2_dphi );
V_att_Tau3_Tau4_deta_.push_back( V_att_Tau3_Tau4_deta );
V_att_Tau3_Tau4_dphi_.push_back( V_att_Tau3_Tau4_dphi );


H_tau_att_genHiggs_M_inv->Fill( V_att_genHiggs_M_inv );
H_tau_att_genA1_M_inv->Fill( V_att_genA1_M_inv );
H_tau_att_genA2_M_inv->Fill( V_att_genA2_M_inv );
H_tau_att_genHiggs_M->Fill( V_att_genHiggs_M );
H_tau_att_genA1_M->Fill( V_att_genA1_M );
H_tau_att_genA2_M->Fill( V_att_genA2_M );
H_tau_att_dR_A1_A2->Fill( V_att_dR_A1_A2 );
H_tau_att_dR_H_A1->Fill( V_att_dR_H_A1 );
H_tau_att_dR_H_A2->Fill( V_att_dR_H_A2 );
H_tau_att_dR_A1_Tau1->Fill( V_att_dR_A1_Tau1 );
H_tau_att_dR_A1_Tau2->Fill( V_att_dR_A1_Tau2 );
H_tau_att_dR_A2_Tau3->Fill( V_att_dR_A2_Tau3 );
H_tau_att_dR_A2_Tau4->Fill( V_att_dR_A2_Tau4 );
H_tau_att_dR_Tau1_Tau2->Fill( V_att_dR_Tau1_Tau2 );
H_tau_att_dR_Tau3_Tau4->Fill( V_att_dR_Tau3_Tau4 );

H_tau_att_H_pt->Fill( V_att_H_pt );
H_tau_att_A1_pt->Fill( V_att_A1_pt );
H_tau_att_A2_pt->Fill( V_att_A2_pt );
H_tau_att_Tau1_pt->Fill( V_att_Tau1_pt );
H_tau_att_Tau2_pt->Fill( V_att_Tau2_pt );
H_tau_att_Tau3_pt->Fill( V_att_Tau3_pt );
H_tau_att_Tau4_pt->Fill( V_att_Tau4_pt );
H_tau_att_H_eta->Fill( V_att_H_eta );
H_tau_att_A1_eta->Fill( V_att_A1_eta );
H_tau_att_A2_eta->Fill( V_att_A2_eta );
H_tau_att_Tau1_eta->Fill( V_att_Tau1_eta );
H_tau_att_Tau2_eta->Fill( V_att_Tau2_eta );
H_tau_att_Tau3_eta->Fill( V_att_Tau3_eta );
H_tau_att_Tau4_eta->Fill( V_att_Tau4_eta );
H_tau_att_H_phi->Fill( V_att_H_phi );
H_tau_att_A1_phi->Fill( V_att_A1_phi );
H_tau_att_A2_phi->Fill( V_att_A2_phi );
H_tau_att_Tau1_phi->Fill( V_att_Tau1_phi );
H_tau_att_Tau2_phi->Fill( V_att_Tau2_phi );
H_tau_att_Tau3_phi->Fill( V_att_Tau3_phi );
H_tau_att_Tau4_phi->Fill( V_att_Tau4_phi );

H_tau_att_Tau1_Tau2_deta->Fill( V_att_Tau1_Tau2_deta );
H_tau_att_Tau1_Tau2_dphi->Fill( V_att_Tau1_Tau2_dphi );
H_tau_att_Tau3_Tau4_deta->Fill( V_att_Tau3_Tau4_deta );
H_tau_att_Tau3_Tau4_dphi->Fill( V_att_Tau3_Tau4_dphi );

H_tau_att_Tau1_Tau2_dphi_deta->Fill(V_att_Tau1_Tau2_dphi,V_att_Tau1_Tau2_deta);
H_tau_att_Tau3_Tau4_dphi_deta->Fill(V_att_Tau3_Tau4_dphi,V_att_Tau3_Tau4_deta);

RHTree->Fill();

}

V_att_genHiggs_M_inv_.clear();
V_att_genA1_M_inv_.clear();
V_att_genA2_M_inv_.clear();
V_att_genHiggs_M_.clear();
V_att_genA1_M_.clear();
V_att_genA2_M_.clear();
V_att_dR_A1_A2_.clear();
V_att_dR_H_A1_.clear();
V_att_dR_H_A2_.clear();
V_att_dR_A1_Tau1_.clear();
V_att_dR_A1_Tau2_.clear();
V_att_dR_A2_Tau3_.clear();
V_att_dR_A2_Tau4_.clear();
V_att_dR_Tau1_Tau2_.clear();
V_att_dR_Tau3_Tau4_.clear();

V_att_H_pt_.clear();
V_att_A1_pt_.clear();
V_att_A2_pt_.clear();
V_att_Tau1_pt_.clear();
V_att_Tau2_pt_.clear();
V_att_Tau3_pt_.clear();
V_att_Tau4_pt_.clear();
V_att_H_eta_.clear();
V_att_A1_eta_.clear();
V_att_A2_eta_.clear();
V_att_Tau1_eta_.clear();
V_att_Tau2_eta_.clear();
V_att_Tau3_eta_.clear();
V_att_Tau4_eta_.clear();
V_att_H_phi_.clear();
V_att_A1_phi_.clear();
V_att_A2_phi_.clear();
V_att_Tau1_phi_.clear();
V_att_Tau2_phi_.clear();
V_att_Tau3_phi_.clear();
V_att_Tau4_phi_.clear();

V_att_Tau1_Tau2_deta_.clear();
V_att_Tau1_Tau2_dphi_.clear();
V_att_Tau3_Tau4_deta_.clear();
V_att_Tau3_Tau4_dphi_.clear();


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
   if (debug) std::cout << "  >>>>>> Total events selected events <<<<<  "<<npassed_event<<"/"<<ntotal_event<<std::endl;
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
