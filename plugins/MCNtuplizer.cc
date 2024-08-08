#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Span.h"

// Gen event info (pt hat)
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// Fast jet for additional clustering
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

// root include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TVector2.h"

// l1trigger
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

// Reco collections
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// Pt regression
#include "NtuplizerGenPart/L1ScoutingAnalyzer/interface/PtRegression.h"

#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include <map>
#include <string>

class MCNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MCNtuplizer(const edm::ParameterSet&);
  ~MCNtuplizer() {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;
  // Fill with regressed pT
  void getRegressedPt(PtRegression* ptRegression, std::vector<float> &Jet_pt_regressed, std::vector<float> Jet_pt, std::vector<float> Jet_eta, std::vector<float> Jet_phi, float totalEt, float missingEt, float missingEtPhi, float totalHt, float missingHt, float missingHtPhi, int towerCount, std::vector<float> EGamma_pt, std::vector<float> EGamma_eta, std::vector<float> EGamma_phi, std::vector<float> Muon_pt, std::vector<float> Muon_eta, std::vector<float> Muon_phi);

  // Tokens for gen objects
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genFatJetsToken_;
  // Tokens for l1t objects
  edm::EDGetTokenT<l1t::JetBxCollection> jetToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> egToken_;
  edm::EDGetTokenT<l1t::EtSumBxCollection> etSumToken_;
  edm::EDGetTokenT<l1t::TauBxCollection> tauToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection> muonToken_;
  // Tokens for reco objects
  edm::EDGetTokenT<pat::JetCollection> recoJetToken_;
  edm::EDGetTokenT<pat::JetCollection> recoJetPuppiToken_;
  edm::EDGetTokenT<pat::METCollection> recoMetToken_;
  edm::EDGetTokenT<pat::METCollection> recoMetPuppiToken_;

  // Tree that contains info per bunch crossing
  TTree* tree;

  //Gen info
  Float_t genPtHat;
  Int_t nGenPart;
  Int_t nGenJet;
  Int_t nGenFatJet;
  vector<Float16_t> GenPart_pt;
  vector<Float16_t> GenPart_eta;
  vector<Float16_t> GenPart_phi;
  vector<Float16_t> GenPart_e;
  vector<Float16_t> GenPart_mass;
  vector<Int_t> GenPart_pdgId;
  vector<Int_t> GenPart_charge;
  vector<Int_t> GenPart_genPartIdxMother;
  vector<Int_t> GenPart_MotherpdgId;
  
  vector<Float16_t> GenJet_pt;
  vector<Float16_t> GenJet_eta;
  vector<Float16_t> GenJet_phi;
  vector<Float16_t> GenJet_e;

  //AK8 gen jets
  vector<Float16_t> GenFatJet_pt;
  vector<Float16_t> GenFatJet_eta;
  vector<Float16_t> GenFatJet_phi;
  vector<Float16_t> GenFatJet_e;

  //Reco info
  //Jets
  Int_t nJet;
  vector<Float16_t> Jet_pt;
  // Jet et is the regressed value
  vector<Float16_t> Jet_pt_regressed;
  vector<Float16_t> Jet_eta;
  vector<Float16_t> Jet_phi;
  vector<Float16_t> Jet_e;
  vector<Int_t> Jet_qual;
  vector<Int_t> Jet_partonFlav;
  //vector<Int_t[4]> Jet_puDonutEt;

  //AK8 jets made from L1 jets
  Int_t nFatJet;
  vector<Float16_t> FatJet_pt;
  vector<Float16_t> FatJet_eta;
  vector<Float16_t> FatJet_phi;
  vector<Float16_t> FatJet_e;
  vector<Float16_t> FatJet_m;
  vector<Int_t> FatJet_partonFlav;

  //EGamma
  Int_t nEGamma;
  vector<Float16_t> EGamma_pt;
  vector<Float16_t> EGamma_eta;
  vector<Float16_t> EGamma_phi;
  vector<Float16_t> EGamma_e;
  vector<Int_t> EGamma_Iso;

  //Muons
  Int_t nMuon;
  vector<Float16_t> Muon_pt;
  vector<Float16_t> Muon_eta;
  vector<Float16_t> Muon_phi;
  vector<Float16_t> Muon_e;
  vector<Int_t> Muon_qual;
  vector<Int_t> Muon_hwCharge;
  vector<Float16_t> Muon_etaAtVtx;
  vector<Float16_t> Muon_phiAtVtx;
  vector<Int_t> Muon_hwDXY;

  //et Sums
  float totalEt, totalHt, missingEt, missingEtPhi, missingHt, missingHtPhi;
  int towerCount;

  //PT regression 
  std::string model_path;
  PtRegression* ptRegression;

  //Reco jets
  Int_t nRecoJet;
  vector<Float16_t> RecoJet_pt;
  vector<Float16_t> RecoJet_eta;
  vector<Float16_t> RecoJet_phi;
  vector<Float16_t> RecoJet_e;
  vector<Float16_t> RecoJet_mass;

  //Puppi jets
  Int_t nRecoJetPuppi;
  vector<Float16_t> RecoJetPuppi_pt;
  vector<Float16_t> RecoJetPuppi_eta;
  vector<Float16_t> RecoJetPuppi_phi;
  vector<Float16_t> RecoJetPuppi_e;
  vector<Float16_t> RecoJetPuppi_mass;

  //Reco MET
  Float16_t RecoMet_pt;
  Float16_t RecoMet_phi;

  //Puppi MET
  Float16_t RecoMetPuppi_pt;
  Float16_t RecoMetPuppi_phi;

};

MCNtuplizer::MCNtuplizer(const edm::ParameterSet& iPset)
  :genEventInfoToken_(consumes<GenEventInfoProduct>(iPset.getParameter<edm::InputTag>("genEventInfoTag"))), 
  genParticlesToken_(consumes<reco::GenParticleCollection>(iPset.getParameter<edm::InputTag>("genParticlesTag"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iPset.getParameter<edm::InputTag>("genJetsTag"))),
  genFatJetsToken_(consumes<reco::GenJetCollection>(iPset.getParameter<edm::InputTag>("genFatJetsTag"))),
  jetToken_(consumes<l1t::JetBxCollection>(iPset.getParameter<edm::InputTag>("jetsTag"))),
  egToken_(consumes<l1t::EGammaBxCollection>(iPset.getParameter<edm::InputTag>("eGammasTag"))),
  etSumToken_(consumes<l1t::EtSumBxCollection>(iPset.getParameter<edm::InputTag>("etSumsTag"))),
  tauToken_(consumes<l1t::TauBxCollection>(iPset.getParameter<edm::InputTag>("tausTag"))),
  muonToken_(consumes<l1t::MuonBxCollection>(iPset.getParameter<edm::InputTag>("muonsTag"))),
  recoJetToken_(consumes<pat::JetCollection>(iPset.getParameter<edm::InputTag>("recoJetsTag"))),
  recoJetPuppiToken_(consumes<pat::JetCollection>(iPset.getParameter<edm::InputTag>("recoJetsPuppiTag"))),
  recoMetToken_(consumes<pat::METCollection>(iPset.getParameter<edm::InputTag>("recoMetTag"))),
  recoMetPuppiToken_(consumes<pat::METCollection>(iPset.getParameter<edm::InputTag>("recoMetPuppiTag")))
{

  // the root file service to handle the output file
  edm::Service<TFileService> fs;

  //Initialize pt regression
  //model_path = std::string("/eos/user/p/ppradeep/L1Scouting/May24v2/CMSSW_14_1_0_pre1/src/NtuplizerGenPart/L1ScoutingAnalyzer/models/jet_regression_model_linear.model");
  model_path = std::string("/eos/user/p/ppradeep/L1Scouting/May24v2/CMSSW_14_1_0_pre1/src/NtuplizerGenPart/L1ScoutingAnalyzer/models/jet_regression_model_new.model");
  ptRegression = new PtRegression(model_path);

  // Create the TTree
  tree = fs->make<TTree>("Events", "Events_bx");

  //Gen info
  tree->Branch("genPtHat", &genPtHat, "genPtHat/F");
  tree->Branch("nGenPart", &nGenPart, "nGenPart/I");
  tree->Branch("GenPart_pt", &GenPart_pt);
  tree->Branch("GenPart_eta", &GenPart_eta);
  tree->Branch("GenPart_phi", &GenPart_phi);
  tree->Branch("GenPart_e", &GenPart_e);
  tree->Branch("GenPart_mass", &GenPart_mass);
  tree->Branch("GenPart_pdgId", &GenPart_pdgId);
  tree->Branch("GenPart_charge", &GenPart_charge);
  tree->Branch("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
  tree->Branch("GenPart_MotherpdgId", &GenPart_MotherpdgId);

  tree->Branch("nGenJet", &nGenJet, "nGenJet/I");
  tree->Branch("GenJet_pt", &GenJet_pt);
  tree->Branch("GenJet_eta", &GenJet_eta);
  tree->Branch("GenJet_phi", &GenJet_phi);
  tree->Branch("GenJet_e", &GenJet_e);

  tree->Branch("nGenFatJet", &nGenFatJet, "nGenFatJet/I");
  tree->Branch("GenFatJet_pt", &GenFatJet_pt);
  tree->Branch("GenFatJet_eta", &GenFatJet_eta);
  tree->Branch("GenFatJet_phi", &GenFatJet_phi);
  tree->Branch("GenFatJet_e", &GenFatJet_e);
  
  // Jets
  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("Jet_pt", &Jet_pt);
  tree->Branch("Jet_pt_regressed", &Jet_pt_regressed);
  tree->Branch("Jet_eta", &Jet_eta);
  tree->Branch("Jet_phi", &Jet_phi);
  tree->Branch("Jet_e", &Jet_e);
  tree->Branch("Jet_qual", &Jet_qual);
  tree->Branch("Jet_partonFlav", &Jet_partonFlav);

  // AK8 jets
  tree->Branch("nFatJet", &nFatJet, "nFatJet/I");
  tree->Branch("FatJet_pt", &FatJet_pt);
  tree->Branch("FatJet_eta", &FatJet_eta);
  tree->Branch("FatJet_phi", &FatJet_phi);
  tree->Branch("FatJet_e", &FatJet_e);
  tree->Branch("FatJet_m", &FatJet_m);
  tree->Branch("FatJet_partonFlav", &FatJet_partonFlav);

  // EGammas
  tree->Branch("nEGamma", &nEGamma, "nEGamma/I");
  tree->Branch("EGamma_pt", &EGamma_pt);
  tree->Branch("EGamma_eta", &EGamma_eta);
  tree->Branch("EGamma_phi", &EGamma_phi);
  tree->Branch("EGamma_e", &EGamma_e);
  tree->Branch("EGamma_Iso", &EGamma_Iso);
  
  // Muons
  tree->Branch("nMuon", &nMuon, "nMuon/I");
  tree->Branch("Muon_pt", &Muon_pt);
  tree->Branch("Muon_eta", &Muon_eta);
  tree->Branch("Muon_phi", &Muon_phi);
  tree->Branch("Muon_e", &Muon_e);
  tree->Branch("Muon_qual", &Muon_qual);
  tree->Branch("Muon_hwCharge", &Muon_hwCharge);
  tree->Branch("Muon_etaAtVtx", &Muon_etaAtVtx);
  tree->Branch("Muon_phiAtVtx", &Muon_phiAtVtx);
  tree->Branch("Muon_hwDXY", &Muon_hwDXY);

  // Sums
  tree->Branch("etSum", &totalEt);
  tree->Branch("htSum", &totalHt);
  tree->Branch("etMiss", &missingEt);
  tree->Branch("etMissPhi", &missingEtPhi);
  tree->Branch("htMiss", &missingHt);
  tree->Branch("htMissPhi", &missingHtPhi);
  tree->Branch("towerCount", &towerCount);

  // Reco jet
  tree->Branch("nRecoJet", &nRecoJet, "nRecoJet/I");
  tree->Branch("RecoJet_pt", &RecoJet_pt);
  tree->Branch("RecoJet_eta", &RecoJet_eta);
  tree->Branch("RecoJet_phi", &RecoJet_phi);
  tree->Branch("RecoJet_e", &RecoJet_e);
  tree->Branch("RecoJet_mass", &RecoJet_mass);

  // Reco jet puppi
  tree->Branch("nRecoJetPuppi", &nRecoJetPuppi, "nRecoJetPuppi/I");
  tree->Branch("RecoJetPuppi_pt", &RecoJetPuppi_pt);
  tree->Branch("RecoJetPuppi_eta", &RecoJetPuppi_eta);
  tree->Branch("RecoJetPuppi_phi", &RecoJetPuppi_phi);
  tree->Branch("RecoJetPuppi_e", &RecoJetPuppi_e);
  tree->Branch("RecoJetPuppi_mass", &RecoJetPuppi_mass);

  // Reco MET
  tree->Branch("RecoMet_pt", &RecoMet_pt, "RecoMet_pt/F");
  tree->Branch("RecoMet_phi", &RecoMet_phi, "RecoMet_phi/F");
  
  // Reco MET Puppi
  tree->Branch("RecoMetPuppi_pt", &RecoMetPuppi_pt, "RecoMetPuppi_pt/F");
  tree->Branch("RecoMetPuppi_phi", &RecoMetPuppi_phi, "RecoMetPuppi_phi/F");

}

void MCNtuplizer::getRegressedPt(PtRegression* ptRegression, std::vector<float> &Jet_pt_regressed, std::vector<float> Jet_pt, std::vector<float> Jet_eta, std::vector<float> Jet_phi, float totalEt, float missingEt, float missingEtPhi, float totalHt, float missingHt, float missingHtPhi, int towerCount, std::vector<float> EGamma_pt, std::vector<float> EGamma_eta, std::vector<float> EGamma_phi, std::vector<float> Muon_pt, std::vector<float> Muon_eta, std::vector<float> Muon_phi){
  // Apply regression to L1 jets (already loaded in class initialisation)
  // This gets called in the analyze function
  for (UInt_t i=0; i<Jet_pt.size(); i++){
    float regressed_pt;
    float regressed_pt_scale;

    //Define inputs for isolation
    TLorentzVector EGamma_sum;
    TLorentzVector Muon_sum;
    EGamma_sum.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    Muon_sum.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    float Jet_relEGammaIso = 0;
    float Jet_relMuIso = 0;
    
    for(UInt_t j=0; j<EGamma_pt.size(); j++){
      TLorentzVector EGamma_vector;
      EGamma_vector.SetPtEtaPhiM(EGamma_pt[j], EGamma_eta[j], EGamma_phi[j], 0.0);
      float dphi = TVector2::Phi_mpi_pi(Jet_phi[i] - EGamma_phi[j]);
      float deta = Jet_eta[i] - EGamma_eta[j];
      float dR = std::sqrt(dphi*dphi + deta*deta);
      if(dR < 0.4){
        EGamma_sum += EGamma_vector;
      }
    }
    for(UInt_t j=0; j<Muon_pt.size(); j++){
      TLorentzVector Muon_vector;
      Muon_vector.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.1055);
      float dphi = TVector2::Phi_mpi_pi(Jet_phi[i] - Muon_phi[j]);
      float deta = Jet_eta[i] - Muon_eta[j];
      float dR = std::sqrt(dphi*dphi + deta*deta);
      if(dR < 0.4){
        Muon_sum += Muon_vector;
      }
    }

    float Jet_pt_log = TMath::Log(Jet_pt[i]);

    //Get relative isolation
    Jet_relEGammaIso = EGamma_sum.Pt() / Jet_pt[i];
    Jet_relMuIso = Muon_sum.Pt() / Jet_pt[i];

    //std::cout << "Jet pT: " << Jet_pt[i] << ", Jet pT log: " << Jet_pt_log << ", Jet relEGammaIso: " << Jet_relEGammaIso << ", Jet relMuIso: " << Jet_relMuIso << "\n";

    //Define inputs for regression
    std::vector<float> inputs{Jet_pt_log, Jet_eta[i], Jet_phi[i], Jet_relMuIso, Jet_relEGammaIso, totalEt, missingEt, missingEtPhi, totalHt, missingHt, missingHtPhi, float(towerCount)};

    //Skip jets with pT > 1000 GeV
    if(Jet_pt[i] > 1000.){
      regressed_pt = Jet_pt[i];
    }
    else{
      regressed_pt_scale = ptRegression->get_regressed_pt(inputs);
      regressed_pt = Jet_pt[i] * std::exp(regressed_pt_scale);
      //std::cout << "Jet " << i << " pT: " << Jet_pt[i] << " Regressed target scale: " << regressed_pt_scale << " Regressed et: " << regressed_pt << "\n";
    }
    Jet_pt_regressed.push_back(regressed_pt);
  }
}


void MCNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup&){
  edm::Handle<GenEventInfoProduct> genEventInfo;
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<reco::GenJetCollection> genFatJets;
  edm::Handle<l1t::JetBxCollection> jetsCollection;
  edm::Handle<l1t::EGammaBxCollection> egammasCollection;
  edm::Handle<l1t::EtSumBxCollection> etSumsCollection;
  edm::Handle<l1t::TauBxCollection> tausCollection;
  edm::Handle<l1t::MuonBxCollection> muonsCollection;
  edm::Handle<pat::JetCollection> recoJets;
  edm::Handle<pat::JetCollection> recoJetsPuppi;
  edm::Handle<pat::METCollection> recoMet;
  edm::Handle<pat::METCollection> recoMetPuppi;

  iEvent.getByToken(genEventInfoToken_, genEventInfo);
  iEvent.getByToken(genParticlesToken_, genParticles);
  iEvent.getByToken(genJetsToken_, genJets);
  iEvent.getByToken(genFatJetsToken_, genFatJets);
  iEvent.getByToken(jetToken_, jetsCollection);
  iEvent.getByToken(egToken_, egammasCollection);
  iEvent.getByToken(etSumToken_, etSumsCollection);
  iEvent.getByToken(tauToken_, tausCollection);
  iEvent.getByToken(muonToken_, muonsCollection);
  iEvent.getByToken(recoJetToken_, recoJets);
  iEvent.getByToken(recoJetPuppiToken_, recoJetsPuppi);
  iEvent.getByToken(recoMetToken_, recoMet);
  iEvent.getByToken(recoMetPuppiToken_, recoMetPuppi);

  GenPart_pt.clear();
  GenPart_eta.clear();
  GenPart_phi.clear();
  GenPart_e.clear();
  GenPart_mass.clear();
  GenPart_pdgId.clear();
  GenPart_charge.clear();
  GenPart_genPartIdxMother.clear();
  GenPart_MotherpdgId.clear();

  GenJet_pt.clear();
  GenJet_eta.clear();
  GenJet_phi.clear();
  GenJet_e.clear();

  GenFatJet_pt.clear();
  GenFatJet_eta.clear();
  GenFatJet_phi.clear();
  GenFatJet_e.clear();

  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_e.clear();
  Muon_qual.clear();
  Muon_hwCharge.clear();
  Muon_etaAtVtx.clear();
  Muon_phiAtVtx.clear();
  Muon_hwDXY.clear();

  Jet_pt.clear();
  Jet_pt_regressed.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_e.clear();
  Jet_qual.clear();
  Jet_partonFlav.clear();

  FatJet_pt.clear();
  FatJet_eta.clear();
  FatJet_phi.clear();
  FatJet_e.clear();
  FatJet_m.clear();
  FatJet_partonFlav.clear();
  
  EGamma_pt.clear();
  EGamma_eta.clear();
  EGamma_phi.clear();
  EGamma_e.clear();
  EGamma_Iso.clear();

  RecoJet_pt.clear();
  RecoJet_eta.clear();
  RecoJet_phi.clear();
  RecoJet_e.clear();
  RecoJet_mass.clear();

  RecoJetPuppi_pt.clear();
  RecoJetPuppi_eta.clear();
  RecoJetPuppi_phi.clear();
  RecoJetPuppi_e.clear();
  RecoJetPuppi_mass.clear();

  RecoMet_pt = -999;
  RecoMet_phi = -999;

  RecoMetPuppi_pt = -999;
  RecoMetPuppi_phi = -999;

  //Cluster AK8 jets using jets in fastjet
  vector<fastjet::PseudoJet> FatJetInputs;
  vector<fastjet::PseudoJet> GhostJetInputs;
  double FatJetR = 1.0;
  double JetR = 0.4;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, FatJetR);
  fastjet::JetDefinition jet_def_ghost(fastjet::antikt_algorithm, JetR);

  //Gen event pT hat
  genPtHat = (genEventInfo->hasBinningValues() ? genEventInfo->binningValues()[0] : 0.0);

  //Gen particles
  nGenPart = 0;
  if(genParticles.isValid()){
    for (reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
      //Skip particles that don't pass pT cut (>5 unless muons for which >3 is set)
      //bool fill = false;
      //if(abs(genParticle->pdgId()) == 13){
      //  if(genParticle->pt() > 3) fill = true;
      //}
      //else{
      //  if(genParticle->pt() > 5) fill = true;
      //}
      //if(!fill) continue;
      GenPart_pt.push_back(genParticle->pt());
      GenPart_eta.push_back(genParticle->eta());
      GenPart_phi.push_back(genParticle->phi());
      GenPart_e.push_back(genParticle->energy());
      GenPart_mass.push_back(genParticle->mass());
      GenPart_pdgId.push_back(genParticle->pdgId());
      GenPart_charge.push_back(genParticle->charge());
      //genFatJetInputs.push_back(fastjet::PseudoJet(genParticle->px(), genParticle->py(), genParticle->pz(), genParticle->energy()));
      //Find mother id
      int genPartIdxMother = (genParticle->numberOfMothers()>0) ? (genParticle->motherRef(0)).key() : -1;
      int MotherpdgId = (genPartIdxMother != -1) ? (GenPart_pdgId[genPartIdxMother]) : 0;
      GenPart_genPartIdxMother.push_back(genPartIdxMother);
      GenPart_MotherpdgId.push_back(MotherpdgId);
      nGenPart++;
    }
  }


  //Gen jets
  nGenJet = 0;
  if(genJets.isValid()){
    for (reco::GenJetCollection::const_iterator genJet = genJets->begin(); genJet != genJets->end(); genJet++){
      //Only store gen jets with pT > 10 GeV
      //if(genJet->pt() < 10) continue;
      GenJet_pt.push_back(genJet->pt());
      GenJet_eta.push_back(genJet->eta());
      GenJet_phi.push_back(genJet->phi());
      GenJet_e.push_back(genJet->energy());
      nGenJet++;
    }
  }
  
  //Gen fat jets
  nGenFatJet = 0;
  if(genFatJets.isValid()){
    for (reco::GenJetCollection::const_iterator genFatJet = genFatJets->begin(); genFatJet != genFatJets->end(); genFatJet++){
      //Only store gen jets with pT > 10 GeV
      //if(genJet->pt() < 10) continue;
      GenFatJet_pt.push_back(genFatJet->pt());
      GenFatJet_eta.push_back(genFatJet->eta());
      GenFatJet_phi.push_back(genFatJet->phi());
      GenFatJet_e.push_back(genFatJet->energy());
      nGenFatJet++;
    }
  }

  //Muons
  nMuon = 0;

  for (l1t::MuonBxCollection::const_iterator muon = muonsCollection->begin(0); muon != muonsCollection->end(0); muon++){
    Muon_pt.push_back(muon->pt());
    Muon_eta.push_back(muon->eta());
    Muon_phi.push_back(muon->phi());
    Muon_e.push_back(muon->energy());
    Muon_qual.push_back(muon->hwQual());
    Muon_hwCharge.push_back(muon->hwCharge());
    Muon_etaAtVtx.push_back(muon->etaAtVtx());
    Muon_phiAtVtx.push_back(muon->phiAtVtx());
    Muon_hwDXY.push_back(muon->hwDXY());
    nMuon++;
  }

  nJet = 0;

  for (l1t::JetBxCollection::const_iterator jet = jetsCollection->begin(0); jet != jetsCollection->end(0); jet++){
    Jet_pt.push_back(jet->pt());
    Jet_eta.push_back(jet->eta());
    Jet_phi.push_back(jet->phi());
    Jet_e.push_back(jet->energy());
    Jet_qual.push_back(jet->hwQual());
    //Feed jets into fastjet for clustering
    TLorentzVector jet_vector;
    jet_vector.SetPtEtaPhiM(jet->pt(), jet->eta(), jet->phi(), 0.0);

    //Simple dR matching with gen partons
    int partonFlav = 0;
    std::vector<int> partonFlavs{};
    std::vector<float> partondRs{};
    for(int i = 0; i < nGenPart; i++){
      //Ignore if not a quark
      if((std::abs(GenPart_pdgId[i]) > 5) && (GenPart_pdgId[i]!=21)) continue;
      float dphi = TVector2::Phi_mpi_pi(jet->phi() - GenPart_phi[i]);
      float deta = jet->eta() - GenPart_eta[i];
      float dR = std::sqrt(dphi*dphi + deta*deta);
      if(dR < 0.4){
        partonFlavs.push_back(GenPart_pdgId[i]);
        partondRs.push_back(dR);
      }
    }
    //Assign flavour based on hierarchy 
    if(partonFlavs.size() > 0){
      if(std::find(partonFlavs.begin(), partonFlavs.end(), 5) != partonFlavs.end()){
        partonFlav = 5;
      }
      else if(std::find(partonFlavs.begin(), partonFlavs.end(), -5) != partonFlavs.end()){
        partonFlav = -5;
      }
      else if(std::find(partonFlavs.begin(), partonFlavs.end(), 4) != partonFlavs.end()){
        partonFlav = 4;
      }
      else if(std::find(partonFlavs.begin(), partonFlavs.end(), -4) != partonFlavs.end()){
        partonFlav = -4;
      }
      else{
        //Pick the smallest dR one
        int minElementIndex = std::min_element(partondRs.begin(), partondRs.end()) - partondRs.begin();
        partonFlav = partonFlavs[minElementIndex];
      }
      //std::cout << "Jet " << nJet << " matched to parton with flavour " << partonFlav << "\n";
    }
    Jet_partonFlav.push_back(partonFlav);

    //Get fat jets
    FatJetInputs.push_back(fastjet::PseudoJet(jet_vector.Px(), jet_vector.Py(), jet_vector.Pz(), jet_vector.E()));
    //Ghost clustering (constituents)
    fastjet::PseudoJet ghostJet(jet_vector.Px(), jet_vector.Py(), jet_vector.Pz(), jet_vector.E());
    ghostJet.set_user_index(nJet);
    GhostJetInputs.push_back(ghostJet);
    nJet++;
  }

  
  //Add gen particles to ghost jet inputs
  if(genParticles.isValid()){
    int nGenGhost = 0;
    for (reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
      TLorentzVector gen_vector;
      gen_vector.SetPtEtaPhiM(genParticle->pt(), genParticle->eta(), genParticle->phi(), genParticle->mass());
      fastjet::PseudoJet ghostJet(gen_vector.Px(), gen_vector.Py(), gen_vector.Pz(), gen_vector.E());
      //ghostJet.set_user_index(nJet+nGenGhost);
      ghostJet *= 1e-18;
      GhostJetInputs.push_back(ghostJet);
      //if (nGenGhost > 5) break;
      nGenGhost++;
    }
  }

  /*
  std::cout << "Number of reco jet inputs: " << nJet << "\n";
  std::cout << "Number of reclustered jet inputs: " << GhostJetInputs.size() << "\n";

  std::cout << "Number of reclustered jet inputs: " << GhostJetInputs.size() << "\n";

  //Recluster, match to reco and if matched, check for parton flavour
  fastjet::ClusterSequence cs_ghost(GhostJetInputs, jet_def_ghost);
  std::cout << "Clustering done" << "\n";
  vector<fastjet::PseudoJet> GhostJets = fastjet::sorted_by_pt(cs_ghost.inclusive_jets());

  //Check size
  std::cout << "Jet clustering done: " << "\n";
  std::cout << "Starting jets: " << nJet << "\n";
  std::cout << "Reclustered jets: " << GhostJets.size() << "\n";
  */


  //Cluster and store fat jets
  fastjet::ClusterSequence cs(FatJetInputs, jet_def);
  vector<fastjet::PseudoJet> FatJets = fastjet::sorted_by_pt(cs.inclusive_jets());

  nFatJet = 0;
  for (UInt_t i=0; i < FatJets.size(); i++){
    fastjet::PseudoJet fatjet = FatJets[i];
    auto constituents = fatjet.constituents();
    //Skip fat jets with only one constituent
    if(constituents.size() < 2) continue;
    FatJet_pt.push_back(FatJets[i].pt());
    FatJet_eta.push_back(FatJets[i].eta());
    FatJet_phi.push_back(FatJets[i].phi());
    FatJet_e.push_back(FatJets[i].E());
    FatJet_m.push_back(FatJets[i].m());
    //Get parton flavour of fat jet (dR of 1.0)
    int partonFlav = 0;
    std::vector<int> partonFlavs{};
    std::vector<float> partondRs{};
    for(int j = 0; j < nGenPart; j++){
      //Ignore if not a quark
      if((std::abs(GenPart_pdgId[j]) > 5) && (GenPart_pdgId[j]!=21)) continue;
      float dphi = TVector2::Phi_mpi_pi(FatJets[i].phi() - GenPart_phi[j]);
      float deta = FatJets[i].eta() - GenPart_eta[j];
      float dR = std::sqrt(dphi*dphi + deta*deta);
      if(dR < 1.0){
        partonFlavs.push_back(GenPart_pdgId[j]);
        partondRs.push_back(dR);
      }
    }
    //Assign flavour based on hierarchy
    if(partonFlavs.size() > 0){
      if(std::find(partonFlavs.begin(), partonFlavs.end(), 5) != partonFlavs.end()){
        partonFlav = 5;
      }
      else if(std::find(partonFlavs.begin(), partonFlavs.end(), -5) != partonFlavs.end()){
        partonFlav = -5;
      }
      else if(std::find(partonFlavs.begin(), partonFlavs.end(), 4) != partonFlavs.end()){
        partonFlav = 4;
      }
      else if(std::find(partonFlavs.begin(), partonFlavs.end(), -4) != partonFlavs.end()){
        partonFlav = -4;
      }
      else{
        //Pick the smallest dR one
        int minElementIndex = std::min_element(partondRs.begin(), partondRs.end()) - partondRs.begin();
        partonFlav = partonFlavs[minElementIndex];
      }
      //std::cout << "FatJet " << nFatJet << " matched to parton with flavour " << partonFlav << "\n";
    }
    FatJet_partonFlav.push_back(partonFlav);
    nFatJet++;
  }

  nEGamma = 0;

  for (l1t::EGammaBxCollection::const_iterator egamma = egammasCollection->begin(0); egamma != egammasCollection->end(0); egamma++){
    EGamma_pt.push_back(egamma->pt());
    EGamma_eta.push_back(egamma->eta());
    EGamma_phi.push_back(egamma->phi());
    EGamma_e.push_back(egamma->energy());
    EGamma_Iso.push_back(egamma->hwIso());
    nEGamma++;
  }

  for (l1t::EtSumBxCollection::const_iterator etSum = etSumsCollection->begin(0); etSum != etSumsCollection->end(0); etSum++){
    switch (etSum->getType()){
      case l1t::EtSum::EtSumType::kTotalEt:
        totalEt = etSum->pt();
        break;
      
      case l1t::EtSum::EtSumType::kTotalHt:
        totalHt = etSum->pt();
        break;
      
      case l1t::EtSum::EtSumType::kMissingEt:
        missingEt = etSum->pt();
        missingEtPhi = etSum->phi();
        break;
      
      case l1t::EtSum::EtSumType::kMissingHt:
        missingHt = etSum->pt();
        missingHtPhi = etSum->phi();
        break;
      
      case l1t::EtSum::EtSumType::kTowerCount:
        towerCount = etSum->hwPt();
        break;
      
      default:
        continue;
    }
  }

  nRecoJet = 0;
  if(recoJets.isValid()){
    for (pat::JetCollection::const_iterator recoJet = recoJets->begin(); recoJet != recoJets->end(); recoJet++){
      if(recoJet->pt() < 10) continue;
      RecoJet_pt.push_back(recoJet->pt());
      RecoJet_eta.push_back(recoJet->eta());
      RecoJet_phi.push_back(recoJet->phi());
      RecoJet_e.push_back(recoJet->energy());
      RecoJet_mass.push_back(recoJet->mass());
      nRecoJet++;
    }
  }

  // Apply regression to L1 jets (already loaded in class initialisation)
  //auto ptRegression = std::make_unique<PtRegression>(model_path);

  /*
  //Only jet variables and event level sums
  for (UInt_t i=0; i < Jet_pt.size(); i++){
    float regressed_pt;
    float regressed_pt_scale;
    std::vector<float> inputs{Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_e[i], totalEt, missingEt, missingEtPhi, totalHt, missingHt, missingHtPhi, float(towerCount)};
    // Skip jets with pT < 30 and > 1000 GeV
    //if(Jet_pt[i] < 30. || Jet_pt[i] > 1000.){
    if(Jet_pt[i] > 1000.){
      regressed_pt = Jet_pt[i];
    }
    else{
      regressed_pt_scale = ptRegression->get_regressed_pt(inputs);
      regressed_pt = Jet_pt[i] * regressed_pt_scale;
      std::cout << "Jet " << i << " pT: " << Jet_pt[i] << " Regressed scale: " << regressed_pt_scale << " Regressed et: " << regressed_pt << "\n";
    }
    Jet_pt_regressed.push_back(regressed_pt);
  }
  */
  

  //Using jet, egamma, muon and event sum information
  
  for (UInt_t i=0; i<Jet_pt.size(); i++){
    float regressed_pt;
    float regressed_pt_scale;

    //Define inputs for isolation
    TLorentzVector EGamma_sum;
    TLorentzVector Muon_sum;
    EGamma_sum.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    Muon_sum.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    float Jet_relEGammaIso = 0;
    float Jet_relMuIso = 0;
    
    for(UInt_t j=0; j<EGamma_pt.size(); j++){
      TLorentzVector EGamma_vector;
      EGamma_vector.SetPtEtaPhiM(EGamma_pt[j], EGamma_eta[j], EGamma_phi[j], 0.0);
      float dphi = TVector2::Phi_mpi_pi(Jet_phi[i] - EGamma_phi[j]);
      float deta = Jet_eta[i] - EGamma_eta[j];
      float dR = std::sqrt(dphi*dphi + deta*deta);
      if(dR < 0.4){
        EGamma_sum += EGamma_vector;
      }
    }
    for(UInt_t j=0; j<Muon_pt.size(); j++){
      TLorentzVector Muon_vector;
      Muon_vector.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.1055);
      float dphi = TVector2::Phi_mpi_pi(Jet_phi[i] - Muon_phi[j]);
      float deta = Jet_eta[i] - Muon_eta[j];
      float dR = std::sqrt(dphi*dphi + deta*deta);
      if(dR < 0.4){
        Muon_sum += Muon_vector;
      }
    }

    float Jet_pt_log = TMath::Log(Jet_pt[i]);

    //Get relative isolation
    Jet_relEGammaIso = EGamma_sum.Pt() / Jet_pt[i];
    Jet_relMuIso = Muon_sum.Pt() / Jet_pt[i];

    //std::cout << "Jet pT: " << Jet_pt[i] << ", Jet pT log: " << Jet_pt_log << ", Jet relEGammaIso: " << Jet_relEGammaIso << ", Jet relMuIso: " << Jet_relMuIso << "\n";

    //Define inputs for regression
    std::vector<float> inputs{Jet_pt_log, Jet_eta[i], Jet_phi[i], Jet_relMuIso, Jet_relEGammaIso, totalEt, missingEt, missingEtPhi, totalHt, missingHt, missingHtPhi, float(towerCount)};

    //Skip jets with pT > 1000 GeV
    if(Jet_pt[i] > 1000.){
      regressed_pt = Jet_pt[i];
    }
    else{
      regressed_pt_scale = ptRegression->get_regressed_pt(inputs);
      regressed_pt = Jet_pt[i] * std::exp(regressed_pt_scale);
      //std::cout << "Jet " << i << " pT: " << Jet_pt[i] << " Regressed target scale: " << regressed_pt_scale << " Regressed et: " << regressed_pt << "\n";
    }
    Jet_pt_regressed.push_back(regressed_pt);
  }
  


  nRecoJetPuppi = 0;
  if(recoJetsPuppi.isValid()){
    for (pat::JetCollection::const_iterator recoJetPuppi = recoJetsPuppi->begin(); recoJetPuppi != recoJetsPuppi->end(); recoJetPuppi++){
      if(recoJetPuppi->pt() < 10) continue;
      RecoJetPuppi_pt.push_back(recoJetPuppi->pt());
      RecoJetPuppi_eta.push_back(recoJetPuppi->eta());
      RecoJetPuppi_phi.push_back(recoJetPuppi->phi());
      RecoJetPuppi_e.push_back(recoJetPuppi->energy());
      RecoJetPuppi_mass.push_back(recoJetPuppi->mass());
      nRecoJetPuppi++;
    }
  }

  if(recoMet.isValid()){
    RecoMet_pt = recoMet->begin()->pt();
    RecoMet_phi = recoMet->begin()->phi();
  }
  
  if(recoMetPuppi.isValid()){
    RecoMetPuppi_pt = recoMetPuppi->begin()->pt();
    RecoMetPuppi_phi = recoMetPuppi->begin()->phi();
  }

  tree->Fill();
}

void MCNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ method called once each job just before starting event loop  ------------
void MCNtuplizer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void MCNtuplizer::endJob() {
}

DEFINE_FWK_MODULE(MCNtuplizer);
