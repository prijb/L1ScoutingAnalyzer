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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

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
  edm::EDGetTokenT<pat::JetCollection> recoFatJetToken_;
  edm::EDGetTokenT<pat::ElectronCollection> recoElectronToken_;
  edm::EDGetTokenT<pat::PhotonCollection> recoPhotonToken_;
  edm::EDGetTokenT<pat::MuonCollection> recoMuonToken_;
  edm::EDGetTokenT<pat::METCollection> recoMetToken_;

  // Tree that contains info per bunch crossing
  TTree* tree;

  // Event info
  ULong64_t event;
  Int_t run;
  Int_t lumi;

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
  vector<Int_t> GenPart_statusFlags;

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

  //Reco jets (PUPPI)
  Int_t nRecoJet;
  vector<Float16_t> RecoJet_pt;
  vector<Float16_t> RecoJet_eta;
  vector<Float16_t> RecoJet_phi;
  vector<Float16_t> RecoJet_e;
  vector<Float16_t> RecoJet_mass;
  vector<Float16_t> RecoJet_partonFlav;
  //ParticleNet scores
  vector<Float16_t> RecoJet_btagPNetB;
  vector<Float16_t> RecoJet_btagPNetCvL;

  //Reco Fat jets (AK8)
  Int_t nRecoFatJet;
  vector<Float16_t> RecoFatJet_pt;
  vector<Float16_t> RecoFatJet_eta;
  vector<Float16_t> RecoFatJet_phi;
  vector<Float16_t> RecoFatJet_e;
  vector<Float16_t> RecoFatJet_mass;

  //Reco electrons
  Int_t nRecoElectron;
  vector<Float16_t> RecoElectron_pt;
  vector<Float16_t> RecoElectron_eta;
  vector<Float16_t> RecoElectron_phi;
  vector<Float16_t> RecoElectron_e;
  vector<Int_t> RecoElectron_charge;

  //Reco photons
  Int_t nRecoPhoton;
  vector<Float16_t> RecoPhoton_pt;
  vector<Float16_t> RecoPhoton_eta;
  vector<Float16_t> RecoPhoton_phi;
  vector<Float16_t> RecoPhoton_e;

  //Reco muons
  Int_t nRecoMuon;
  vector<Float16_t> RecoMuon_pt;
  vector<Float16_t> RecoMuon_eta;
  vector<Float16_t> RecoMuon_phi;
  vector<Float16_t> RecoMuon_e;
  vector<Int_t> RecoMuon_charge;

  //Reco MET
  Float16_t RecoMet_pt;
  Float16_t RecoMet_phi;

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
  recoFatJetToken_(consumes<pat::JetCollection>(iPset.getParameter<edm::InputTag>("recoFatJetsTag"))),
  recoElectronToken_(consumes<pat::ElectronCollection>(iPset.getParameter<edm::InputTag>("recoElectronsTag"))),
  recoPhotonToken_(consumes<pat::PhotonCollection>(iPset.getParameter<edm::InputTag>("recoPhotonsTag"))),
  recoMuonToken_(consumes<pat::MuonCollection>(iPset.getParameter<edm::InputTag>("recoMuonsTag"))),
  recoMetToken_(consumes<pat::METCollection>(iPset.getParameter<edm::InputTag>("recoMetTag")))
{

  // the root file service to handle the output file
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("Events", "Events_bx");

  // Event info
  tree->Branch("event", &event, "event/l");
  tree->Branch("run", &run, "run/I");
  tree->Branch("luminosityBlock", &lumi, "luminosityBlock/I");

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
  tree->Branch("GenPart_statusFlags", &GenPart_statusFlags);

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
  tree->Branch("RecoJet_partonFlav", &RecoJet_partonFlav);
  tree->Branch("RecoJet_btagPNetB", &RecoJet_btagPNetB);
  tree->Branch("RecoJet_btagPNetCvL", &RecoJet_btagPNetCvL);

  // Reco fat jet
  tree->Branch("nRecoFatJet", &nRecoFatJet, "nRecoFatJet/I");
  tree->Branch("RecoFatJet_pt", &RecoFatJet_pt);
  tree->Branch("RecoFatJet_eta", &RecoFatJet_eta);
  tree->Branch("RecoFatJet_phi", &RecoFatJet_phi);
  tree->Branch("RecoFatJet_e", &RecoFatJet_e);
  tree->Branch("RecoFatJet_mass", &RecoFatJet_mass);

  // Reco electrons
  tree->Branch("nRecoElectron", &nRecoElectron, "nRecoElectron/I");
  tree->Branch("RecoElectron_pt", &RecoElectron_pt);
  tree->Branch("RecoElectron_eta", &RecoElectron_eta);
  tree->Branch("RecoElectron_phi", &RecoElectron_phi);
  tree->Branch("RecoElectron_e", &RecoElectron_e);
  tree->Branch("RecoElectron_charge", &RecoElectron_charge);

  // Reco photons
  tree->Branch("nRecoPhoton", &nRecoPhoton, "nRecoPhoton/I");
  tree->Branch("RecoPhoton_pt", &RecoPhoton_pt);
  tree->Branch("RecoPhoton_eta", &RecoPhoton_eta);
  tree->Branch("RecoPhoton_phi", &RecoPhoton_phi);
  tree->Branch("RecoPhoton_e", &RecoPhoton_e);

  // Reco muons
  tree->Branch("nRecoMuon", &nRecoMuon, "nRecoMuon/I");
  tree->Branch("RecoMuon_pt", &RecoMuon_pt);
  tree->Branch("RecoMuon_eta", &RecoMuon_eta);
  tree->Branch("RecoMuon_phi", &RecoMuon_phi);
  tree->Branch("RecoMuon_e", &RecoMuon_e);
  tree->Branch("RecoMuon_charge", &RecoMuon_charge);

  // Reco MET
  tree->Branch("RecoMet_pt", &RecoMet_pt, "RecoMet_pt/F");
  tree->Branch("RecoMet_phi", &RecoMet_phi, "RecoMet_phi/F");

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
  edm::Handle<pat::JetCollection> recoFatJets;
  edm::Handle<pat::ElectronCollection> recoElectrons;
  edm::Handle<pat::PhotonCollection> recoPhotons;
  edm::Handle<pat::MuonCollection> recoMuons;
  edm::Handle<pat::METCollection> recoMet;

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
  iEvent.getByToken(recoFatJetToken_, recoFatJets);
  iEvent.getByToken(recoElectronToken_, recoElectrons);
  iEvent.getByToken(recoPhotonToken_, recoPhotons);
  iEvent.getByToken(recoMuonToken_, recoMuons);
  iEvent.getByToken(recoMetToken_, recoMet);

  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.luminosityBlock();

  GenPart_pt.clear();
  GenPart_eta.clear();
  GenPart_phi.clear();
  GenPart_e.clear();
  GenPart_mass.clear();
  GenPart_pdgId.clear();
  GenPart_charge.clear();
  GenPart_genPartIdxMother.clear();
  GenPart_MotherpdgId.clear();
  GenPart_statusFlags.clear();

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
  RecoJet_partonFlav.clear();
  RecoJet_btagPNetB.clear();
  RecoJet_btagPNetCvL.clear();

  RecoFatJet_pt.clear();
  RecoFatJet_eta.clear();
  RecoFatJet_phi.clear();
  RecoFatJet_e.clear();
  RecoFatJet_mass.clear();

  RecoElectron_pt.clear();
  RecoElectron_eta.clear();
  RecoElectron_phi.clear();
  RecoElectron_e.clear();
  RecoElectron_charge.clear();

  RecoPhoton_pt.clear();
  RecoPhoton_eta.clear();
  RecoPhoton_phi.clear();
  RecoPhoton_e.clear();

  RecoMuon_pt.clear();
  RecoMuon_eta.clear();
  RecoMuon_phi.clear();
  RecoMuon_e.clear();
  RecoMuon_charge.clear();

  RecoMet_pt = -999;
  RecoMet_phi = -999;

  //Cluster AK8 jets using jets in fastjet
  vector<fastjet::PseudoJet> FatJetInputs;
  double FatJetR = 1.0;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, FatJetR);

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
      GenPart_statusFlags.push_back(
            genParticle->statusFlags().isLastCopyBeforeFSR()                  * 16384 +
            genParticle->statusFlags().isLastCopy()                           * 8192  +
            genParticle->statusFlags().isFirstCopy()                          * 4096  +
            genParticle->statusFlags().fromHardProcessBeforeFSR()             * 2048  +
            genParticle->statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +
            genParticle->statusFlags().isHardProcessTauDecayProduct()         * 512   +
            genParticle->statusFlags().fromHardProcess()                      * 256   +
            genParticle->statusFlags().isHardProcess()                        * 128   +
            genParticle->statusFlags().isDirectHadronDecayProduct()           * 64    +
            genParticle->statusFlags().isDirectPromptTauDecayProduct()        * 32    +
            genParticle->statusFlags().isDirectTauDecayProduct()              * 16    +
            genParticle->statusFlags().isPromptTauDecayProduct()              * 8     +
            genParticle->statusFlags().isTauDecayProduct()                    * 4     +
            genParticle->statusFlags().isDecayedLeptonHadron()                * 2     +
            genParticle->statusFlags().isPrompt()                             * 1
      );
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
    nJet++;
  }

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
      //if(recoJet->pt() < 10) continue;
      RecoJet_pt.push_back(recoJet->pt());
      RecoJet_eta.push_back(recoJet->eta());
      RecoJet_phi.push_back(recoJet->phi());
      RecoJet_e.push_back(recoJet->energy());
      RecoJet_mass.push_back(recoJet->mass());

      //Simple dR matching with gen partons (for  jets)
      int partonFlav = 0;
      std::vector<int> partonFlavs{};
      std::vector<float> partondRs{};
      for(int i = 0; i < nGenPart; i++){
        //Ignore if not a quark
        if((std::abs(GenPart_pdgId[i]) > 5) && (GenPart_pdgId[i]!=21)) continue;
        float dphi = TVector2::Phi_mpi_pi(recoJet->phi() - GenPart_phi[i]);
        float deta = recoJet->eta() - GenPart_eta[i];
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
      RecoJet_partonFlav.push_back(partonFlav);
      float pNetbDisc = -1;
      float pNetcDisc = -1;
      if(recoJet->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsAll") > 0.){
        pNetbDisc = recoJet->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsAll");
      }
      if(recoJet->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsL") > 0.){
        pNetcDisc = recoJet->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsL");
      }
      RecoJet_btagPNetB.push_back(pNetbDisc);
      RecoJet_btagPNetCvL.push_back(pNetcDisc);
      /*
      //Testing b jet discriminator score
      float pNetbDisc = recoJet->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsAll");
      float pNetcDisc = recoJet->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsL");
      if(pNetbDisc > 0. || pNetcDisc > 0.){
        std::cout << "Reco jet: " << nRecoJet << "\n";
        std::cout << "Genmatched parton flavour " << partonFlav << "\n";
        std::cout << "Particle Net b vs udscg " << pNetbDisc << "\n";
        std::cout << "Particle Net c vs udsg " << pNetcDisc << "\n";
        std::cout << "\n";
      }
      */
      nRecoJet++;
    }
  }

  nRecoFatJet = 0;
  if(recoFatJets.isValid()){
    for (pat::JetCollection::const_iterator recoFatJet = recoFatJets->begin(); recoFatJet != recoFatJets->end(); recoFatJet++){
      //if(recoFatJet->pt() < 10) continue;
      RecoFatJet_pt.push_back(recoFatJet->pt());
      RecoFatJet_eta.push_back(recoFatJet->eta());
      RecoFatJet_phi.push_back(recoFatJet->phi());
      RecoFatJet_e.push_back(recoFatJet->energy());
      RecoFatJet_mass.push_back(recoFatJet->mass());
      nRecoFatJet++;
    }
  }

  nRecoElectron = 0;
  if(recoElectrons.isValid()){
    for (pat::ElectronCollection::const_iterator recoElectron = recoElectrons->begin(); recoElectron != recoElectrons->end(); recoElectron++){
      RecoElectron_pt.push_back(recoElectron->pt());
      RecoElectron_eta.push_back(recoElectron->eta());
      RecoElectron_phi.push_back(recoElectron->phi());
      RecoElectron_e.push_back(recoElectron->energy());
      RecoElectron_charge.push_back(recoElectron->charge());
      nRecoElectron++;
    }
  }

  nRecoPhoton = 0;
  if(recoPhotons.isValid()){
    for (pat::PhotonCollection::const_iterator recoPhoton = recoPhotons->begin(); recoPhoton != recoPhotons->end(); recoPhoton++){
      RecoPhoton_pt.push_back(recoPhoton->pt());
      RecoPhoton_eta.push_back(recoPhoton->eta());
      RecoPhoton_phi.push_back(recoPhoton->phi());
      RecoPhoton_e.push_back(recoPhoton->energy());
      nRecoPhoton++;
    }
  }

  nRecoMuon = 0;
  if(recoMuons.isValid()){
    for (pat::MuonCollection::const_iterator recoMuon = recoMuons->begin(); recoMuon != recoMuons->end(); recoMuon++){
      RecoMuon_pt.push_back(recoMuon->pt());
      RecoMuon_eta.push_back(recoMuon->eta());
      RecoMuon_phi.push_back(recoMuon->phi());
      RecoMuon_e.push_back(recoMuon->energy());
      RecoMuon_charge.push_back(recoMuon->charge());
      nRecoMuon++;
    }
  }

  if(recoMet.isValid()){
    RecoMet_pt = recoMet->begin()->pt();
    RecoMet_phi = recoMet->begin()->phi();
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
