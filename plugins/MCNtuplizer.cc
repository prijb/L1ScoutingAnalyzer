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

// root include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

// l1trigger
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

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

  // Tokens for l1t objects
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  edm::EDGetTokenT<l1t::JetBxCollection> jetToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> egToken_;
  edm::EDGetTokenT<l1t::EtSumBxCollection> etSumToken_;
  edm::EDGetTokenT<l1t::TauBxCollection> tauToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection> muonToken_;

  // Tree that contains info per bunch crossing
  TTree* tree;

  //Gen info
  Float_t genPtHat;
  Int_t nGenPart;
  Int_t nGenJet;
  vector<Float16_t> GenPart_pt;
  vector<Float16_t> GenPart_eta;
  vector<Float16_t> GenPart_phi;
  vector<Float16_t> GenPart_e;
  vector<Float16_t> GenPart_mass;
  vector<Int_t> GenPart_pdgId;
  vector<Int_t> GenPart_charge;
  
  vector<Float16_t> GenJet_pt;
  vector<Float16_t> GenJet_eta;
  vector<Float16_t> GenJet_phi;
  vector<Float16_t> GenJet_e;
  
  //Jets
  Int_t nJet;
  vector<Float16_t> Jet_pt;
  vector<Float16_t> Jet_eta;
  vector<Float16_t> Jet_phi;
  vector<Float16_t> Jet_e;
  vector<Int_t> Jet_qual;
  //vector<Int_t[4]> Jet_puDonutEt;

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

};

MCNtuplizer::MCNtuplizer(const edm::ParameterSet& iPset)
  :genEventInfoToken_(consumes<GenEventInfoProduct>(iPset.getParameter<edm::InputTag>("genEventInfoTag"))), 
  genParticlesToken_(consumes<reco::GenParticleCollection>(iPset.getParameter<edm::InputTag>("genParticlesTag"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iPset.getParameter<edm::InputTag>("genJetsTag"))),
  jetToken_(consumes<l1t::JetBxCollection>(iPset.getParameter<edm::InputTag>("jetsTag"))),
  egToken_(consumes<l1t::EGammaBxCollection>(iPset.getParameter<edm::InputTag>("eGammasTag"))),
  etSumToken_(consumes<l1t::EtSumBxCollection>(iPset.getParameter<edm::InputTag>("etSumsTag"))),
  tauToken_(consumes<l1t::TauBxCollection>(iPset.getParameter<edm::InputTag>("tausTag"))),
  muonToken_(consumes<l1t::MuonBxCollection>(iPset.getParameter<edm::InputTag>("muonsTag")))
{

  // the root file service to handle the output file
  edm::Service<TFileService> fs;

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

  tree->Branch("nGenJet", &nGenJet, "nGenJet/I");
  tree->Branch("GenJet_pt", &GenJet_pt);
  tree->Branch("GenJet_eta", &GenJet_eta);
  tree->Branch("GenJet_phi", &GenJet_phi);
  tree->Branch("GenJet_e", &GenJet_e);
  
  // Jets
  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("Jet_pt", &Jet_pt);
  tree->Branch("Jet_eta", &Jet_eta);
  tree->Branch("Jet_phi", &Jet_phi);
  tree->Branch("Jet_e", &Jet_e);
  tree->Branch("Jet_qual", &Jet_qual);

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

}

void MCNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup&){
  edm::Handle<GenEventInfoProduct> genEventInfo;
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<l1t::JetBxCollection> jetsCollection;
  edm::Handle<l1t::EGammaBxCollection> egammasCollection;
  edm::Handle<l1t::EtSumBxCollection> etSumsCollection;
  edm::Handle<l1t::TauBxCollection> tausCollection;
  edm::Handle<l1t::MuonBxCollection> muonsCollection;

  iEvent.getByToken(genEventInfoToken_, genEventInfo);
  iEvent.getByToken(genParticlesToken_, genParticles);
  iEvent.getByToken(genJetsToken_, genJets);
  iEvent.getByToken(jetToken_, jetsCollection);
  iEvent.getByToken(egToken_, egammasCollection);
  iEvent.getByToken(etSumToken_, etSumsCollection);
  iEvent.getByToken(tauToken_, tausCollection);
  iEvent.getByToken(muonToken_, muonsCollection);

  GenPart_pt.clear();
  GenPart_eta.clear();
  GenPart_phi.clear();
  GenPart_e.clear();
  GenPart_mass.clear();
  GenPart_pdgId.clear();
  GenPart_charge.clear();

  GenJet_pt.clear();
  GenJet_eta.clear();
  GenJet_phi.clear();
  GenJet_e.clear();

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

  EGamma_pt.clear();
  EGamma_eta.clear();
  EGamma_phi.clear();
  EGamma_e.clear();
  EGamma_Iso.clear();


  //Gen event pT hat
  genPtHat = (genEventInfo->hasBinningValues() ? genEventInfo->binningValues()[0] : 0.0);

  //Gen particles
  nGenPart = 0;
  if(genParticles.isValid()){
    for (reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
      //Skip particles that don't pass pT cut (>5 unless muons for which >3 is set)
      bool fill = false;
      if(abs(genParticle->pdgId()) == 13){
        if(genParticle->pt() > 3) fill = true;
      }
      else{
        if(genParticle->pt() > 5) fill = true;
      }
      if(!fill) continue;
      GenPart_pt.push_back(genParticle->pt());
      GenPart_eta.push_back(genParticle->eta());
      GenPart_phi.push_back(genParticle->phi());
      GenPart_e.push_back(genParticle->energy());
      GenPart_mass.push_back(genParticle->mass());
      GenPart_pdgId.push_back(genParticle->pdgId());
      GenPart_charge.push_back(genParticle->charge());
      nGenPart++;
    }
  }

  //Gen jets
  nGenJet = 0;
  if(genJets.isValid()){
    for (reco::GenJetCollection::const_iterator genJet = genJets->begin(); genJet != genJets->end(); genJet++){
      //Only store gen jets with pT > 10 GeV
      if(genJet->pt() < 10) continue;
      GenJet_pt.push_back(genJet->pt());
      GenJet_eta.push_back(genJet->eta());
      GenJet_phi.push_back(genJet->phi());
      GenJet_e.push_back(genJet->energy());
      nGenJet++;
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
    nJet++;
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