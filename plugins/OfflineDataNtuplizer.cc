//This looks at MINIAOD data
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

// HLT Scouting
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

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

class OfflineDataNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit OfflineDataNtuplizer(const edm::ParameterSet&);
  ~OfflineDataNtuplizer(){};
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

  //Tokens for L1T objects
  edm::EDGetTokenT<l1t::JetBxCollection> jetToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> egToken_;
  edm::EDGetTokenT<l1t::EtSumBxCollection> etSumToken_;
  edm::EDGetTokenT<l1t::TauBxCollection> tauToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection> muonToken_;
  //Tokens for Reco objects
  edm::EDGetTokenT<pat::JetCollection> recoJetToken_;
  edm::EDGetTokenT<pat::JetCollection> recoFatJetToken_;
  edm::EDGetTokenT<pat::ElectronCollection> recoElectronToken_;
  edm::EDGetTokenT<pat::PhotonCollection> recoPhotonToken_;
  edm::EDGetTokenT<pat::MuonCollection> recoMuonToken_;
  edm::EDGetTokenT<pat::METCollection> recoMetToken_;

  //Tree that contains info per event
  TTree* tree;

  //L1 info
  //Jets
  Int_t nJet;
  vector<Float16_t> Jet_pt;
  vector<Float16_t> Jet_eta;
  vector<Float16_t> Jet_phi;
  vector<Float16_t> Jet_e;
  vector<Int_t> Jet_qual;

  //AK8 jets made from L1 jets
  Int_t nFatJet;
  vector<Float16_t> FatJet_pt;
  vector<Float16_t> FatJet_eta;
  vector<Float16_t> FatJet_phi;
  vector<Float16_t> FatJet_e;
  vector<Float16_t> FatJet_m;

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

OfflineDataNtuplizer::OfflineDataNtuplizer(const edm::ParameterSet& iPset)
  :jetToken_(consumes<l1t::JetBxCollection>(iPset.getParameter<edm::InputTag>("jetsTag"))),
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

  edm::Service<TFileService> fs;

  //Create the TTree
  tree = fs->make<TTree>("Events", "Events");

  //L1 info
  tree->Branch("nJet", &nJet);
  tree->Branch("Jet_pt", &Jet_pt);
  tree->Branch("Jet_eta", &Jet_eta);
  tree->Branch("Jet_phi", &Jet_phi);
  tree->Branch("Jet_e", &Jet_e);
  tree->Branch("Jet_qual", &Jet_qual);

  tree->Branch("nFatJet", &nFatJet);
  tree->Branch("FatJet_pt", &FatJet_pt);
  tree->Branch("FatJet_eta", &FatJet_eta);
  tree->Branch("FatJet_phi", &FatJet_phi);
  tree->Branch("FatJet_e", &FatJet_e);
  tree->Branch("FatJet_m", &FatJet_m);

  tree->Branch("nEGamma", &nEGamma);
  tree->Branch("EGamma_pt", &EGamma_pt);
  tree->Branch("EGamma_eta", &EGamma_eta);
  tree->Branch("EGamma_phi", &EGamma_phi);
  tree->Branch("EGamma_e", &EGamma_e);
  tree->Branch("EGamma_Iso", &EGamma_Iso);    

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

  tree->Branch("etSum", &totalEt);
  tree->Branch("htSum", &totalHt);
  tree->Branch("etMiss", &missingEt);
  tree->Branch("etMissPhi", &missingEtPhi);
  tree->Branch("htMiss", &missingHt);
  tree->Branch("htMissPhi", &missingHtPhi);
  tree->Branch("towerCount", &towerCount);

  //Reco info
  tree->Branch("nRecoJet", &nRecoJet);
  tree->Branch("RecoJet_pt", &RecoJet_pt);
  tree->Branch("RecoJet_eta", &RecoJet_eta);
  tree->Branch("RecoJet_phi", &RecoJet_phi);
  tree->Branch("RecoJet_e", &RecoJet_e);
  tree->Branch("RecoJet_mass", &RecoJet_mass);
  tree->Branch("RecoJet_btagPNetB", &RecoJet_btagPNetB);
  tree->Branch("RecoJet_btagPNetCvL", &RecoJet_btagPNetCvL);

  tree->Branch("nRecoFatJet", &nRecoFatJet);
  tree->Branch("RecoFatJet_pt", &RecoFatJet_pt);
  tree->Branch("RecoFatJet_eta", &RecoFatJet_eta);
  tree->Branch("RecoFatJet_phi", &RecoFatJet_phi);
  tree->Branch("RecoFatJet_e", &RecoFatJet_e);
  tree->Branch("RecoFatJet_mass", &RecoFatJet_mass);

  tree->Branch("nRecoElectron", &nRecoElectron);
  tree->Branch("RecoElectron_pt", &RecoElectron_pt);
  tree->Branch("RecoElectron_eta", &RecoElectron_eta);
  tree->Branch("RecoElectron_phi", &RecoElectron_phi);
  tree->Branch("RecoElectron_e", &RecoElectron_e);
  tree->Branch("RecoElectron_charge", &RecoElectron_charge);
  
  tree->Branch("nRecoPhoton", &nRecoPhoton);
  tree->Branch("RecoPhoton_pt", &RecoPhoton_pt);
  tree->Branch("RecoPhoton_eta", &RecoPhoton_eta);
  tree->Branch("RecoPhoton_phi", &RecoPhoton_phi);
  tree->Branch("RecoPhoton_e", &RecoPhoton_e);
  
  tree->Branch("nRecoMuon", &nRecoMuon);
  tree->Branch("RecoMuon_pt", &RecoMuon_pt);
  tree->Branch("RecoMuon_eta", &RecoMuon_eta);
  tree->Branch("RecoMuon_phi", &RecoMuon_phi);
  tree->Branch("RecoMuon_e", &RecoMuon_e);
  tree->Branch("RecoMuon_charge", &RecoMuon_charge);

  tree->Branch("RecoMet_pt", &RecoMet_pt);
  tree->Branch("RecoMet_phi", &RecoMet_phi);  
}

void OfflineDataNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<l1t::JetBxCollection> jetsCollection;
  edm::Handle<l1t::EGammaBxCollection> egammasCollection;
  edm::Handle<l1t::EtSumBxCollection> etSumsCollection;
  edm::Handle<l1t::TauBxCollection> tausCollection;
  edm::Handle<l1t::MuonBxCollection> muonsCollection;
  //Reco
  edm::Handle<pat::JetCollection> recoJets;
  edm::Handle<pat::JetCollection> recoFatJets;
  edm::Handle<pat::ElectronCollection> recoElectrons;
  edm::Handle<pat::PhotonCollection> recoPhotons;
  edm::Handle<pat::MuonCollection> recoMuons;
  edm::Handle<pat::METCollection> recoMet;

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

  //L1 info
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_e.clear();
  Jet_qual.clear();

  FatJet_pt.clear();
  FatJet_eta.clear();
  FatJet_phi.clear();
  FatJet_e.clear();
  FatJet_m.clear();

  EGamma_pt.clear();
  EGamma_eta.clear();
  EGamma_phi.clear();
  EGamma_e.clear();
  EGamma_Iso.clear();

  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_e.clear();
  Muon_qual.clear();
  Muon_hwCharge.clear();
  Muon_etaAtVtx.clear();
  Muon_phiAtVtx.clear();
  Muon_hwDXY.clear();

  //Reco info
  RecoJet_pt.clear();
  RecoJet_eta.clear();
  RecoJet_phi.clear();
  RecoJet_e.clear();
  RecoJet_mass.clear();
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

  //Fill L1 info
  nJet = 0;
  
  for (l1t::JetBxCollection::const_iterator jet = jetsCollection->begin(0); jet != jetsCollection->end(0); ++jet) {
    Jet_pt.push_back(jet->pt());
    Jet_eta.push_back(jet->eta());
    Jet_phi.push_back(jet->phi());
    Jet_e.push_back(jet->energy());
    Jet_qual.push_back(jet->hwQual());
    
    //Fat jet clustering
    TLorentzVector jet_vector;
    jet_vector.SetPtEtaPhiM(jet->pt(), jet->eta(), jet->phi(), 0.0);
    FatJetInputs.push_back(fastjet::PseudoJet(jet_vector.Px(), jet_vector.Py(), jet_vector.Pz(), jet_vector.E()));
    
    nJet++;
  }

  //Cluster AK8 jets
  fastjet::ClusterSequence cs(FatJetInputs, jet_def);
  vector<fastjet::PseudoJet> FatJets = fastjet::sorted_by_pt(cs.inclusive_jets());

  nFatJet = 0;
  for (UInt_t i=0; i < FatJets.size(); i++){
    fastjet::PseudoJet fatjet = FatJets[i];
    auto constituents = fatjet.constituents();
    //Skip jets with less than 2 constituents
    if (constituents.size() < 2) continue;
    FatJet_pt.push_back(FatJets[i].pt());
    FatJet_eta.push_back(FatJets[i].eta());
    FatJet_phi.push_back(FatJets[i].phi());
    FatJet_e.push_back(FatJets[i].E());
    FatJet_m.push_back(FatJets[i].m());
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
  
  //Fill Reco info
  nRecoJet = 0;

  if(recoJets.isValid()){
    for (pat::JetCollection::const_iterator recoJet = recoJets->begin(); recoJet != recoJets->end(); recoJet++){
      RecoJet_pt.push_back(recoJet->pt());
      RecoJet_eta.push_back(recoJet->eta());
      RecoJet_phi.push_back(recoJet->phi());
      RecoJet_e.push_back(recoJet->energy());
      RecoJet_mass.push_back(recoJet->mass());
      
      //PNet scores
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
      
      nRecoJet++;
    }
  }

  nRecoFatJet = 0;

  if(recoFatJets.isValid()){
    for (pat::JetCollection::const_iterator recoFatJet = recoFatJets->begin(); recoFatJet != recoFatJets->end(); recoFatJet++){
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
    for (pat::ElectronCollection::const_iterator electron = recoElectrons->begin(); electron != recoElectrons->end(); electron++){
      RecoElectron_pt.push_back(electron->pt());
      RecoElectron_eta.push_back(electron->eta());
      RecoElectron_phi.push_back(electron->phi());
      RecoElectron_e.push_back(electron->energy());
      RecoElectron_charge.push_back(electron->charge());
      nRecoElectron++;
    }
  }

  nRecoPhoton = 0;

  if(recoPhotons.isValid()){
    for (pat::PhotonCollection::const_iterator photon = recoPhotons->begin(); photon != recoPhotons->end(); photon++){
      RecoPhoton_pt.push_back(photon->pt());
      RecoPhoton_eta.push_back(photon->eta());
      RecoPhoton_phi.push_back(photon->phi());
      RecoPhoton_e.push_back(photon->energy());
      nRecoPhoton++;
    }
  }

  nRecoMuon = 0;

  if(recoMuons.isValid()){
    for (pat::MuonCollection::const_iterator muon = recoMuons->begin(); muon != recoMuons->end(); muon++){
      RecoMuon_pt.push_back(muon->pt());
      RecoMuon_eta.push_back(muon->eta());
      RecoMuon_phi.push_back(muon->phi());
      RecoMuon_e.push_back(muon->energy());
      RecoMuon_charge.push_back(muon->charge());
      nRecoMuon++;
    }
  }

  if(recoMet.isValid()){
    RecoMet_pt = recoMet->begin()->pt();
    RecoMet_phi = recoMet->begin()->phi();
  }

  tree->Fill();
}

void OfflineDataNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ method called once each job just before starting event loop  ------------
void OfflineDataNtuplizer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void OfflineDataNtuplizer::endJob() {}

DEFINE_FWK_MODULE(OfflineDataNtuplizer);  