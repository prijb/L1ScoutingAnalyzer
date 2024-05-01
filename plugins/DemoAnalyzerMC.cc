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

class DemoAnalyzerMC : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DemoAnalyzerMC(const edm::ParameterSet&);
  ~DemoAnalyzerMC() {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

  // Tokens for l1t objects
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<l1t::JetBxCollection> jetToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> egToken_;
  edm::EDGetTokenT<l1t::EtSumBxCollection> etSumToken_;
  edm::EDGetTokenT<l1t::TauBxCollection> tauToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection> muonToken_;

  // Tree that contains info per bunch crossing
  TTree* tree;

  Float_t genPtHat;
  
  //Jets
  Int_t nJet;
  vector<Float16_t> Jet_pt;
  vector<Float16_t> Jet_eta;
  vector<Float16_t> Jet_phi;
  vector<Float16_t> Jet_e;
  vector<Int_t> Jet_qual;
  vector<Int_t> Jet_towerIEta;
  vector<Int_t> Jet_towerIPhi;
  vector<Int_t> Jet_rawEt;
  vector<Int_t> Jet_seedEt;
  vector<Int_t> Jet_puEt;
  //vector<Int_t[4]> Jet_puDonutEt;

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

};

DemoAnalyzerMC::DemoAnalyzerMC(const edm::ParameterSet& iPset)
  :genEventInfoToken_(consumes<GenEventInfoProduct>(iPset.getParameter<edm::InputTag>("genEventInfoTag"))), 
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
  tree->Branch("genPtHat", &genPtHat, "genPtHat/F");
  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("Jet_pt", &Jet_pt);
  tree->Branch("Jet_eta", &Jet_eta);
  tree->Branch("Jet_phi", &Jet_phi);
  tree->Branch("Jet_e", &Jet_e);
  tree->Branch("Jet_qual", &Jet_qual);
  tree->Branch("Jet_towerIEta", &Jet_towerIEta);
  tree->Branch("Jet_towerIPhi", &Jet_towerIPhi);
  tree->Branch("Jet_rawEt", &Jet_rawEt);
  tree->Branch("Jet_seedEt", &Jet_seedEt);
  tree->Branch("Jet_puEt", &Jet_puEt);
  //tree->Branch("Jet_puDonutEt", &Jet_puDonutEt);
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

}

void DemoAnalyzerMC::analyze(const edm::Event& iEvent, const edm::EventSetup&){
  edm::Handle<GenEventInfoProduct> genEventInfo;
  edm::Handle<l1t::JetBxCollection> jetsCollection;
  edm::Handle<l1t::EGammaBxCollection> egammasCollection;
  edm::Handle<l1t::EtSumBxCollection> etsumsCollection;
  edm::Handle<l1t::TauBxCollection> tausCollection;
  edm::Handle<l1t::MuonBxCollection> muonsCollection;

  iEvent.getByToken(genEventInfoToken_, genEventInfo);
  iEvent.getByToken(jetToken_, jetsCollection);
  iEvent.getByToken(egToken_, egammasCollection);
  iEvent.getByToken(etSumToken_, etsumsCollection);
  iEvent.getByToken(tauToken_, tausCollection);
  iEvent.getByToken(muonToken_, muonsCollection);

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
  Jet_towerIEta.clear();
  Jet_towerIPhi.clear();
  Jet_rawEt.clear();
  Jet_seedEt.clear();
  Jet_puEt.clear();
  //Jet_puDonutEt.clear();

  //Gen event pT hat
  genPtHat = 0;
  if(genEventInfo.isValid()){
    genPtHat = genEventInfo->binningValues()[0];
  }

  nMuon = 0;

  for (const auto& muon : *muonsCollection) {
    Muon_pt.push_back(muon.pt());
    Muon_eta.push_back(muon.eta());
    Muon_phi.push_back(muon.phi());
    Muon_e.push_back(muon.energy());
    Muon_qual.push_back(muon.hwQual());
    Muon_hwCharge.push_back(muon.hwCharge());
    Muon_etaAtVtx.push_back(muon.etaAtVtx());
    Muon_phiAtVtx.push_back(muon.phiAtVtx());
    Muon_hwDXY.push_back(muon.hwDXY());
    nMuon++;
  }

  nJet = 0;

  for (const auto& jet : *jetsCollection) {
    Jet_pt.push_back(jet.pt());
    Jet_eta.push_back(jet.eta());
    Jet_phi.push_back(jet.phi());
    Jet_e.push_back(jet.energy());
    Jet_qual.push_back(jet.hwQual());
    Jet_towerIEta.push_back(jet.towerIEta());
    Jet_towerIPhi.push_back(jet.towerIPhi());
    Jet_rawEt.push_back(jet.rawEt());
    Jet_seedEt.push_back(jet.seedEt());
    Jet_puEt.push_back(jet.puEt());
    //Jet_puDonutEt.push_back(jet.puDonutEt());
    nJet++;
  }

  tree->Fill();
}

void DemoAnalyzerMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ method called once each job just before starting event loop  ------------
void DemoAnalyzerMC::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoAnalyzerMC::endJob() {
}

DEFINE_FWK_MODULE(DemoAnalyzerMC);