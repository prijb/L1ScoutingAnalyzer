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

// root include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

// L1 scouting 
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCalo.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "L1TriggerScouting/Utilities/interface/conversion.h"
#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"

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

using namespace l1ScoutingRun3;

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DemoAnalyzer(const edm::ParameterSet&);
  ~DemoAnalyzer() {}
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

  void processDataBx(
    unsigned orbitNum,
    unsigned bx,
    const edm::Handle<MuonOrbitCollection>& muonsCollection,
    const edm::Handle<JetOrbitCollection>& jetsCollection,
    const edm::Handle<EGammaOrbitCollection>& eGammasCollection,
    const edm::Handle<TauOrbitCollection>& tausCollection,
    const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection
  );

  // tokens for scouting data
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Muon>> muonsTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Jet>> jetsTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::EGamma>> eGammasTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Tau>> tausTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::BxSums>> bxSumsTokenData_;

  // l1t standard data format
  std::vector<l1t::Jet> l1jets_;
  std::vector<l1t::EGamma> l1egs_;
  std::vector<l1t::EtSum> l1sums_;
  std::vector<l1t::Tau> l1taus_;
  std::vector<l1t::Muon> l1muons_;

  // map containing TH1D histograms
  std::map<std::string, TH1D*> m_1dhist_;

  //Tree that contains info per bunch crossing
  TTree* tree;

  Int_t orbit;
  Int_t bxid;

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

DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iPSet)
  : muonsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("muonsTag"))),
    jetsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("jetsTag"))),
    eGammasTokenData_(consumes(iPSet.getParameter<edm::InputTag>("eGammasTag"))),
    tausTokenData_(consumes(iPSet.getParameter<edm::InputTag>("tausTag"))),
    bxSumsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("bxSumsTag")))
  {

  // the root file service to handle the output file
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("Events", "Events_bx");

  tree->Branch("orbit", &orbit, "orbit/I");
  tree->Branch("bx", &bxid, "bx/I");
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
  
  // init internal containers for l1 objects
  l1muons_.reserve(8);
  l1jets_.reserve(12);
  l1egs_.reserve(12);
  l1taus_.reserve(12);
  l1sums_.reserve(12);
  
  // create histogram subdir
  TFileDirectory histoSubDir = fs->mkdir("histograms");

  // Init histograms
  m_1dhist_["MuonBxOcc"] = histoSubDir.make<TH1D>("MuonBxOcc", "BX in orbit with at least one muon", 3566, -0.5, 3565.5);
  m_1dhist_["Muon_JetVeto_BxOcc"] = histoSubDir.make<TH1D>("Muon_JetVeto_BxOcc", "BX in orbit with at least one muon and no jet in BX-1", 3566, -0.5, 3565.5);
  m_1dhist_["MuonPt"] = histoSubDir.make<TH1D>("MuonPt", "MuonPt", 512, 0, 256);
  m_1dhist_["numJetsBx_wMuon"] = histoSubDir.make<TH1D>("numJetsBx_wMuon", "Jet multiplicity in BX with >1Muon", 13, -0.5, 12.5);

}

void DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<OrbitCollection<l1ScoutingRun3::Muon>> muonsCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::Jet>> jetsCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::EGamma>> eGammasCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::Tau>> tausCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::BxSums>> bxSumsCollection; 

  iEvent.getByToken(muonsTokenData_, muonsCollection); 
  iEvent.getByToken(jetsTokenData_, jetsCollection); 
  iEvent.getByToken(eGammasTokenData_, eGammasCollection); 
  iEvent.getByToken(tausTokenData_, tausCollection); 
  iEvent.getByToken(bxSumsTokenData_, bxSumsCollection); 

  // get orbit number orbit
  unsigned orbitNum = iEvent.id().event();

  // process all BX in orbit containing at least a Muon
  // getFilledBxs() returns the list of filled BX in the muon orbit collection
  /*
  for (const unsigned& bx : muonsCollection->getFilledBxs()) {
    processDataBx(
        orbitNum,
        bx,
        muonsCollection,
        jetsCollection,
        eGammasCollection,
        tausCollection,
        bxSumsCollection
    ); 
  }
  */

  //process all BX in orbit containing at least a Jet
  for (const unsigned& bx : jetsCollection->getFilledBxs()) {
    processDataBx(
        orbitNum,
        bx,
        muonsCollection,
        jetsCollection,
        eGammasCollection,
        tausCollection,
        bxSumsCollection
    ); 
  }

  //process all BX regardless of object count
  /*
  processDataBx(
        orbitNum,
        bx,
        muonsCollection,
        jetsCollection,
        eGammasCollection,
        tausCollection,
        bxSumsCollection
    );
  */
  
 }

void DemoAnalyzer::processDataBx(
    unsigned orbitNum,
    unsigned bx,
    const edm::Handle<MuonOrbitCollection>& muonsCollection,
    const edm::Handle<JetOrbitCollection>& jetsCollection,
    const edm::Handle<EGammaOrbitCollection>& eGammasCollection,
    const edm::Handle<TauOrbitCollection>& tausCollection,
    const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection
  ) {

    orbit = orbitNum;
    bxid = bx;

    // get iterator for the current BX
    const auto& jets = jetsCollection->bxIterator(bx);
    const auto& eGammas = eGammasCollection->bxIterator(bx);
    const auto& taus = tausCollection->bxIterator(bx);
    const auto& bxSums = bxSumsCollection->bxIterator(bx);
    const auto& muons = muonsCollection->bxIterator(bx);

    // convert scouting objects to l1t::objects for semplicity
    // Note: Scouting objects are stored in hw quantities, if only a subset
    // of the features is needed, the function in L1TriggerScouting/Utilities/interface/conversion.h
    // can be used to get physical quantities.
    // for example, the momentum of a muon can be obtained with ugmt::fPt(hwPt);
    l1jets_.clear();
    l1egs_.clear();
    l1taus_.clear();
    l1sums_.clear();
    l1muons_.clear();

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


    int muon_counter = 0;
    nMuon = muons.size();

    for (const auto& muon : muons) {
      l1t::Muon l1muon = getL1TMuon(muon);
      l1muons_.emplace_back(getL1TMuon(muon));
      Muon_pt.push_back(l1muon.pt());
      Muon_eta.push_back(l1muon.eta());
      Muon_phi.push_back(l1muon.phi());
      Muon_e.push_back(l1muon.energy());
      Muon_qual.push_back(l1muon.hwQual());
      Muon_hwCharge.push_back(l1muon.hwCharge());
      Muon_etaAtVtx.push_back(l1muon.etaAtVtx());
      Muon_phiAtVtx.push_back(l1muon.phiAtVtx());
      Muon_hwDXY.push_back(l1muon.hwDXY());
      muon_counter++;
    }

    int jet_counter = 0;
    nJet = jets.size();

    for (const auto& jet : jets) {
      l1t::Jet l1jet = getL1TJet(jet);
      l1jets_.emplace_back(getL1TJet(jet));
      Jet_pt.push_back(l1jet.pt());
      Jet_eta.push_back(l1jet.eta());
      Jet_phi.push_back(l1jet.phi());
      Jet_e.push_back(l1jet.energy());
      Jet_qual.push_back(l1jet.hwQual());
      Jet_towerIEta.push_back(l1jet.towerIEta());
      Jet_towerIPhi.push_back(l1jet.towerIPhi());
      Jet_rawEt.push_back(l1jet.rawEt());
      Jet_seedEt.push_back(l1jet.seedEt());
      Jet_puEt.push_back(l1jet.puEt());
      //Jet_puDonutEt.push_back(l1jet.puDonutEt());
      jet_counter++;
    }
    for (const auto& egamma : eGammas) {
      l1egs_.emplace_back(getL1TEGamma(egamma));
    }
    for (const auto& tau : taus) {
      l1taus_.emplace_back(getL1TTau(tau));
    }    

    // store some of the sums
    if (bxSums.size()>0){
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTotalEt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTotalHt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kMissingEt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kMissingHt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTowerCount));
    }

    // fill the tree
    tree->Fill();

    // fill histograms
    m_1dhist_["MuonBxOcc"]->Fill(bx);

    for (const auto& muon: l1muons_){
      m_1dhist_["MuonPt"]->Fill(muon.pt());
    }

    // collections are sorted based on the object Pt. For exampel, the leading muon Pt
    // can be obtained with l1muons_[0].pt()
    
    // number of jets in bx
    m_1dhist_["numJetsBx_wMuon"]->Fill(jets.size());

    // require no jets in the previous BX
    // check that we are not in the first BX. There could be halo muons from non colliding
    // bunches / cosmics.
    if (bx>1){
      // we can access the jets (if any) with jetsCollection->bxIterator(bx-1);
      // for this example, we will just check how many of them there are in bx-1
      if (jetsCollection->getBxSize(bx-1)>0){
        m_1dhist_["Muon_JetVeto_BxOcc"]->Fill(bx); 

        // fill more hists, do something else....
      }
    }

}

void DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
 
// ------------ method called once each job just before starting event loop  ------------
void DemoAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoAnalyzer::endJob() {
}

DEFINE_FWK_MODULE(DemoAnalyzer);
