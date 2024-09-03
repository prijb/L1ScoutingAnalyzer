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
#include "DataFormats/Math/interface/deltaR.h"

// root include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

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

// Fast jet for additional clustering
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include <map>
#include <string>

using namespace l1ScoutingRun3;

class DataNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DataNtuplizer(const edm::ParameterSet&);
  ~DataNtuplizer() {}
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
    const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection
  );

  // Get list of bxs that have at least one object in them
  void getFilledBx(
    const edm::Handle<MuonOrbitCollection>& muonsCollection,
    const edm::Handle<JetOrbitCollection>& jetsCollection,
    const edm::Handle<EGammaOrbitCollection>& eGammasCollection,
    const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection,
    std::vector<unsigned> &filledBxVec
  );

  // tokens for scouting data
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Muon>> muonsTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Jet>> jetsTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::EGamma>> eGammasTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::BxSums>> bxSumsTokenData_;

  //Bool for accounting for online bx selection
  bool onlineSelection_;
  edm::EDGetTokenT<std::vector<unsigned>> selectedBxToken_;

  // Histograms
  std::map<std::string, TH1D*> m_1dhist_;
  
  // Tree that contains info per bunch crossing
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
  vector<Float16_t> Muon_ptUnconstrained;
  vector<Float16_t> Muon_etaAtVtx;
  vector<Float16_t> Muon_phiAtVtx;
  vector<Int_t> Muon_hwDXY;
  vector<Int_t> Muon_tfIndex;

  //et Sums
  float totalEt, totalHt, missingEt, missingEtPhi, missingHt, missingHtPhi;
  int towerCount;
 };

DataNtuplizer::DataNtuplizer(const edm::ParameterSet& iPSet)
  : muonsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("muonsTag"))),
    jetsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("jetsTag"))),
    eGammasTokenData_(consumes(iPSet.getParameter<edm::InputTag>("eGammasTag"))),
    bxSumsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("bxSumsTag")))
  {
  
  // the root file service to handle the output file
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("Events", "Events_bx");

  // Jets
  tree->Branch("orbit", &orbit, "orbit/I");
  tree->Branch("bx", &bxid, "bx/I");
  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("Jet_pt", &Jet_pt);
  tree->Branch("Jet_eta", &Jet_eta);
  tree->Branch("Jet_phi", &Jet_phi);
  tree->Branch("Jet_e", &Jet_e);
  tree->Branch("Jet_qual", &Jet_qual);

  // AK8 jets
  tree->Branch("nFatJet", &nFatJet, "nFatJet/I");
  tree->Branch("FatJet_pt", &FatJet_pt);
  tree->Branch("FatJet_eta", &FatJet_eta);
  tree->Branch("FatJet_phi", &FatJet_phi);
  tree->Branch("FatJet_e", &FatJet_e);
  tree->Branch("FatJet_m", &FatJet_m);
  
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
  tree->Branch("Muon_ptUnconstrained", &Muon_ptUnconstrained);
  tree->Branch("Muon_etaAtVtx", &Muon_etaAtVtx);
  tree->Branch("Muon_phiAtVtx", &Muon_phiAtVtx);
  tree->Branch("Muon_hwDXY", &Muon_hwDXY);
  tree->Branch("Muon_tfIndex", &Muon_tfIndex);

  // Sums
  tree->Branch("etSum", &totalEt);
  tree->Branch("htSum", &totalHt);
  tree->Branch("etMiss", &missingEt);
  tree->Branch("etMissPhi", &missingEtPhi);
  tree->Branch("htMiss", &missingHt);
  tree->Branch("htMissPhi", &missingHtPhi);
  tree->Branch("towerCount", &towerCount);
  
  // create histogram subdir
  TFileDirectory histoSubDir = fs->mkdir("histograms");
  
  //Some additional dimuon and dijet histograms
  m_1dhist_["DimuonMass"] = histoSubDir.make<TH1D>("DimuonMass", "Opposite charge dimuons", 15000, 0., 150.);
  m_1dhist_["DijetMass"] = histoSubDir.make<TH1D>("DijetMass", "Dijet mass", 200, 0., 1000.);
  m_1dhist_["DijetMassSel"] = histoSubDir.make<TH1D>("DijetMassSel", "Dijet mass (mu matched)", 1000, 0., 1000.);

  // Get list of bxs that pass selections
  if (onlineSelection_){
    selectedBxToken_ = consumes<std::vector<unsigned>>(iPSet.getParameter<edm::InputTag>("selectedBxTag"));
  }

}

void DataNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<OrbitCollection<l1ScoutingRun3::Muon>> muonsCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::Jet>> jetsCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::EGamma>> eGammasCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::BxSums>> bxSumsCollection; 

  iEvent.getByToken(muonsTokenData_, muonsCollection); 
  iEvent.getByToken(jetsTokenData_, jetsCollection); 
  iEvent.getByToken(eGammasTokenData_, eGammasCollection); 
  iEvent.getByToken(bxSumsTokenData_, bxSumsCollection); 

  // List of bunch crossings passing selections
  std::vector<unsigned> bxList;
  if (onlineSelection_){
    edm::Handle<std::vector<unsigned>> selectedBx;
    iEvent.getByToken(selectedBxToken_, selectedBx);
    bxList = *selectedBx;
  }
  else{
    getFilledBx(muonsCollection, jetsCollection, eGammasCollection, bxSumsCollection, bxList);
  }

  // get orbit number orbit
  unsigned orbitNum = iEvent.id().event();

  //Process all bunch crossings from the list
  for (const unsigned& bx : bxList) {
    processDataBx(
        orbitNum,
        bx,
        muonsCollection,
        jetsCollection,
        eGammasCollection,
        bxSumsCollection
    ); 
  }
  
  
 }

void DataNtuplizer::processDataBx(
  unsigned orbitNum,
  unsigned bx,
  const edm::Handle<MuonOrbitCollection>& muonsCollection,
  const edm::Handle<JetOrbitCollection>& jetsCollection,
  const edm::Handle<EGammaOrbitCollection>& eGammasCollection,
  const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection
  )
  {

  orbit = orbitNum;
  bxid = bx;

  // get iterator for the current BX
  const auto& jets = jetsCollection->bxIterator(bx);
  const auto& eGammas = eGammasCollection->bxIterator(bx);
  const auto& bxSums = bxSumsCollection->bxIterator(bx);
  const auto& muons = muonsCollection->bxIterator(bx);

  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_e.clear();
  Muon_qual.clear();
  Muon_hwCharge.clear();
  Muon_ptUnconstrained.clear();
  Muon_etaAtVtx.clear();
  Muon_phiAtVtx.clear();
  Muon_hwDXY.clear();
  Muon_tfIndex.clear();

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
  
  //Cluster AK8 jets using jets in fastjet
  vector<fastjet::PseudoJet> FatJetInputs;
  double FatJetR = 1.0;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, FatJetR);

  int muon_counter = 0;
  nMuon = muons.size();

  for (const auto& muon : muons) {
    l1t::Muon l1muon = getL1TMuon(muon);
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
    Jet_pt.push_back(l1jet.pt());
    Jet_eta.push_back(l1jet.eta());
    Jet_phi.push_back(l1jet.phi());
    Jet_e.push_back(l1jet.energy());
    Jet_qual.push_back(l1jet.hwQual());
    //Feed jets into fastjet for clustering
    TLorentzVector jet_vector;
    jet_vector.SetPtEtaPhiM(l1jet.pt(), l1jet.eta(), l1jet.phi(), 0.0);
    //Get fat jets
    FatJetInputs.push_back(fastjet::PseudoJet(jet_vector.Px(), jet_vector.Py(), jet_vector.Pz(), jet_vector.E()));


    jet_counter++;
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
    nFatJet++;
  }


  int egamma_counter = 0;
  nEGamma = eGammas.size();

  for (const auto& egamma : eGammas) {
    l1t::EGamma l1egamma = getL1TEGamma(egamma);
    EGamma_pt.push_back(l1egamma.pt());
    EGamma_eta.push_back(l1egamma.eta());
    EGamma_phi.push_back(l1egamma.phi());
    EGamma_e.push_back(l1egamma.energy());
    EGamma_Iso.push_back(l1egamma.hwIso());
    egamma_counter++;
  }


  // store some of the sums
  if (bxSums.size()>0){
    l1t::EtSum l1sum;

    l1sum = getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTotalEt);
    totalEt = l1sum.pt();

    l1sum = getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTotalHt);
    totalHt = l1sum.pt();

    l1sum = getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kMissingEt);
    missingEt = l1sum.pt();
    missingEtPhi = l1sum.phi();

    l1sum = getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kMissingHt);
    missingHt = l1sum.pt();
    missingHtPhi = l1sum.phi();

    l1sum = getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTowerCount);
    towerCount = l1sum.hwPt();
  }

  //Processing with filled elements (jet-muon matching etc)
  if (nMuon==2){
    math::PtEtaPhiMLorentzVector mu1(Muon_pt[0], Muon_eta[0], Muon_phi[0], 0.1055);
    math::PtEtaPhiMLorentzVector mu2(Muon_pt[1], Muon_eta[1], Muon_phi[1], 0.1055);
    if(Muon_hwCharge[0]!=Muon_hwCharge[1]){
      m_1dhist_["DimuonMass"]->Fill((mu1+mu2).M());
    }

  }

  if (nJet==2){
    math::PtEtaPhiMLorentzVector jet1(Jet_pt[0], Jet_eta[0], Jet_phi[0], 0.);
    math::PtEtaPhiMLorentzVector jet2(Jet_pt[1], Jet_eta[1], Jet_phi[1], 0.);
    m_1dhist_["DijetMass"]->Fill((jet1+jet2).M());
    
    //Matching jets with muons
    bool isMuMatched{false};
    if (nMuon>1){
      for (int i=0; i<nMuon; i++){
        math::PtEtaPhiMLorentzVector mu(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.1055);
        if(reco::deltaR(mu, jet1)<0.4 && reco::deltaR(mu, jet2)<0.4){
          isMuMatched = true;
          break;
        }
      }
    }
    //Add dijet mass if muon matched
    if (isMuMatched){
      m_1dhist_["DijetMassSel"]->Fill((jet1+jet2).M());
    }
  }

  // fill the tree
  tree->Fill();
}

// Get the list of filled bunch crossings for a given orbit
void DataNtuplizer::getFilledBx(
  const edm::Handle<MuonOrbitCollection>& muonsCollection,
  const edm::Handle<JetOrbitCollection>& jetsCollection,
  const edm::Handle<EGammaOrbitCollection>& eGammasCollection,
  const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection,
  std::vector<unsigned> &filledBxVec
  )
  {
  // Get the list of filled bunch crossings for each object
  std::vector<unsigned> muonBx = muonsCollection->getFilledBxs();
  std::vector<unsigned> jetBx = jetsCollection->getFilledBxs();
  std::vector<unsigned> eGammaBx = eGammasCollection->getFilledBxs();
  std::vector<unsigned> bxSumsBx = bxSumsCollection->getFilledBxs();

  // Get the intersection of the bunch crossings
  std::set<unsigned> filledBxSet;

  for(const unsigned & bx: muonsCollection->getFilledBxs()){
    filledBxSet.insert(bx);
  }
  for(const unsigned & bx: jetsCollection->getFilledBxs()){
    filledBxSet.insert(bx);
  }
  for(const unsigned & bx: eGammasCollection->getFilledBxs()){
    filledBxSet.insert(bx);
  }

  std::copy(filledBxSet.begin(), filledBxSet.end(), std::back_inserter(filledBxVec));
}

void DataNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
 
// ------------ method called once each job just before starting event loop  ------------
void DataNtuplizer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DataNtuplizer::endJob() {
}

DEFINE_FWK_MODULE(DataNtuplizer);
