// This code is a barebones EDAnalyzer to just test the reweighting output
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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
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

// Added classes
#include "L1ScoutingAnalyzer/L1ScoutingAnalyzer/interface/QCDWeightCalc.h"

class TestReweightAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TestReweightAnalyzer(const edm::ParameterSet&);
  ~TestReweightAnalyzer(){}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void beginJob() override;
    virtual void endJob() override;

    // Tokens for gen info
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PileupInfoToken_;
    
    // Tree 
    TTree* tree;

    // Weight
    float bxFreq = 30E6;
    std::string qcdWeightFile;
    QCDWeightCalc qcdWeightCalc;
    float weight;
    float weight_v2;
};

TestReweightAnalyzer::TestReweightAnalyzer(const edm::ParameterSet& iConfig):
    genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    PileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupInfo"))),
    qcdWeightFile(iConfig.getParameter<std::string>("qcdWeightFile")),
    qcdWeightCalc(iConfig.getParameter<std::string>("qcdWeightFile"), bxFreq)   
{
    edm::Service<TFileService> fs;
    fs->file().cd();
    tree = new TTree("putree","putree");
    //tree = fs->make<TTree>("putree","putree");

    tree->Branch("weight",&weight,"weight/F");
    tree->Branch("weight_v2",&weight_v2,"weight_v2/F");

    std::cout << "QCD filename: " << qcdWeightFile << std::endl;
}

void TestReweightAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_,genEventInfo);
    weight = 1.0;
    weight_v2 = 0.8;
    tree->Fill();
}

void TestReweightAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

void TestReweightAnalyzer::beginJob() {
}

void TestReweightAnalyzer::endJob() {
}

DEFINE_FWK_MODULE(TestReweightAnalyzer);