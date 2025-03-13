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
#include "DataFormats/Common/interface/TriggerResults.h"

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
    edm::EDGetTokenT<double> PuWeightToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    
    // Tree 
    TTree* tree;

    // Weight
    float bxFreq = 30E6;
    std::string qcdWeightFile;
    QCDWeightCalc qcdWeightCalc;
    float weight;
    float weight_v2;
    std::vector<PileupSummaryInfo> pileupInfoIntime;
};

TestReweightAnalyzer::TestReweightAnalyzer(const edm::ParameterSet& iConfig):
    genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    PileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupInfo"))),
    PuWeightToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("PuWeight"))),
    triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
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
    // Clear vectors
    pileupInfoIntime.clear();

    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_,genEventInfo);
    // Add pileup info 
    edm::Handle<std::vector<PileupSummaryInfo>> PileupInfo;
    iEvent.getByToken(PileupInfoToken_,PileupInfo);

    weight = 1.0;
    weight_v2 = 0.8;

    // Retrieve the hard scale from GenEventInfo
    float hardPtHat = genEventInfo->qScale();

    // Loop over pileup info and store only those for bunch crossing == 0 (in-time)
    for (const auto& pu : *PileupInfo) {
        if (pu.getBunchCrossing() == 0) {
        pileupInfoIntime.push_back(pu);
        }
    }

    // Pileup pT hats (sorted in decreasing order)
    std::vector<float> puPtHats;
    if (pileupInfoIntime.size() > 0) {
        puPtHats = pileupInfoIntime[0].getPU_pT_hats();
    }
    std::sort(puPtHats.begin(), puPtHats.end(), std::greater<float>());

    // Debug: Replace puPtHats with a dummy singular event with pT hat = 1
    puPtHats.clear();
    //puPtHats.push_back(1);


    // Ignore triggers
    bool pass_em = false;
    bool pass_mu = false;

    // Calculate weights
    float calcweight = qcdWeightCalc.weight(hardPtHat, puPtHats, pass_em, pass_mu);
    weight_v2 = calcweight;

    // Debug print
    std::cout << "Hard pT hat: " << hardPtHat << std::endl;
    std::cout << puPtHats.size() << " PU pT hats: ";
    std::cout << "highest pT hat of " << puPtHats[0];
    //for (const auto& puPtHat : puPtHats) {
    //    std::cout << puPtHat << " ";
    //}
    std::cout << std::endl;
    std::cout << "Weight: " << calcweight << std::endl;


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