# Example L1Trigger Scouting Analyzer

This repo contains a simple example of an EDAnalyzer used to analyze L1T Scouting data collected during 2023.
It will produce a ROOT file containing 3 histograms:
1. The BX occupancy in the orbit
2. Pt distribution of the muons (example of how to manipulate objects)
3. The BX occupancy such that at least a muon is present and there are no jets in the previous BX (example of how to use objects from a different BX)

## 2023 L1 Scouting data
Data available in DAS, more info in [this](https://mattermost.web.cern.ch/cms-l1-scouting/channels/town-square) slide.
```
file dataset=/L1ScoutUGMTCALO/Run2023C-v1/RAW run=368636
```

Event content:
```
$ edmDumpEventContent root://cms-xrd-global.cern.ch//store/data/Run2023C/L1ScoutUGMTCALO/RAW/v1/000/368/636/00000/run368636_ls0400.root
Type                                  Module             Label     Process   
-----------------------------------------------------------------------------
OrbitCollection<l1ScoutingRun3::BxSums>    "CaloUnpacker"     ""        "SCPU"    
OrbitCollection<l1ScoutingRun3::EGamma>    "CaloUnpacker"     ""        "SCPU"    
OrbitCollection<l1ScoutingRun3::Jet>    "CaloUnpacker"     ""        "SCPU"    
OrbitCollection<l1ScoutingRun3::Muon>    "GmtUnpacker"      ""        "SCPU"    
OrbitCollection<l1ScoutingRun3::Tau>    "CaloUnpacker"     ""        "SCPU"
```

Each event contains data from one orbit, i.e. a collection of BX, and are stored in an [OrbitCollection](https://github.com/cms-sw/cmssw/blob/master/DataFormats/L1Scouting/interface/OrbitCollection.h).

## Instructions for the demo

```
# create a project area
cmsrel CMSSW_14_1_0_pre1 

cd CMSSW_14_1_0_pre1/src
cmsenv

mkdir Demo
cd Demo

git clone https://github.com/Mmiglio/L1ScoutingAnalyzer DemoAnalyzer

cd DemoAnalyzer
scram b

# init proxy
voms-proxy-init --voms cms

# process 10k orbits
cmsRun python/demo_cfg.py inFile=root://cms-xrd-global.cern.ch//store/data/Run2023C/L1ScoutUGMTCALO/RAW/v1/000/368/636/00000/run368636_ls0400.root outFile=test.root numOrbits=10000
```

