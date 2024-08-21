import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys

options = VarParsing.VarParsing ('analysis')

options.register ("numOrbits",
                  -1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  "Number of orbits/events to process")

options.register ("inFile",
                  "file:",
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.string,
                  "Path to the input file")

options.register ("outFile",
                  "file:output.root",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Path of the output file")

options.register ("onlineSelection",
                  "",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Data online selection string")

options.register ("isData",
                  True,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "Whether input is data or MC")

options.parseArguments()

process = cms.Process( "SCNTUPLIZER" )

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(options.numOrbits)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(True)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
)

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(options.inFile)
)

#Skip duplicate check if MC
if not options.isData:
  process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# Choice of analyzer depends on whether the file is Data or MC
if options.isData:
  if options.onlineSelection != "":
    process.scNtuplizer = cms.EDAnalyzer("DataNtuplizer",
      muonsTag      = cms.InputTag("FinalBxSelector", "Muon"),
      jetsTag       = cms.InputTag("FinalBxSelector", "Jet"),
      eGammasTag    = cms.InputTag("FinalBxSelector", "EGamma"),
      bxSumsTag     = cms.InputTag("FinalBxSelector", "EtSum"),
      onlineSelection = cms.untracked.bool(True),
      selectedBxTag = cms.InputTag(options.onlineSelection, "SelBx")
    )
  else:
      process.scNtuplizer = cms.EDAnalyzer("DataNtuplizer",
      muonsTag      = cms.InputTag("l1ScGmtUnpacker", "Muon"),
      jetsTag       = cms.InputTag("l1ScCaloUnpacker", "Jet"),
      eGammasTag    = cms.InputTag("l1ScCaloUnpacker", "EGamma"),
      bxSumsTag     = cms.InputTag("l1ScCaloUnpacker", "EtSum"),
      onlineSelection = cms.untracked.bool(False)
    )

else:
  process.scNtuplizer = cms.EDAnalyzer("MCNtuplizer",
    genEventInfoTag        = cms.InputTag("generator"),
    genParticlesTag        = cms.InputTag("prunedGenParticles"),
    genJetsTag             = cms.InputTag("slimmedGenJets"),
    genFatJetsTag          = cms.InputTag("slimmedGenJetsAK8"),
    muonsTag      = cms.InputTag("gmtStage2Digis", "Muon"),
    jetsTag       = cms.InputTag("caloStage2Digis", "Jet"),
    eGammasTag    = cms.InputTag("caloStage2Digis", "EGamma"),
    tausTag       = cms.InputTag("caloStage2Digis", "Tau"),
    etSumsTag     = cms.InputTag("caloStage2Digis", "EtSum"),
    recoJetsTag   = cms.InputTag("slimmedJets"),
    recoJetsPuppiTag = cms.InputTag("slimmedJetsPuppi"),
    recoMetTag    = cms.InputTag("slimmedMETs"),
    recoMetPuppiTag = cms.InputTag("slimmedMETsPuppi"),
    #regressionPath = cms.FileInPath("L1ScoutingAnalyzer/L1ScoutingAnalyzer/data/jet_regression_model_new.model")
    regressionPath = cms.FileInPath("L1ScoutingAnalyzer/L1ScoutingAnalyzer/data/jet_regression_leading_subleading.model")
  )

process.p = cms.Path(
  process.scNtuplizer
)
