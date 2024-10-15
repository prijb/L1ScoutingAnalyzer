import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys

options = VarParsing.VarParsing ('analysis')

options.register  ("numOrbits",
                  -1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  "Number of events to process")

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

#Skip product if it's not found
#process.options = cms.untracked.PSet(
#  TryToContinue = cms.untracked.vstring('ProductNotFound')
#)

if options.isData:
  process.scNtuplizer = cms.EDAnalyzer("OfflineDataNtuplizer",
    muonsTag      = cms.InputTag("gmtStage2Digis", "Muon"),
    jetsTag       = cms.InputTag("caloStage2Digis", "Jet"),
    eGammasTag    = cms.InputTag("caloStage2Digis", "EGamma"),
    tausTag       = cms.InputTag("caloStage2Digis", "Tau"),
    etSumsTag     = cms.InputTag("caloStage2Digis", "EtSum"),
    recoJetsTag = cms.InputTag("slimmedJetsPuppi"),
    recoFatJetsTag = cms.InputTag("slimmedJetsAK8"),
    recoElectronsTag = cms.InputTag("slimmedElectrons"),
    recoPhotonsTag = cms.InputTag("slimmedPhotons"),
    recoMuonsTag = cms.InputTag("slimmedMuons"),
    recoMetTag = cms.InputTag("slimmedMETsPuppi"),
  )
else:
  print("Only Data is supported for this configuration")


process.p = cms.Path(
  process.scNtuplizer
)