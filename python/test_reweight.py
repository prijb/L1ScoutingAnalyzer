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


process.testReweight = cms.EDAnalyzer("TestReweightAnalyzer",
    genEventInfo = cms.InputTag("generator"),
    PileupInfo = cms.InputTag("slimmedAddPileupInfo"),
    PuWeight = cms.InputTag("stitchingWeight"),
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    qcdWeightFile = cms.string("data/qcd_with_weights.json")
)

process.p = cms.Path(
    process.testReweight
)