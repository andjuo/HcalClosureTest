import FWCore.ParameterSet.Config as cms

process = cms.Process('ANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']

#load the response corrections calculator
process.load('HcalClosureTest.Analyzers.calcrespcorrdijets_cfi')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

# run over files

process.calcrespcorrdijets.rootHistFilename = cms.string('dijet_tree_5809400A.root')
process.calcrespcorrdijets.doCaloJets = cms.bool(False)
#process.calcrespcorrdijets.doPFJets = cms.bool(False)
#process.calcrespcorrdijets.doGenJets = cms.bool(False)
#process.calcrespcorrdijets.debug = cms.untracked.bool(True)
#process.calcrespcorrdijets.maxDeltaEta = cms.double(0.5)
#process.calcrespcorrdijets.minJetEt = cms.double(10.)
#process.calcrespcorrdijets.minSumJetEt = cms.double(20.)
#process.calcrespcorrdijets.maxThirdJetEt = cms.double(20.)

#import FWCore.Python.FileUtils as FileUtils
#readFiles = cms.untracked.vstring( FileUtils.loadListFromFile ('Pion_Pt-50.list') )
#process.source = cms.Source ("PoolSource",fileNames = readFiles)
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:selectionQCD.root'))
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:selectionGPt170to300_5809400A.root'))

#print readFiles

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(1000)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# timing
#process.Timing = cms.Service('Timing')

process.p = cms.Path(process.calcrespcorrdijets)
