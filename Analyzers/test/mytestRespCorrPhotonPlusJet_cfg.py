import FWCore.ParameterSet.Config as cms
process = cms.Process('ANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']
#process.GlobalTag.globaltag='START53_V7G::All' # latest for Summer12 MC

#process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets = process.kt6PFJets.clone(rParam = 0.6, doRhoFastjet = True)          

#load the response corrections calculator
process.load('HcalClosureTest.Analyzers.calcrespcorrphotonplusjet_cfi')

process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

# run over files
process.calcrespcorrphotonplusjet.rootHistFilename = cms.string('PhoJet_tree_nonCHS_120to170.root')

process.calcrespcorrphotonplusjet.doCaloJets = cms.bool(False)
process.calcrespcorrphotonplusjet.doPFJets = cms.bool(True)
process.calcrespcorrphotonplusjet.doGenJets = cms.bool(True)
process.calcrespcorrphotonplusjet.photonTriggers = cms.vstring(
    'HLT_Photon20_CaloIdVL_IsoL','HLT_Photon30_CaloIdVL_IsoL',
    'HLT_Photon50_CaloIdVL_IsoL','HLT_Photon75_CaloIdVL_IsoL',
    'HLT_Photon90_CaloIdVL_IsoL','HLT_Photon135',
    'HLT_Photon150','HLT_Photon160')
#process.calcrespcorrphotonplusjet.photonTriggers = cms.vstring()

# CMSSW 7_X_Y needs ak4 jets instead of ak5 jets!
process.calcrespcorrphotonplusjet.pfJetCollName = cms.string('ak5PFJets')
#process.calcrespcorrphotonplusjet.pfJetCollName = cms.string('ak5PFJetsL1') # L1 corrected jet? # needs production
process.calcrespcorrphotonplusjet.pfJetCorrName = cms.string('ak5PFL1FastL2L3')
process.calcrespcorrphotonplusjet.photonJetDPhiMin = cms.double(0.)


# Load file list
# Summer12_DR53X production G_Pt_XtoY
import FWCore.Utilities.FileUtils as FileUtils
#listFileName='fileinfo_GJet/makepy_Summer12_DR53X_G_Pt_170to300.txt'
#listFileName='selection_tmp.txt'
listFileName='selection_keepAlive.txt'
listFileName='fileInfo_pionGun.txt'
listFileName='fileInfo_RelVal_5_3_14_PhotonJets.txt'
mylist = FileUtils.loadListFromFile(listFileName)
#mylist.extend( FileUtils.loadListFromFile(listFileName) )
readFiles = cms.untracked.vstring( *mylist )

#readFiles = cms.untracked.vstring(
#  '/store/relval/CMSSW_5_3_14/RelValPhotonJets_Pt_10/GEN-SIM-RECO/START53_LV6-v1/00000/5847CE87-FB60-E311-A45A-0025905A6134.root',
#  '/store/relval/CMSSW_5_3_14/RelValPhotonJets_Pt_10/GEN-SIM-RECO/START53_LV6-v1/00000/C82EE6E4-4D60-E311-8621-0025905A4964.root'
#)

process.source = cms.Source("PoolSource", 
                            fileNames= readFiles
#                            fileNames= cms.untracked.vstring(
#        'file:selectionGPt_470to800_3k.root'
# )

)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.threshold             = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(1000)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# timing
#process.Timing = cms.Service('Timing')

process.p = cms.Path(process.calcrespcorrphotonplusjet)

# Not sure if needed
#process.output = \
#    cms.OutputModule("PoolOutputModule",
#                     outputCommands = cms.untracked.vstring("drop *"),
#                     fileName = cms.untracked.string('dummy_output.root'),
#                     )
#process.out = cms.EndPath(process.output)
