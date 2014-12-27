# copy a sample of data file

import FWCore.ParameterSet.Config as cms

process = cms.Process("MYCOPY")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

gjetsFiles = cms.untracked.vstring(
  '/store/relval/CMSSW_5_3_14/RelValPhotonJets_Pt_10/GEN-SIM-RECO/START53_LV6-v1/00000/5847CE87-FB60-E311-A45A-0025905A6134.root',
  '/store/relval/CMSSW_5_3_14/RelValPhotonJets_Pt_10/GEN-SIM-RECO/START53_LV6-v1/00000/C82EE6E4-4D60-E311-8621-0025905A4964.root'
)

qcdFiles = cms.untracked.vstring(
 '/store/relval/CMSSW_5_3_14/RelValQCD_FlatPt_15_3000/GEN-SIM-RECO/START53_LV3_Feb14-v1/00000/D2D13174-A495-E311-8372-0025905A60F4.root',
 '/store/relval/CMSSW_5_3_14/RelValQCD_FlatPt_15_3000/GEN-SIM-RECO/START53_LV3_Feb14-v1/00000/EAB6B353-9B95-E311-B000-00304867915A.root'
)

# Load file list
# Summer12_DR53X production G_Pt_XtoY
import FWCore.Utilities.FileUtils as FileUtils
listFileName='fileinfo_GJet/makepy_Summer12_DR53X_G_Pt_170to300.txt'
#listFileName='selection_tmp.txt'
mylist = FileUtils.loadListFromFile(listFileName)
#mylist.extend( FileUtils.loadListFromFile(listFileName) )
gjetsFiles = cms.untracked.vstring( *mylist )



process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(
#        '/store/mc/Summer12_DR53X/G_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/0035BB1E-950E-E211-B5FE-002481E101DA.root'
#        )
                            fileNames= gjetsFiles
                            #fileNames= qcdFiles
)

process.copyAll = cms.OutputModule("PoolOutputModule",
          fileName = cms.untracked.string("selection.root") )
process.printEventNumber = cms.OutputModule("AsciiOutputModule")
process.p = cms.EndPath(process.copyAll + process.printEventNumber)
