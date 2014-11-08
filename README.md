HcalClosureTest
===============

Code to calibrate CMS's HCAL in eta using dijet balance. Created by J. P. Chou (Brown) https://github.com/johnpaulchou/usercode/tree/master/HcalClosureTest and updated by David G. Sheffield (Rutgers).

Runs in CMSSW_5_3_18.

# Creating the tree

Running Analyzers/test/testRespCorrDiJets_cfg.py creates a tree containing the information needed for dijet balancing. It is set by default to run over CaloJets, PFJets, and GenJets. Any of these can be shut off by adding the lines

```
process.calcrespcorrdijets.doCaloJets = cms.bool(False)
process.calcrespcorrdijets.doPFJets   = cms.bool(False)
process.calcrespcorrdijets.doGenJets  = cms.bool(False)
```

When using PFJets, the dataset being run over must be RECO as AOD does not save the necessary PF block information. When using GenJets, the dataset must contain SIM information.

## Testing the results

The source code Analyzers/test/testRespCorrDiJetsTree.cc creates an executable CMSSW_5_3_18/test/slc5_amd64_gcc462/testRespCorrDiJetsTree that creates plots to validate the tree.

# Getting response corrections

The source code  DataFormat/test/runCaloJetCorr.cc creates an excecutable CMSSW_5_3_18/test/slc5_amd64_gcc462/runPFJetCorr that takes the tree and determines the response corrections.

## Testing the results

The source code DataFormat/test/testResponseCorr.cc creates an executable CMSSW_5_3_18/test/slc5_amd64_gcc462/testResponseCorr that creates plots of E<sub>reconstructed</sub>/E<sub>gen</sub> for response corrections of unity and response corrections found in the last step to compare their performance.
