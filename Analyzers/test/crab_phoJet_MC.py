from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'process_Summer12_GPt_120to170'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'mytestRespCorrPhotonPlusJet_cfg.py'
config.JobType.outputFiles = [ 'PhoJet_tree.root' ]

config.section_("Data")
#config.Data.inputDataset = '/GenericTTbar/HC-CMSSW_5_3_1_START53_V5-v1/GEN-SIM-RECO'
config.Data.inputDataset = '/G_Pt-120to170_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'

# how many files to analyze in 1 job
config.Data.unitsPerJob = 1

# not mandatory -- total number of units to process
config.Data.totalUnits = 1

#config.Data.outLFN = '/store/user/<subdir>' # or '/store/group/<subdir>'
config.Data.outLFN = '/store/user/andriusj/test'
config.Data.publication = False
config.Data.publishDataName = 'photonJet'

# If ignoreLocality=true, run anywhere, not just where the data is
config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T2_FR_CCIN2P3'
