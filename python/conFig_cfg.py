import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('skipEvents',
    default=0,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
options.register('processMode',
    default='GenLevel',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "process mode: Genlevel by default")
options.parseArguments()


process = cms.Process("FEVTAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v4')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      options.inputFiles
      )
    , skipEvents = cms.untracked.uint32(options.skipEvents)
    )
print (" >> Loaded",len(options.inputFiles),"input files from list.")
# process.fevt.mode = cms.string(options.processMode)
# print (" >> Processing as:",(process.fevt.mode))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
    )



# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#
# process.source = cms.Source("PoolSource",
#                                 fileNames = cms.untracked.vstring(
#             'file:/uscms/home/bbbam/nobackup/CMSSW_10_6_20/src/MCProduction/E2E-HToAATo4Tau/GEN_HToAAToTauTau_M13_2018UL.root'
#             # 'file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGeneration/gen_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/crab_gen_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/230510_053835/0000/GEN_HToAAToTauTau_M10_2018UL_90.root'
#                 )
#                             )

# process.TFileService = cms.Service("TFileService",
#                 fileName = cms.string("GenInfo_only_M13.root")
#                             )

process.fevt = cms.EDAnalyzer('GenAnalyzer',
   # genParticles    = cms.untracked.InputTag('genParticles',"","GEN"),
   genParticles    = cms.InputTag('genParticles',"",""),
   ak8GenJets    = cms.InputTag('ak8GenJets',"",""),
   genMetTrue    = cms.InputTag('genMetTrue',"",""),
                              )

process.p = cms.Path(process.fevt)
