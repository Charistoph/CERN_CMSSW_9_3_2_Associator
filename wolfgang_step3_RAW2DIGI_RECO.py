# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent RECOSIM --runUnscheduled -s RAW2DIGI,RECO --datatier GEN-SIM-RECO --geometry DB:Extended --filein file:step2.root --fileout file:step3.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

#process.Tracer = cms.Service("Tracer")

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/adamw/MultiElectron_FlatPt5To100/Run2_2017_phase1_2017_realistic_RAW_DIGI/171020_115233/0000/step2_1.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.tpSelection = cms.EDFilter("TrackingParticleSelector",
    chargedOnly = cms.bool(True),
    intimeOnly = cms.bool(False),
    lip = cms.double(280),
    maxPhi = cms.double(3.2),
    maxRapidity = cms.double(2.5),
    minHit = cms.int32(0),
    minPhi = cms.double(-3.2),
    minRapidity = cms.double(-2.5),
    pdgId = cms.vint32(),
    ptMin = cms.double(1.0),
    ptMax = cms.double(999999.0),
    signalOnly = cms.bool(False),
    src = cms.InputTag("mix","MergedTrackTruth"),
    stableOnly = cms.bool(False),
    tip = cms.double(120)
)
process.contentAnalyzer = cms.EDAnalyzer("EventContentAnalyzer")
process.MyTrackAssociator = cms.EDAnalyzer("MyTrackAssociator",
    AbsoluteNumberOfHits = cms.bool(False),
    Cut_RecoToSim = cms.double(0.75),
    SimToRecoDenominator = cms.string('reco'), # either "sim" or "reco"
    Quality_SimToReco = cms.double(0.5),
    Purity_SimToReco = cms.double(0.75),
    ThreeHitTracksAreSpecial = cms.bool(True),
    PixelHitWeight = cms.double(1.0),
    # Das hier wird vom cmsRun nicht gefunden.
    associatePixel = cms.bool(True),
    associateStrip = cms.bool(True),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    useClusterTPAssociation = cms.bool(True),
    # Das stimmt sicher noch nicht.
    cluster2TPSrc = cms.InputTag("tpClusterProducer")
)
process.analyzer_step = cms.Path(#process.contentAnalyzer
    process.tpSelection+process.MyTrackAssociator)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
#    outputCommands = cms.untracked.vstring('drop *', 'keep *_electronGsfTracks_*_*'),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.TFileService = cms.Service("TFileService", fileName = cms.string("output_gsf_associator.root") )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step, process.analyzer_step,
                                process.endjob_step,process.RECOSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

from FWCore.ParameterSet.Utilities import convertToUnscheduled
process = convertToUnscheduled(process)
