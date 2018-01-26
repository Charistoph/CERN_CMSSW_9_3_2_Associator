# Setup

# Login
ssh cbernkop@lxplus.cern.ch

# Initiate Grid
voms-proxy-init -voms cms

# Initiate  CMSSW
cd CMSSW_9_3_2/src/
cmsenv
cd -

# Run associator_run_Cfg
cd CMSSW_9_3_2/src/
scram b
cd -
cmsRun associator_run_Cfg.py

# Access Output on Grid
root -l root://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/adamw/MultiElectron_FlatPt5To100/Run2_2017_phase1_2017_realistic_RAW_DIGI/171020_115233/0000/step2_4.root

# Copy Output files from Grid
xrdcp root://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/adamw/MultiElectron_FlatPt5To100/crab_MyMultiElectronFlatPt5To100_RECO/171129_085129/0000/output_gsf_associator_1.root .

# Test output of Associator in Python
start python in command line
import ROOT
tf = ROOT.TFile("output_gsf_associator.root")
tf.ls()
tf.Get("MyTrackAssociator").cd()
tree = ROOT.gDirectory.Get("track_associator_tree")
tree.Scan()

