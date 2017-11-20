#! /bin/bash
#BSUB -q 1nh
#
echo "Application started"
cd $WORKDIR
#cd /tmp/cbernkop
cd /afs/cern.ch/user/c/cbernkop/Associator/CMSSW_9_3_2/src/
echo "`scramv1 runtime -sh`"
eval `scramv1 runtime -sh`
echo "`scram b`"
eval `scram b`
echo "1 complete: cmsenv"
cd $WORKDIR
cmsRun /afs/cern.ch/user/c/cbernkop/Associator/associator_run_Cfg.py
#python /afs/cern.ch/user/c/cbernkop/Associator/print_associator_tracks.py
echo "2 complete: cmsRun"
#mkdir /afs/cern.ch/work/c/cbernkop/condor_output/assoc_output_new
#echo "3 complete: mkdir"
cp output_gsf_associator.root /afs/cern.ch/work/c/cbernkop/condor_output/assoc_output_20171120/output_gsf_associator.root
echo "3 complete: cp output_gsf_associator.root"
cp step3.root /afs/cern.ch/work/c/cbernkop/condor_output/assoc_output_20171120/step3.root
echo "4 complete: cp step3.root"

