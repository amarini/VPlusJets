#!/bin/bash

#this file will print out a multicrab configuration to run on the full ReReco data

cat <<EOF
# section for multicrab: now has just the template crab.cfg, but more
# keys might appear in the future
[MULTICRAB]
#cfg=crab.cfg

# Section [COMMON] is common for all datasets
# General idea: you define all the parameter in the template (crab.cfg), 
# but you might want to change the template values for all dataset.
# The general syntax is that you first put the crab.cfg [SECTION] and
# the the crab.cfg [key], with a "." in between, exactly as you would do
# to pass to CRAB keys via command line.

[COMMON]

# This determines the direcory where the CRAB log files and CMSSW output files will go.
# It will be USER.ui_working_dir/section_name/
# where section_name is the corresponding  section "[xyz]" that you give below.
CRAB.jobtype = cmssw
CRAB.scheduler = remoteGlidein
CRAB.use_server = 0

CMSSW.output_file = PATZJetsExpress.root
#CMSSW.total_number_of_lumis = -1

USER.return_data = 0
USER.eMail = amarini@cern.ch
USER.copy_data = 1
USER.se_black_list = T2_US_Florida,T3_US_Colorado,T2_US_Nebraska,T2_EE
USER.ce_black_list = T2_US_Florida,T3_US_Colorado,T2_US_Nebraska,T2_EE
USER.ui_working_dir = Data
USER.storage_element = srm-eoscms.cern.ch
USER.storage_path=/srm/v2/server?SFN=/eos/cms/store/user/amarini/zjets_unmerged_V00-10_RC
#user_remote_dir=group/phys_higgs/amarini/MC/Summer12

CMSSW.lumis_per_job = 120
CMSSW.total_number_of_lumis = -1

EOF






for i in "/DoubleMuParked/Run2012C-22Jan2013-v1/AOD" "/DoubleMu/Run2012A-22Jan2013-v1/AOD" "/DoubleMuParked/Run2012D-22Jan2013-v1/AOD" "/DoubleMuParked/Run2012B-22Jan2013-v1/AOD" "/DoubleMuParked/Run2012B-22Jan2013-v1/AOD" "/DoubleElectron/Run2012B-22Jan2013-v1/AOD" "/DoubleElectron/Run2012C-22Jan2013-v1/AOD" "/DoubleElectron/Run2012D-22Jan2013-v1/AOD" "/DoubleElectron/Run2012A-22Jan2013-v1/AOD" "/MuEG/Run2012D-22Jan2013-v1/AOD" "/MuEG/Run2012A-22Jan2013-v1/AOD" "/MuEG/Run2012B-22Jan2013-v1/AOD" "/MuEG/Run2012C-22Jan2013-v1/AOD" "/Photon/Run2012A-22Jan2013-v1/AOD" "/SinglePhoton/Run2012B-22Jan2013-v1/AOD" "/SinglePhoton/Run2012C-22Jan2013-v1/AOD" "/SinglePhotonParked/Run2012D-22Jan2013-v1/AOD" ; do 

#echo $i ; 

echo "[$(echo $i | sed 's:/::' | sed 's:/:_:g')]"
echo "CMSSW.datasetpath=$i"
echo "CMSSW.lumi_mask = /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"
echo "CMSSW.pset = PATZJetsExpress_cfg.py"
echo "USER.user_remote_dir = $(echo $i | sed 's:/:_:g' | sed 's:_:/:' )"
echo 
echo


done
