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
USER.ui_working_dir = Mc
#USER.storage_element = srm-eoscms.cern.ch
#USER.storage_path=/srm/v2/server?SFN=/eos/cms/store/user/amarini/zjets_V00-14
USER.storage_element = T2_CH_CSCS 

CMSSW.number_of_jobs = 2500
CMSSW.total_number_of_events = -1

EOF

#for i in "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/DYJetsToLL_PtZ-100_TuneZ2star_8TeV_ext-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM" "/DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV_ext-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM" "/DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV_ext-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM" "/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM" "/DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa/Summer12_DR53X-PU_S10_START53_V7C-v2/AODSIM" "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_RD1_START53_V7N-v1/AODSIM"
for i in "/GJets_HT-40To100_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM" "/GJets_HT-100To200_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM" "/GJets_HT-200To400_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/GJets_HT-400ToInf_8TeV-madgraph_v3/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM"
#for i in "/QCD_Pt_20_30_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_30_80_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_80_170_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_170_250_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_250_350_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_350_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM" "/QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU25bx25_START53_V19D-v1/AODSIM" "/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU25bx50_START53_V19D-v1/AODSIM" "/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU25bx25_START53_V19D-v1/AODSIM" "/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU25bx50_START53_V19D-v1/AODSIM" "/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM" "/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU25bx25_START53_V19D-v1/AODSIM" "/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU25bx50_START53_V19D-v1/AODSIM" "/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM"
#for i in "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM"
do
#echo $i ; 

echo "[$(echo $i | sed 's:/::' | sed 's:/:_:g')]"
echo "CMSSW.datasetpath=$i"
if echo "$i" | grep "RD1" &>/dev/null ; then
	echo "CMSSW.pset = PATZJetsExpressRD_cfg.py"
else
	echo "CMSSW.pset = PATZJetsExpress_cfg.py"
fi
#echo "USER.user_remote_dir = $(echo $i | sed 's:/:_:g' | sed 's:_:/:' )"
echo "USER.user_remote_dir = zjets_V00_14/$(echo $i | sed 's:/:_:g' | sed 's:_:/:' )"
echo 
echo

done
