#!/bin/bash

DIR=JSON
mkdir -p $DIR

for file in DoubleElectron_Run2012A-22Jan2013-v1_AOD_DONE.root DoubleElectron_Run2012A-22Jan2013-v1_AOD_missingLumis.root DoubleElectron_Run2012A-22Jan2013-v1_AOD_notDone.root DoubleElectron_Run2012B-22Jan2013-v1_AOD.root DoubleElectron_Run2012B-22Jan2013-v1_AOD_missingLumis.root DoubleElectron_Run2012C-22Jan2013-v1_AOD.root DoubleElectron_Run2012C-22Jan2013-v1_AOD_missingLumis.root DoubleElectron_Run2012D-22Jan2013-v1_AOD.root DoubleMuParked_Run2012B-22Jan2013-v1_AOD.root DoubleMuParked_Run2012C-22Jan2013-v1_AOD.root DoubleMuParked_Run2012D-22Jan2013-v1_AOD.root DoubleMu_Run2012A-22Jan2013-v1_AOD.root MuEG_Run2012A-22Jan2013-v1_AOD.root MuEG_Run2012B-22Jan2013-v1_AOD.root MuEG_Run2012C-22Jan2013-v1_AOD.root MuEG_Run2012C-22Jan2013-v1_AOD_missingLumis.root MuEG_Run2012D-22Jan2013-v1_AOD.root Photon_Run2012A-22Jan2013-v1_AOD_v2.root SinglePhotonParked_Run2012D-22Jan2013-v1_AOD_v2_DONE.root SinglePhotonParked_Run2012D-22Jan2013-v1_AOD_v2_missingLumis.root SinglePhotonParked_Run2012D-22Jan2013-v1_AOD_v2_notDone.root SinglePhoton_Run2012B-22Jan2013-v1_AOD.root SinglePhoton_Run2012C-22Jan2013-v1_AOD.root ; do

echo root://eoscms///store/user/amarini/zjets_V00-12/$file > /tmp/amarini/file.txt
echo Computing lumi for file $file

python lumi.py  -f /tmp/amarini/file.txt > $DIR/${file%%.root}.txt

done

