#!/bin/bash

## OLD
##pileupCalc.py -i json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt --inputLumiJSON pileup_JSON_DCSONLY_190389-208686_corr.txt --calcMode true --minBiasXsec 69400 --maxPileupBin 100 --numPileupBins 100  MyPileup.root

#compute JSON FILES

JSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"
PILEUP="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-208686_corr.txt"

filterJSON.py --min 190600  --max 197000 ${JSON} --output json/Cert_RunAB.json
filterJSON.py --min 197000  --max 203750 ${JSON} --output json/Cert_RunC.json
filterJSON.py --min 203750  --max 209000 ${JSON} --output json/Cert_RunD.json



for RDJSON in json/Cert_RunAB.json json/Cert_RunC.json json/Cert_RunD.json ; do

RUN="`echo $RDJSON | sed 's:json/Cert_::g' | sed 's:\.json::g'`"
pileupCalc.py -i $RDJSON --inputLumiJSON $PILEUP --calcMode true --minBiasXsec 69300 --maxPileupBin 100 --numPileupBins 100  MyPileup_${RUN}.root

pileupCalc.py -i $RDJSON --inputLumiJSON $PILEUP --calcMode true --minBiasXsec 72765 --maxPileupBin 100 --numPileupBins 100  MyPileup_${RUN}_UP.root

pileupCalc.py -i $RDJSON --inputLumiJSON $PILEUP --calcMode true --minBiasXsec 65835 --maxPileupBin 100 --numPileupBins 100  MyPileup_${RUN}_DN.root

pixelLumiCalc.py -i $RDJSON overview | tail -n 3

done

