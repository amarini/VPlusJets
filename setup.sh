#!/bin/bash


#### Commands to be issued before checking out the code

#export "SCRAM_ARCH=slc5_amd64_gcc462"
#cmsrel CMSSW_5_3_10_patch2
#cd CMSSW_5_3_10_patch2/src
#cmsenv

#mkdir amarini; cd amarini; git clone git@github.com:amarini/VPlusJets -b V00-14;


get_git() ## usage get_git tag src dest
{
	cd $CMSSW_BASE/src
        tag=$1 && shift
        src=$1 && shift
        dst=$1 && shift
        wd=$(pwd)

        echo
        echo "Checking out $tag $src into $dst"

        if [[ ! -d $dst ]]; then
         git clone $src $dst
        fi
        cd $dst
        git fetch origin
        git checkout $tag
        
        cd $wd
        echo "done"
        echo
        
        packs="$packs $dst"
}


echo "Going into $CMSSW_BASE"
cd $CMSSW_BASE/src

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/SuperClusterFootprintRemoval --
#Random cone
# V01-05 bug fix for RC

get_git V01-05 https://github.com/peruzzim/SCFootprintRemoval.git PFIsolation/SCFootprintRemoval

#------------------PATCH TO HAVE 20 - 30 in RC isolation-----------------------------------------------#
cd $CMSSW_BASE/src/PFIsolation/SCFootprintRemoval/src
patch SuperClusterFootprintRemoval.cc <<EOF
--- SuperClusterFootprintRemoval.cc	2013-10-31 11:21:39.000000000 +0100
+++ /afs/cern.ch/user/a/amarini/work/ProductionAugust2013/CMSSW_5_3_9/src/PFIsolation/SCFootprintRemoval/src/SuperClusterFootprintRemoval.cc	2013-10-29 11:42:04.000000000 +0100
@@ -400,13 +400,13 @@
   bool found=false;
 
   for (reco::PFJetCollection::const_iterator jet=jetHandle->begin(); jet!=jetHandle->end(); jet++){
-    if (jet->pt()<20) continue;
+    if (jet->pt()<30) continue;
     float dR = reco::deltaR(eta,phi,jet->eta(),jet->phi());
     if (dR<mindR) found=true;
   }
 
   for (reco::PhotonCollection::const_iterator pho=photonHandle->begin(); pho!=photonHandle->end(); pho++){
-    if (pho->pt()<10) continue;
+    if (pho->pt()<20) continue;
     float dR = reco::deltaR(eta,phi,pho->eta(),pho->phi());
     if (dR<mindR) found=true;
   }
EOF
##################
cd $CMSSW_BASE/src

#these commands are needed by cms-cvs-history 
git init .
git fetch cmssw-main

## PU Jet ID
# taken from h2gglobe. Used by PU-ID

get_git V00-03-04 https://github.com/h2gglobe/External.git CMGTools/External 

##################

## Energy Regression by Josh
#get_git hggpaperV6 git@github.com:bendavid/GBRLikelihood HiggsAnalysis/GBRLikelihood
# after speaking with Josh he suggested the master. I put the sha1 of the
# current head. For EGTools, current head is still the tag.

get_git 20836b5 git@github.com:bendavid/GBRLikelihood HiggsAnalysis/GBRLikelihood
get_git hggpaperV6 git@github.com:bendavid/GBRLikelihoodEGTools HiggsAnalysis/GBRLikelihoodEGTools

########################### EG TOOLS  ########################
# are these related to H->gg variables? Globe is using head.
#cvs co -r V00-00-21 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
get_git master https://github.com/h2gglobe/EGammaAnalysisTools EGamma/EGammaAnalysisTools

git cms-cvs-history import V00-03-13 CommonTools/ParticleFlow
git cms-cvs-history import V15-01-11 RecoParticleFlow/PFProducer

############################ QG Tagger #######################
get_git v1-2-6 git@github.com:amarini/QuarkGluonTagger.git QuarkGluonTagger

scram b -j 4

cd amarini/VPlusJets

