#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1I.h"

#include <stdio.h>
#include <stdlib.h>

#include <string>
using namespace std;

//#include "pileup_Data_v2.C"

int AddPUWeight(	const char * FileName="ntuple.root", //File name
			const char * HistoName="",//histo PU
			const char* Directory="accepted", //Directory name in the root file 
			const char* TreeName="events", // Tree Name in the directory Chosen
			int DeltaPU=0,
			int TrueDistr=1, // 1 for true 0 for false (observed)
			const char * PUFile="pileup/all.json.pileup.root",
			const char * branchName="PUWeight"
		   	)
{
	//declare a temporary string variable
	char str[255];	
	//open file
	TFile *f=TFile::Open(FileName,"UPDATE");
	if(f==NULL){
		fprintf(stderr,"%s: No such file or directory\n",FileName);
		return 1;
		}
	//cd to directory
	f->cd(Directory);

	//Constructing TreeName and open it
	sprintf(str,"%s/%s",Directory,TreeName);
	TTree *t=(TTree*)f->Get(str);
	if(t==NULL){
		fprintf(stderr,"%s: No such tree\n",TreeName);
		return 3;
		}
	//Set Existing branches
	double eventWeight;
	t->SetBranchAddress("eventWeight",&eventWeight);
	int puINT;
	if(TrueDistr)t->SetBranchAddress("puTrueINT",&puINT);
	else t->SetBranchAddress("puINT",&puINT);
	//Getting PU Profile
	TH1F *PU;
	if(HistoName[0]=='\0'){
		PU=new TH1F("PU","PU",100,0,100);
		PU->Sumw2();
		fprintf(stderr,"Generating Histogram from the tree\n");
		t->Draw("puINT>>PU","eventWeight","norm");
		}else{
		sprintf(str,"%s/%s",Directory,HistoName);
		fprintf(stderr,"Getting Histogram %s\n",HistoName);
		PU=(TH1F*)f->Get(str);
		if(PU==NULL)return 2;
		PU->Scale(1./PU->Integral("width"));
		}
     TH1D* pileup=NULL;
	if(DeltaPU==0) {TFile *f3=TFile::Open(  PUFile )    ;pileup=(TH1D*)f3->Get("pileup")->Clone();;} 
	if(DeltaPU>0)  {TFile *f3=TFile::Open(  PUFile )    ;pileup=(TH1D*)f3->Get("pileup")->Clone();;} 
	if(DeltaPU<0)  {TFile *f3=TFile::Open(  PUFile )    ;pileup=(TH1D*)f3->Get("pileup")->Clone();;} 
	if(pileup==NULL) {fprintf(stderr,"NO PILEUP FILE or HISTO %s \n",PUFile);}
	f->cd(Directory);
	fprintf(stderr,"Scaling pileup data\n");
   pileup->Scale(1./pileup->Integral("width"));
	//
	//Declaring variables
	double PUWeight;
	//Creating an empty branch in the tree
	TBranch *b;
	if(DeltaPU==0)	   { b=t->Branch( Form("%s",branchName)       ,&PUWeight,Form("%s/D",branchName)       );fprintf(stderr,"BranchName=%s\n",branchName)       ;}
	else if (DeltaPU>0){ b=t->Branch( Form("%sSysUp",branchName)  ,&PUWeight,Form("%sSysUp/D",branchName)  );fprintf(stderr,"BranchName=%sSysUp\n",branchName)  ;}
	else if (DeltaPU<0){ b=t->Branch( Form("%sSysDown",branchName),&PUWeight,Form("%sSysDown/D",branchName));fprintf(stderr,"BranchName=%sSysDown\n",branchName);}

	//Getting the Number of entries in the tree
	long long int NumberEntries=t->GetEntries();
	//looping on the entries in order to add the correct number of entries in the branch
	for(long long int i=0;i<NumberEntries;i++){
		t->GetEntry(i);
		if(PU->GetBinContent(PU->FindBin(puINT))==0) PUWeight=0;
		PUWeight = eventWeight*(double)pileup->GetBinContent(pileup->FindBin(puINT))/(double)PU->GetBinContent(PU->FindBin(   TMath::Max(puINT,0)   ));
		
		//cerr<<PUWeight<<endl;
		b->Fill();
		}
	//Write the Tree (With OverWrite Option)
	t->Write("",TObject::kOverwrite);
	//Close the file
	f->Close();
	//Print a message on stdout
	printf("Done\n");
	return 0;
}
