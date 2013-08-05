#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1I.h"
//unit definition of femtobarn vs picobarn: it may be useful to have plot normalized to 1fb
#define FB 1000.
#define MAXRUN 10

#include <stdio.h>
#include <stdlib.h>

class Run{
public:
	Run(){PU=NULL;PU_UP=NULL;PU_DN=NULL;PU_actual=NULL;};
	~Run(){};//{if(PU) delete PU;if(PU_UP) delete PU_UP; if(PU_DN)delete PU_DN;if(PU_actual) delete PU_actual;};
	long min;
	long max;
	float nEventsProcessed;
	TH1D * PU;
	TH1D * PU_UP;
	TH1D * PU_DN;
	TH1D * PU_actual;
	float CrossSection;
	float lumi;
	float weight(){ 
			return lumi*CrossSection*FB/nEventsProcessed; 
			};
	void Normalize() { 
			PU->Scale(1./PU->Integral("width"));
			PU_UP->Scale(1./PU_UP->Integral("width"));
			PU_DN->Scale(1./PU_DN->Integral("width"));
			PU_actual->Scale(1./PU_actual->Integral("width"));
			};
    	float puReweight(int puTrueInt,int syst=0) {
				
				if(syst==0)
					return PU->GetBinContent(PU->FindBin(puTrueInt)) / PU_actual->GetBinContent( PU_actual->FindBin(puTrueInt) );
				else if(syst>0)
					return PU_UP->GetBinContent(PU->FindBin(puTrueInt)) / PU_actual->GetBinContent( PU_actual->FindBin(puTrueInt) );
				else if(syst<0)
					return PU_DN->GetBinContent(PU->FindBin(puTrueInt)) / PU_actual->GetBinContent( PU_actual->FindBin(puTrueInt) );
				};


};

int AddRDWeight(	const char * FileName="ntuple.root", //File name
			const char* Directory="accepted", //Directory name in the root file 
			const char* TreeName="events", // Tree Name in the directory Chosen
			const double CrossSection=1300, // Cross Section of the Sample
			const char* ProcessedTreeName="processedData", // Tree Name in the directory Chosen with all data processed
			const char * configRun="RunAndLumi.txt",
			const float kFactor=1.0
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
	//Minimum check on Cross Section
	if(CrossSection<0.0){
		fprintf(stderr,"Cross Section (%lf) <0\n",CrossSection);
		return 2;
		}
	//Constructing TreeName and open it
	sprintf(str,"%s/%s",Directory,TreeName);
	TTree *t=(TTree*)f->Get(str);
	if(t==NULL){
		fprintf(stderr,"%s: No such tree\n",TreeName);
		return 3;
		}
	//Constructing HistoName and open it 
	sprintf(str,"%s/%s",Directory,ProcessedTreeName);
	TTree *processed=(TTree*)f->Get(str);
	if(processed==NULL){
		fprintf(stderr,"%s: No such tree\n",ProcessedTreeName);
		return 4;
		}
	vector<Run*> R;
	//Read Lumi from file
	FILE *fr=fopen(configRun,"r");
	if(fr==NULL) {fprintf(stderr,"No Such file\n"); return 5;}
			{
			char buffer[1024];
			while ( fgets (buffer , 1024 , fr) != NULL )
				{
				if (buffer[0]=='#') continue; 
				long min,max;
				float lumi;
				char pileup1[1024],pileup2[1024],pileup3[1024];
				fprintf(stderr,"%s",buffer);
				sscanf(buffer,"%ld %ld %f %s %s %s",&min,&max,&lumi,pileup1,pileup2,pileup3);
				Run *A=new Run();

				A->min=min;
				A->max=max;
				A->lumi=lumi;
				A->CrossSection=CrossSection;

				//DEBUG
				fprintf(stderr,"minrun=%ld, max=%ld lumi=%f xSec=%f\n",min,max,lumi,CrossSection);
				fprintf(stderr,"going to open %s\n",pileup1); 	

				TH1D*h;

				TFile *f_pu0=TFile::Open(pileup1);
					h=(TH1D*)f_pu0->Get("pileup");
					if(h==NULL) fprintf(stderr,"ERROR: no pu histo\n");
					f->cd();
					A->PU=(TH1D*)h->Clone(Form("pileup_%d",int(R.size())));
					f_pu0->Close();
				fprintf(stderr,"going to open %s\n",pileup2); 	//DEBUG
				TFile *f_pu1=TFile::Open(pileup2);
					h=(TH1D*)f_pu1->Get("pileup");
					if(h==NULL) fprintf(stderr,"ERROR: no pu histo\n");
					f->cd();
					A->PU_DN=(TH1D*)h->Clone(Form("pileup_dn_%d",int(R.size())));
					f_pu1->Close();
				fprintf(stderr,"going to open %s\n",pileup3); 	 //DEBUG
				TFile *f_pu2=TFile::Open(pileup3);
					h=(TH1D*)f_pu2->Get("pileup");
					if(h==NULL) fprintf(stderr,"ERROR: no pu histo\n");
					f->cd();
					A->PU_UP=(TH1D*)h->Clone(Form("pileup_up_%d",int(R.size())));
					f_pu2->Close();
				//Get the Number of event Processed per Run
				TH1F* dummy=new TH1F("dummy","dummy",1,-.5,.5);
				processed->Draw("0>>dummy",Form("mcWeight*(%ld<runNum && runNum<%ld)",min,max),"goff");
				A->nEventsProcessed= dummy->GetBinContent(1);
				delete dummy;
				//DEBUG
				fprintf(stderr," nEv=%f",A->nEventsProcessed);
				if(A->PU ==NULL) fprintf(stderr,"ERROR PU\n");
				if(A->PU_UP ==NULL) fprintf(stderr,"ERROR PU UP\n");
				if(A->PU_DN ==NULL) fprintf(stderr,"ERROR PU DN\n");
				//TODO - update ProcessedData with puINTTrue
				TH1D* PU_A=new TH1D("PU_actual","PU_actual",100,0,100);
				t->Draw("puTrueINT>>PU_actual",Form("mcWeight*(%ld<runNum && runNum<%ld)",min,max ), "goff");
				gDirectory->cd();
				A->PU_actual=(TH1D*) PU_A->Clone(Form("pu_actual_%d",int(R.size())));	
				if(A->PU_actual ==NULL) fprintf(stderr,"ERROR PU actual\n");
				delete PU_A;
				fprintf(stderr," pu Integrals: \n");
				fprintf(stderr," pu  =%lf\n",  A->PU   ->Integral() );
				fprintf(stderr," puUP=%lf\n",  A->PU_UP->Integral() );
				fprintf(stderr," puDN=%lf\n",  A->PU_DN->Integral() );
				fprintf(stderr," puA =%lf\n",  A->PU_actual->Integral() );
				//Normalize PU Histos
				A->Normalize();
				R.push_back(A);
				}
			}

	fprintf(stderr,"creating branches\n");	
	//Creating an empty branch in the tree
	double RDWeight,RDWeightSysUp,RDWeightSysDown,RDWeightBare;
	TBranch *b0=t->Branch("RDWeightBare",&RDWeightBare,"RDWeightBare/D");
	TBranch *b1=t->Branch("RDWeight",&RDWeight,"RDWeight/D");
	TBranch *b2=t->Branch("RDWeightSysUp"  ,&RDWeightSysUp  ,"RDWeightSysUp/D");
	TBranch *b3=t->Branch("RDWeightSysDown",&RDWeightSysDown,"RDWeightSysDown/D");
		
	fprintf(stderr,"Set Branch Address\n");
	//setting the branch for the mcWeight
	float mcWeight=1.0;
	t->SetBranchAddress("mcWeight",&mcWeight);

	int puTrueINT;
        t->SetBranchAddress("puTrueINT",&puTrueINT); 
	
	Int_t runNum;
        t->SetBranchAddress("runNum",&runNum); 

	//Getting the Number of entries in the tree
	long long int NumberEntries=t->GetEntries();
	//looping on the entries in order to add the correct number of entries in the branch
	fprintf(stderr,"----------- START LOOP -----------\n");
	for(long long int i=0;i<NumberEntries;i++){
		t->GetEntry(i);
		if( (i&1048575)==1)fprintf(stderr,"entry %lld/%lld\n",i,NumberEntries);
		//find Run
		int iRun=0;bool runFound=false;
		for(iRun=0;iRun<R.size();iRun++) 
			{
			if( R[iRun]->min <= runNum && runNum<=R[iRun]->max){runFound=true; break;}
			}
		if(!runFound){
				fprintf(stderr,"ERROR RUN %d not found\n",runNum);
				RDWeight=0;RDWeightSysUp=0;RDWeightSysDown=0;RDWeightBare=0;
				}
		else	{
				RDWeight=R[iRun]->weight() *R[iRun]->puReweight(puTrueINT) *kFactor;
				RDWeightSysUp=R[iRun]->weight() *R[iRun]->puReweight(puTrueINT,1) *kFactor;
				RDWeightSysDown=R[iRun]->weight() *R[iRun]->puReweight(puTrueINT,-1) *kFactor;
				RDWeightBare=R[iRun]->weight()*kFactor;
			}
		
	
		b0->Fill();
		b1->Fill();
		b2->Fill();
		b3->Fill();
		}
	//Write the Tree (With OverWrite Option)
	f->cd(Directory);
	t->Write("",TObject::kOverwrite);
	//Close the file
	f->Close();
	//Print a message on stdout
	printf("Done\n");
	return 0;
}

