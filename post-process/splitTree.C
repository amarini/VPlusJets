
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"

#include <vector>
#include <set>
#include <string>

using namespace std;

// --- see ${ROOTSYS}/tutorials/tree/copytree.C 

//void deleteWeightBranches(string oldName,string newName)
void splitTree(string oldName,string newName)
{ 

   //Get old file, old tree and set top branch address
   TFile *oldfile = new TFile(oldName.c_str());
   TTree *oldtree = (TTree*)oldfile->Get("accepted/events");
   oldtree->SetBranchStatus("*",1);
   //oldtree->SetBranchStatus("eventWeight",0);
   //oldtree->SetBranchStatus("PUWeight*",0);
   //oldtree->SetBranchStatus("RDWeight*",0);

   set<string> dontCopyAuto;
	dontCopyAuto.insert("events");

   Long64_t nentries=30000;
   Long64_t allEntries=oldtree->GetEntries();

   int nTrees=(allEntries%nentries == 0) ? allEntries/nentries :allEntries/nentries +1;

   for( long iTree = 0;iTree< nTrees ;iTree++){

   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile( Form("%s_%ld.root",newName.c_str(),iTree),"recreate");
   //created accepted dir
   TDirectory *newdir=newfile->mkdir("accepted");
   newfile->cd("accepted");
   //TTree::SetMaxTreeSize(size);
  // TTree *newtree = oldtree->CloneTree();
   	//TTree* CopyTree(const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0)

   TTree *newtree = oldtree->CopyTree("","",nentries,iTree * nentries);
   newtree->Write("",TObject::kOverwrite);
	
	//copy everything else in the accepted directory
	TDirectory *olddir=(TDirectory*)oldfile->Get("accepted");
	TList *l=olddir->GetListOfKeys();
   	newfile->cd("accepted");
	for(int i=0;i<l->GetSize();i++)
		{
		printf("considering item %s\n",l->At(i)->GetName());
		string name=l->At(i)->GetName();
		TKey *key= (TKey*)l->At(i);
		if ( dontCopyAuto.find( name ) != dontCopyAuto.end() ) {
			printf(" --> dont copy auto\n");
			continue;
			}
		newdir->cd();
		printf(" --> going to copy\n");
		olddir->cd();
		TObject* obj  =  key->ReadObj();// ->Clone();	
		if ( obj->InheritsFrom("TTree") ) 
			{
			TTree *t = (TTree*)oldfile->Get( Form("%s/%s",newdir->GetName(),name.c_str()) );	
			if( t== NULL ) printf("ERROR t is NULL\n");
			t->SetBranchStatus("*",1);
			newdir->cd();
			printf("Going to copy in %s\n",gDirectory->GetName());
			TTree *T=t->CloneTree();
			T->Write("",TObject::kOverwrite);
			}
		else if ( obj->InheritsFrom("TH1") || obj->InheritsFrom("TH2")  ) 
			{
         		newdir->cd();
         		obj->Write();
         		delete obj;
			}
		}
   newtree->Write("",TObject::kOverwrite);

   newdir->Write("",TObject::kOverwrite);
   newdir-> SaveSelf(kTRUE);
   newfile->Write("",TObject::kOverwrite);
   delete newfile;
   }
   delete oldfile;
}
