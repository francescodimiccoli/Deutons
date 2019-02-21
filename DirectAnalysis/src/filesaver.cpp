#include "filesaver.h"


void FileSaver::writeObjsInFolder(string folder, bool recreate)
{
	cout<<"*** Updating "<<filename.c_str()<<" file in "<< folder << endl;
	TFile * fileFinalPlots;

	if(recreate) fileFinalPlots=TFile::Open(filename.c_str(), "RECREATE");
	else fileFinalPlots=TFile::Open(filename.c_str(), "UPDATE");

	if (!fileFinalPlots->GetDirectory(folder.c_str()))
		fileFinalPlots->mkdir(folder.c_str());
	fileFinalPlots->cd   (folder.c_str());
	cout<<"vrce: "<<fArr->At(0)<<endl;
	if(fArr->At(0)){
		cout<<fArr->At(0)->GetName()<<endl; 	
		for (int i = 0; i <= fArr->GetLast(); i++){
			cout<<fArr->At(i)<<" "<<fArr->At(i)->GetName()<<endl;
			fArr->At(i)->Write(fArr->At(i)->GetName(),2);
			cout<<"Object saved ..."<<endl;
		}
	}
	fileFinalPlots->Flush();
	fileFinalPlots->Write();
        fileFinalPlots->Close();
        fArr->Clear();
        return;
}

bool FileSaver::CheckFile(){
	bool check;
	TFile * fileFinalPlots;
	fileFinalPlots=TFile::Open(filename.c_str(), "READ");
	check = (fileFinalPlots>0);
	if(check) fileFinalPlots->Close();
	return check;
}

TObject *  FileSaver::Get(std::string name){
	TObject * Obj;
	TFile * fileFinalPlots=TFile::Open(filename.c_str(), "READ");
	if(!fileFinalPlots) return Obj;
	fileFinalPlots->GetObject(name.c_str(),Obj);
	return Obj;	
}


TFile * FileSaver::GetFile(){
	TFile * fileFinalPlots=TFile::Open(filename.c_str(), "READ");
	return fileFinalPlots;
}


