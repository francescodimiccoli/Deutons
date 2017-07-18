class FileSaver {
   public:
      FileSaver(bool Isfinal=false) : filename("") {fArr=new TObjArray(); IsFinal = Isfinal;}
      FileSaver (string fname,bool Isfinal=false) : filename (fname) {fArr=new TObjArray(); IsFinal = Isfinal;}
      void writeObjsInFolder(string folder, bool recreate = false);
      void writeObjs();
      void setName(string fname) {filename=fname;}
      void Add(TObject* obj) {fArr->Add(obj);}
      string setName() {return filename;}
      bool CheckFile();	      
      TObject * Get(std::string name);	
      TFile * GetFile();
	
   private:
      string filename;
      TObjArray* fArr;
      bool IsFinal;
};

void FileSaver::writeObjsInFolder(string folder, bool recreate)
{
	cout<<"*** Updating "<<filename.c_str()<<" file in "<< folder << endl;
	TFile * fileFinalPlots;

	if(recreate) fileFinalPlots=TFile::Open(filename.c_str(), "RECREATE");
	else fileFinalPlots=TFile::Open(filename.c_str(), "UPDATE");

	if (!fileFinalPlots->GetDirectory(folder.c_str()))
		fileFinalPlots->mkdir(folder.c_str());
	fileFinalPlots->cd   (folder.c_str());
	for (int i = 0; i <= fArr->GetLast(); i++){
		fArr->At(i)->Write(fArr->At(i)->GetName(),2);
		cout<<"Object saved ..."<<endl;
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
