#ifndef HISTOBOOKER_H
#define HISTOBOOKER_H

#include <TVector3.h> 
#include <TH1.h>
#include <ParallelFiller.h>
#include "filesaver.h"
#include "DBarReader.h"

using namespace std;

class SingleHisto{

	protected:
	std::string collectionname;
	std::vector<TH1*> Histos;
	Binning bins;
	std::string cut;

	public:
	SingleHisto(std::string CollectionName, int nbinsx, float xmin, float xmax,std::string Cut);
	virtual void FillEventByEventScatter(Variables * vars, float (*var) (Variables * vars),float (*secondvar) (Variables * vars),float (*discr_var) (Variables * vars)){};
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
	void Save(FileSaver finalhisto);

};

class SingleScatter: public SingleHisto{
	public:
	SingleScatter(std::string CollectionName, int nbinsx, float xmin, float xmax,int nbinsy, float ymin, float ymax,std::string Cut);
	virtual void FillEventByEventScatter(Variables * vars, float (*var) (Variables * vars),float (*secondvar) (Variables * vars),float (*discr_var) (Variables * vars));
};

class HistoBinCollection: public SingleHisto {

	public:
	HistoBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax,std::string Cut);
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
};

class ScatterBinCollection: public SingleHisto {

        public:
        ScatterBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax,int nbinsy, float ymin, float ymax,std::string Cut);
	virtual void FillEventByEventScatter(Variables * vars, float (*var) (Variables * vars),float (*secondvar) (Variables * vars),float (*discr_var) (Variables * vars));
};

class BinnedHisto: public SingleHisto{
	public:
	BinnedHisto(std::string Name, Binning Bins,std::string Cut);	
        virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)); 
};

class BinnedScatter: public SingleHisto{
	public:
	BinnedScatter(std::string Name, Binning Bins,std::string Cut);	
        virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)); 
};

class HistoBooker{
	
	typedef float (*GetFillinVariable)  (Variables * vars);
    	typedef float (*GetDiscrimVariable) (Variables * vars);

	private:
	std::vector<SingleHisto*> Histos;
	ParallelFiller<SingleHisto*> Filler;

	public:
	HistoBooker();

	void BookSingleHisto(std::string CollectionName, int nbinsx, float xmin, float xmax,std::string Cut, GetFillinVariable var=0x0, GetDiscrimVariable discr_var=0x0);	
	void BookBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax,std::string Cut,GetFillinVariable var=0x0, GetDiscrimVariable discr_var=0x0);
	void BookSingleScatter(std::string CollectionName, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, std::string Cut, GetFillinVariable var=0x0,GetFillinVariable secondvar=0x0, GetDiscrimVariable discr_var=0x0);
	void BookScatterBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax,int nbinsy, float ymin, float ymax,std::string Cut,GetFillinVariable var=0x0,GetFillinVariable secondvar=0x0, GetDiscrimVariable discr_var=0x0);
	void BookBinnedHisto(std::string CollectionName, Binning Bins,std::string Cut, GetFillinVariable var=0x0, GetDiscrimVariable discr_var=0x0);
        void BookBinnedScatter(std::string CollectionName, Binning Bins,std::string Cut, GetFillinVariable var=0x0, GetDiscrimVariable discr_var=0x0);
	

	void FillEverything(DBarReader reader,Variables * vars);
	void SaveEverything(FileSaver finalhisto);	
};

#endif
