class HistogramBuilder {
protected:
    const TupleVars & tup;
    const CutsManager & cuts;
    virtual void dispatch() = 0;
public:
    void ProcessEvent(const TupleVars & t, const CutsManager & c){
        tup = t;
        cuts = c;
        dispatch();
    }
};

class Proton : public HistogramBuilder {
protected:
    virtual void dispatch(){
        if( tup.GetMCType() != TupleVars::MCTYPE::Proton) return;
        
    }
};
