
class Distributions {
public:
    virtual void Fill(const TupleVars & tup) = 0;
    virtual void Write(TDirectory * dir) = 0;

    static Distributions * TriggerMC(AnalysisBins b, const Calibrations & calibrations);
    static Distributions *    QualMC(AnalysisBins b, const Calibrations & calibrations);
};


