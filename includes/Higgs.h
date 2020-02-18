

#ifndef Higgs_h
#define Higgs_h



#include "GlobalDefs.h"
#include "Luminosity.h"
#include "CrossSection.h"


class Higgs{
private:
public:
    
    int PDFMember;
    string PDFSet;
    
    double mh;
    double muf;
    double mur;
    double E;
    double mZ;
    double mb0,mb1,mb2,mb3;
    double mb_input;
    double vev;
    double conversion;
    double ar;
    double Born;
    double Mb_OS;
    bool perturbative_mass_evolution;
    
    double MCPrecision;
    int MCVerbose;
    
    Luminosity Lumi;
    CrossSection xs;
    
    vector<vector<vector<double> > > Integrals,Errors;
    
    bool integrated;
    
    Higgs()
    {
        perturbative_mass_evolution=false;
        integrated=false;
        
        Mb_OS=4.58;
        ar=0.118/Pi;
        mh=125;
        muf=mh;
        mur=mh;
        E=13000;
        mZ=91.1876;
        mb_input=4.18;
        mb0=mb1=mb2=mb3=mb_input;
        
        conversion=3.893793656e8;
        vev=246.221;
        Born=conversion;
        
        xs.Lumi=&Lumi;
        
        MCPrecision=1e-3;
        MCVerbose=1;
    };
    
    //Set Stuff
    void SetMh(double Q){mh=Q;integrated=false;SetMBOS(Mb_OS);return;};
    void SetMuf(double Q){muf=Q;integrated=false;return;};
    void SetE(double Q){E=Q;integrated=false;return;};
    void SetmZ(double Q){mZ=Q;return;};
    void SetMur(double Q);
    void Setmb(double Q);
    void SetPDF(int mem,string set){PDFMember=mem;PDFSet=set;integrated=false;return;};
    void SetVerbose(int v){MCVerbose=v;return;};
    void SetPrecision(double prec){MCPrecision=prec;integrated=false;return;};
    void SetFONNLL(int in){xs.SetFONNLL(in);integrated=false;return;};
    void SetMBOS(double mbos){Mb_OS=mbos;xs.Lmb=log(mbos*mbos/mh/mh);integrated=false;return;};
    
    
    //Do Stuff
    void IntegrateCrossSection();
    vector<double> IntegrateLuminosity(double zz);
    void ApplyPrefactors(vector<vector<double> > & result,vector<vector<double> > & error);
    
    


};


vector<double> MuREvolution(double muf,double mur,const vector<double> & xs);
vector<double> MassEvolution(double muf,double mur,const vector<double> & xs);

#endif
