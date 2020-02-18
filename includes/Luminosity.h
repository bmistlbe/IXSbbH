


#ifndef Lumi_h
#define Lumi_h


#include "GlobalDefs.h"
#include "PDF.h"

class Luminosity{


private:
public:
    
    
    PDFF pdf;
    double L[8];
    double L0;
    
    Luminosity()
    {
        L0=0;
        for(int i=0;i<8;++i)
            L[i]=0;
    };
    
    
    void SetLuminosity(double z, double x,double tau,double Q);
    void SetBottomLuminosity(double z, double x,double tau,double Q);
    double Get(int channel){return L[channel];};
    
    double Lumi(double z, double x, double tau,double Q,int parton1,int parton2);
    double Lbab(double z, double x, double tau,double Q);
    double Lbg(double z, double x, double tau,double Q);
    double Lbq(double z, double x, double tau,double Q);
    double Lbqb(double z, double x, double tau,double Q);
    double Lbb(double z, double x, double tau,double Q);
    double Lgg(double z, double x, double tau,double Q);
    double Lqqbar(double z, double x, double tau,double Q);
    double Lqg(double z, double x, double tau,double Q);
    
};



#endif
