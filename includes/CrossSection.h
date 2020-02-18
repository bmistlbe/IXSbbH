

#ifndef CrossSection_h
#define CrossSection_h



#include "GlobalDefs.h"
#include "PDF.h"
#include "Luminosity.h"
#include "cuba.h"

#include <complex>



//void Cuhre(int NDIM,int  NCOMP,integrand_t CrossSectionIntegrand,void * xs,int NVEC,double EPSREL,double EPSABS,int VERBOSE ,int MINEVAL,int MAXEVAL,int KEY,const char* cc,const char * cc2,int * nregions,int * neval,int * fail,double * integral,double * error,double * prob);

vector<int> lst(int a1,int a2,int a3);

extern int N3LOExpansion;
class CrossSection{
private:
public:
    
    int channel;
    double L;
    double ar;
    double tau;
    double Q;
    Luminosity * Lumi;
    double Lumi1,Lumi1Err;
    int dist;
    int FONNLL;
    double Lmb;
    
    vector<vector<vector<double> > > xs,error,zero,xs_SV,error_SV;
    vector<vector<int> > pos,SVpos;
    vector<vector<vector<double> > > DistributionCoefs;
    
    
    double ZExp[4][8][4][300][6];
    double ZbExp[4][8][4][300][6];
    double WExp[4][8][4][300][6];
    
    double MCPrecision;
    int MCVerbose;
    
    double sqrtz;
    double zval;

    void SetDistributionCoefs();
    void SetZCoefs();
    void SetWCoefs();
    void SetZbCoefs();
    
    CrossSection()
    {
        MCPrecision=1e-3;
        FONNLL=0;
        Lmb=log(4.58*4.58/125.09/125.09);
        
        L=0;
        ar=0;
        tau=0;
        Q=0;
        
        Lumi1=0;
        Lumi1Err=0;
        
        DistributionCoefs=vector<vector<vector<double> > >  (4,vector<vector<double> > (8,vector<double>(4,0)));
        SetDistributionCoefs();
        
        
        //Order , Initial Stat , Log Power, x Power, x Log Power
        for(int i=0;i<4;i++)
            for(int j=0;j<8;j++)
                for(int k=0;k<4;k++)
                    for(int l=0;l<300;l++)
                        for(int m=0;m<6;m++)
                        {
                            ZExp[i][j][k][l][m]=0;
                            ZbExp[i][j][k][l][m]=0;
                            WExp[i][j][k][l][m]=0;
                        }

        
        //Order , Initial Stat , Log Power
        zero=vector<vector<vector<double> > >  (4,vector<vector<double> > (8,vector<double>(4,0)));
        xs=zero;
        error=zero;
        
        SetDistributionCoefs();
        SetZCoefs();
        SetWCoefs();
        SetZbCoefs();
        
        //Soft - Virtual Positions
        //order , initial state , Log Power
        SVpos.push_back(lst(0,0,0));
        
        SVpos.push_back(lst(1,0,0));
        SVpos.push_back(lst(1,0,1));
        
        SVpos.push_back(lst(2,0,0));
        SVpos.push_back(lst(2,0,1));
        SVpos.push_back(lst(2,0,2));
        
        SVpos.push_back(lst(3,0,0));
        SVpos.push_back(lst(3,0,1));
        SVpos.push_back(lst(3,0,2));
        SVpos.push_back(lst(3,0,3));
        
        
        //order , initial state , Log Power
        //NLO
        pos.push_back(lst(1,0,0));
        pos.push_back(lst(1,0,1));
        pos.push_back(lst(1,1,0));
        pos.push_back(lst(1,1,1));
        //NNLO - b ab
        pos.push_back(lst(2,0,0));
        pos.push_back(lst(2,0,1));
        pos.push_back(lst(2,0,2));
        //NNLO - b g
        pos.push_back(lst(2,1,0));
        pos.push_back(lst(2,1,1));
        pos.push_back(lst(2,1,2));
        //NNLO - b q
        pos.push_back(lst(2,2,0));
        pos.push_back(lst(2,2,1));
        pos.push_back(lst(2,2,2));
        //NNLO - b qbar
        pos.push_back(lst(2,3,0));
        pos.push_back(lst(2,3,1));
        pos.push_back(lst(2,3,2));
        //NNLO - b b
        pos.push_back(lst(2,4,0));
        pos.push_back(lst(2,4,1));
        pos.push_back(lst(2,4,2));
        //*/
        //NNLO - g g
        pos.push_back(lst(2,5,0));
        pos.push_back(lst(2,5,1));
        pos.push_back(lst(2,5,2));
        //*/
        //NNLO - q qbar
        pos.push_back(lst(2,6,0));
        pos.push_back(lst(2,6,1));
        pos.push_back(lst(2,6,2));
        //*/
        //N3LO - b ab
        pos.push_back(lst(3,0,0));
        pos.push_back(lst(3,0,1));
        pos.push_back(lst(3,0,2));
        pos.push_back(lst(3,0,3));
        //N3LO - b g
        pos.push_back(lst(3,1,0));
        pos.push_back(lst(3,1,1));
        pos.push_back(lst(3,1,2));
        pos.push_back(lst(3,1,3));
        //N3LO - b q
        pos.push_back(lst(3,2,0));
        pos.push_back(lst(3,2,1));
        pos.push_back(lst(3,2,2));
        pos.push_back(lst(3,2,3));
        //N3LO - b qbar
        pos.push_back(lst(3,3,0));
        pos.push_back(lst(3,3,1));
        pos.push_back(lst(3,3,2));
        pos.push_back(lst(3,3,3));
        //N3LO - b b
        pos.push_back(lst(3,4,0));
        pos.push_back(lst(3,4,1));
        pos.push_back(lst(3,4,2));
        pos.push_back(lst(3,4,3));
        //N3LO - g g
        pos.push_back(lst(3,5,0));
        pos.push_back(lst(3,5,1));
        pos.push_back(lst(3,5,2));
        pos.push_back(lst(3,5,3));
        //N3LO - q qbar
        pos.push_back(lst(3,6,0));
        pos.push_back(lst(3,6,1));
        pos.push_back(lst(3,6,2));
        pos.push_back(lst(3,6,3));
        //N3LO - q g
        pos.push_back(lst(3,7,0));
        pos.push_back(lst(3,7,1));
        pos.push_back(lst(3,7,2));
        pos.push_back(lst(3,7,3));
        //*/
    };
    
    void SetZb(vector<vector<vector<double> > > & values,const vector<double> & zbPows,const vector<double> & zbLogPows);
    void SetZbN3LO(vector<vector<vector<double> > > & values,const vector<double> & zbPows,const vector<double> & zbLogPows);
    void SetW(vector<vector<vector<double> > > & values,const vector<double> & wPows);
    void SetZ(vector<vector<vector<double> > > & values,const vector<double> & sqrtzPows,const vector<double> & zLogPows);
    
    void Integrate();
    void IntegrateSV();
    void ComputeDummyVariables(double z,vector<vector<vector<double> > > & values);
    void ComputeCoefficients(double z);
    double ComputeTotalXS(const vector<vector<vector<double > > > & vec);
    void Evaluate(double z,double x,vector<vector<vector<double> > > & values);
    vector<double> IntegrateLuminosity();
    void SetFONNLL(int in){FONNLL=in;return;};

    
};
static int ChannelIntegrand(const int * ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
static int DistributionIntegrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
static int Integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
static int LuminosityIntegrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);

#endif
