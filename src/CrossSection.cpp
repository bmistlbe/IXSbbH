#include "../includes/CrossSection.h"


int N3LOExpansion=-1;


vector<int> lst(int a1,int a2,int a3)
{
    vector<int>loc(3,0);
    loc[0]=a1;
    loc[1]=a2;
    loc[2]=a3;
    return loc;
}

void CrossSection::SetDistributionCoefs()
{
#include "../XS/Distributions.txt"
    return;
}

void CrossSection::ComputeDummyVariables(double z,vector<vector<vector<double> > > & values)
{
    sqrtz=sqrt(z);
    values=zero;
    int i=0;
    
    if(z>=0.75)
    {
        vector<double> zbPows(50,0);
        vector<double> zbLogPows(6,0);
        
#pragma omp parallel for private(i) schedule(dynamic)
        for(i=0;i<50;i++)
            zbPows[i]=pow(1.0-z,i);
        zbLogPows[0]=1;
        zbLogPows[1]=log(1-z);
        zbLogPows[2]=pow(log(1-z),2);
        zbLogPows[3]=pow(log(1-z),3);
        zbLogPows[4]=pow(log(1-z),4);
        zbLogPows[5]=pow(log(1-z),5);
        SetZb(values,zbPows,zbLogPows);
    }
    else if (z<0.75&&z>=0.0769231)
    {
        vector<double> wPows(201,0);
        
#pragma omp parallel for private(i) schedule(dynamic)
        for(i=0;i<wPows.size();i++)
            wPows[i]=pow(0.5-z,i);
        SetW(values,wPows);
    }
    else
    {
        vector<double> zPows(101,0);
        vector<double> zLogPows(6,0);
        
#pragma omp parallel for private(i) schedule(dynamic)
        for(i=0;i<zPows.size();i++)
            zPows[i]=pow(z,i-1);
        zLogPows[0]=1;
        zLogPows[1]=log(z);
        zLogPows[2]=pow(log(z),2);
        zLogPows[3]=pow(log(z),3);
        zLogPows[4]=pow(log(z),4);
        zLogPows[5]=pow(log(z),5);
        SetZ(values,zPows,zLogPows);
    }
    //*/
    
    if(N3LOExpansion!=-1)
    {
        vector<double> zbPows(N3LOExpansion,0);
        vector<double> zbLogPows(6,0);
        
#pragma omp parallel for private(i) schedule(dynamic)
        for(i=0;i<zbPows.size();i++)
            zbPows[i]=pow(1.0-z,i);
        zbLogPows[0]=1;
        zbLogPows[1]=log(1-z);
        zbLogPows[2]=pow(log(1-z),2);
        zbLogPows[3]=pow(log(1-z),3);
        zbLogPows[4]=pow(log(1-z),4);
        zbLogPows[5]=pow(log(1-z),5);
        SetZbN3LO(values,zbPows,zbLogPows);
    }
    //*/
    return;
}


void CrossSection::SetZbN3LO(vector<vector<vector<double> > > & values,const vector<double> & zbPows,const vector<double> & zbLogPows)
{
    
    //Order , Initial Stat , Log Power, x Power, x Log Power
    int p;
#pragma omp parallel for private(p) schedule(dynamic)
    for(p=0;p<pos.size();++p)
    {
        if(pos[p][0]<3||pos[p][2]!=0)
            continue;
        values[pos[p][0]][pos[p][1]][pos[p][2]]=0;
        for(int i=0;i<zbPows.size();++i)
        {
            for(int j=0;j<=2*pos[p][0]-1;++j)
            {
                if(ZbExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]==0)
                    continue;
                values[pos[p][0]][pos[p][1]][pos[p][2]]+=ZbExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]*zbPows[i]*zbLogPows[j];
            }
        }
    }
    return;
}


void CrossSection::SetZb(vector<vector<vector<double> > > & values,const vector<double> & zbPows,const vector<double> & zbLogPows)
{
    
    //Order , Initial Stat , Log Power, x Power, x Log Power
    int p;
#pragma omp parallel for private(p) schedule(dynamic)
    for(p=0;p<pos.size();++p)
    {
        for(int i=0;i<zbPows.size();++i)
        {
            for(int j=0;j<=2*pos[p][0]-1;++j)
            {
                if(ZbExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]==0)
                    continue;
                values[pos[p][0]][pos[p][1]][pos[p][2]]+=ZbExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]*zbPows[i]*zbLogPows[j];
            }
        }
    }
    return;
}


void CrossSection::SetZ(vector<vector<vector<double> > > & values,const vector<double> & zPows,const vector<double> & zLogPows)
{
    
    //Order , Initial Stat , Log Power, x Power, x Log Power
    int p;
#pragma omp parallel for private(p) schedule(dynamic)
    for(p=0;p<pos.size();++p)
    {
        for(int i=0;i<zPows.size();++i)
        {
            for(int j=0;j<=2*pos[p][0]-1;++j)
            {
                if(ZExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]==0)
                    continue;
                values[pos[p][0]][pos[p][1]][pos[p][2]]+=ZExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]*zPows[i]*zLogPows[j];
            }
        }
    }
    return;
}



void CrossSection::SetW(vector<vector<vector<double> > > & values,const vector<double> & wPows)
{
    
    //Order , Initial Stat , Log Power, x Power, x Log Power
    int p;
#pragma omp parallel for private(p) schedule(dynamic)
    for(p=0;p<pos.size();++p)
    {
        for(int i=0;i<wPows.size();++i)
        {
            for(int j=0;j<=2*pos[p][0]-1;++j)
            {
                if(WExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]==0)
                    continue;
                values[pos[p][0]][pos[p][1]][pos[p][2]]+=WExp[pos[p][0]][pos[p][1]][pos[p][2]][i][j]*wPows[i];
            }
        }
    }
    return;
}



void CrossSection::IntegrateSV()
{
    int NDIM=2;
    int NCOMP=SVpos.size()+1;
    int NVEC=1;
    double EPSREL= MCPrecision;
    double EPSABS= 1e-8;
    int VERBOSE=MCVerbose;
    int LAST=4;
    int SEED=0;
    int MINEVAL=10000;
    int MAXEVAL=500000000;
    int NSTART=1000;
    int NINCREASE=500;
    int NBATCH=1000;
    int KEY=9;
    double res[NCOMP], err[NCOMP], chi[NCOMP];
    int comp, nregions, neval, fail,cuhre_key=9;
     //*/
    
    Cuhre(NDIM, NCOMP, DistributionIntegrand, this, NVEC,EPSREL, EPSABS, VERBOSE ,MINEVAL, MAXEVAL, KEY, NULL, NULL,&nregions, &neval, &fail, res,err,chi);
    
    xs_SV=zero;
    error_SV=zero;
    for(int i=0;i<SVpos.size();i++)
    {
        xs_SV[SVpos[i][0]][SVpos[i][1]][SVpos[i][2]]=res[i+1];
        error_SV[SVpos[i][0]][SVpos[i][1]][SVpos[i][2]]=err[i+1];
        //cout<<"SV Result: "<<res[i+1] <<" +- "<<err[i+1]<<" +- "<<fabs(err[i+1]/res[i+1])*100<<" %"<<endl;
    }
    //*/
    return;
}


vector<double> CrossSection::IntegrateLuminosity()
{
    
    int NDIM=2;
    int NCOMP=9;
    int NVEC=1;
    double EPSREL= MCPrecision;
    double EPSABS= 1e-8;
    int VERBOSE=MCVerbose;
    int LAST=4;
    int SEED=0;
    int MINEVAL=10000;
    int MAXEVAL=500000000;
    int NSTART=1000;
    int NINCREASE=500;
    int NBATCH=1000;
    int KEY=9;
    double res[NCOMP], err[NCOMP], chi[NCOMP];
    int comp, nregions, neval, fail,cuhre_key=9;
    //*/
    Cuhre(NDIM, NCOMP, LuminosityIntegrand, this, NVEC,EPSREL, EPSABS, VERBOSE ,MINEVAL, MAXEVAL, KEY, NULL, NULL,&nregions, &neval, &fail, res,err,chi);
    
    vector<double> lumis;
    for(int i=0;i<9;i++)
        lumis.push_back(res[i]);
    //*/
    return lumis;
}



void CrossSection::Integrate()
{
    xs_SV=zero;
    error_SV=zero;

    if(FONNLL==0)
        IntegrateSV();
    
    int NDIM=2;
    int NCOMP=pos.size()+1;
    int NVEC=1;
    double EPSREL= MCPrecision;
    double EPSABS= 1e-7;
    int VERBOSE=MCVerbose;
    int LAST=4;
    int SEED=0;
    int MINEVAL=10000;
    int MAXEVAL=500000000;
    int NSTART=1000;
    int NINCREASE=500;
    int NBATCH=1000;
    int KEY=9;
    double res[NCOMP], err[NCOMP], chi[NCOMP];
    int comp, nregions, neval, fail,cuhre_key=9;
    //*/

    xs=zero;
    error=zero;

    Cuhre(NDIM, NCOMP, Integrand, this, NVEC,EPSREL, EPSABS, VERBOSE ,MINEVAL, MAXEVAL, KEY, NULL, NULL,&nregions, &neval, &fail, res,err,chi);
    
    xs=zero;
    error=zero;
    for(int i=0;i<pos.size();++i)
    {
        xs[pos[i][0]][pos[i][1]][pos[i][2]]=res[i+1]-1e-6;
        error[pos[i][0]][pos[i][1]][pos[i][2]]=err[i+1]*err[i+1];
        //cout<<"Integration Result "<<i<<" = "<<res[i+1]<<" +- "<<err[i+1]<<" +- "<<fabs(err[i+1]/res[i+1])*100 <<" %"<<endl;
    }
    //*/

    //Add the distributions
    for(int i=0;i<SVpos.size();++i)
    {
        xs[SVpos[i][0]][SVpos[i][1]][SVpos[i][2]]+=xs_SV[SVpos[i][0]][SVpos[i][1]][SVpos[i][2]];
        error[SVpos[i][0]][SVpos[i][1]][SVpos[i][2]]+=error_SV[SVpos[i][0]][SVpos[i][1]][SVpos[i][2]]*error_SV[SVpos[i][0]][SVpos[i][1]][SVpos[i][2]];
    }
    //*/
    
    //Errors were squared , take the sqrt
    for(int i=0;i<error.size();i++)
        for(int j=0;j<error[i].size();++j)
            for(int k=0;k<error[i][j].size();++k)
                error[i][j][k]=sqrt(error[i][j][k]);
    return;
}



double CrossSection::ComputeTotalXS(const vector<vector<vector<double > > > & vec)
{
    double tot=0;
    for(int i=3;i<vec.size();++i)
        for(int j=0;j<vec[i].size();++j)
            for(int k=0;k<vec[i][j].size();++k)
                tot+=pow(ar,i)*pow(L,k)*vec[i][j][k];
    return tot;
}



void CrossSection::Evaluate(double z,double x,vector<vector<vector<double> > > & values)
{
    values=zero;
    Lumi->SetLuminosity(z,x,tau,Q);
    ComputeDummyVariables(z,values);

    for(int i=0;i<pos.size();++i)
        values[pos[i][0]][pos[i][1]][pos[i][2]]*=Lumi->L[pos[i][1]];
    return;
}

static int ChannelIntegrand(const int * ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
    CrossSection * xs=(CrossSection*) userdata;
    ff[0]=1e-8;
    if(xx[1]<xs->tau/xx[0])
        return 0;
    
    int chan=xs->channel;
    vector<vector<vector<double> > > values;
    xs->Evaluate(xx[0],xx[1],values);
    
    ff[0]+=values[xs->pos[chan][0]][xs->pos[chan][1]][xs->pos[chan][2]];
    return 0;
}


static int Integrand(const int * ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
    CrossSection * xs=(CrossSection*) userdata;
    for(int i=0;i<xs->pos.size()+1;i++)
        ff[i]=1e-6;
    
    if(xx[1]<xs->tau/xx[0])
        return 0;
    
    vector<vector<vector<double> > > values;
    xs->Evaluate(xx[0],xx[1],values);
    for(int i=0;i<xs->pos.size();++i)
        ff[i+1]+=values[xs->pos[i][0]][xs->pos[i][1]][xs->pos[i][2]];
    
    ff[0]=0;
    ff[0]=xs->ComputeTotalXS(values)+1e-6;
    xs=0;
    return 0;
}


static int DistributionIntegrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
    
    CrossSection * xs=(CrossSection*) userdata;
    
    for(int i=0;i<xs->SVpos.size()+1;++i)
        ff[i]=0;
    
    if(xx[1]<xs->tau)
        return 0;
    
    xs->Lumi->SetBottomLuminosity(xx[0],xx[1],xs->tau,xs->Q);
    vector<double> Distributions(7,0);
    
    Distributions[0]=xs->Lumi->L0;
    for(int i=0;i<7;i++)
    {
        double Lz=0;
        if(xx[1]>xs->tau/xx[0])
            Lz=xs->Lumi->L[0];
        Distributions[i+1]=pow(log(1-xx[0]),i)/(1-xx[0])*(Lz-xs->Lumi->L0);
    }
    
    for(int i=0;i<xs->SVpos.size();++i)
        for(int j=0;j<Distributions.size();++j)
            ff[i+1]+=Distributions[j]*xs->DistributionCoefs[xs->SVpos[i][0]][j][xs->SVpos[i][2]];
    ff[0]=ff[xs->SVpos.size()]*pow(xs->L,3)+ff[xs->SVpos.size()-1]*pow(xs->L,2)+ff[xs->SVpos.size()-2]*pow(xs->L,1)+ff[xs->SVpos.size()-3];
    
    xs=0;
    return 0;
}





static int LuminosityIntegrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
    
    CrossSection * xs=(CrossSection*) userdata;
    
    for(int i=0;i<9;++i)
        ff[i]=0;
    
    if(xx[0]<xs->tau)
        return 0;
    xs->Lumi->SetLuminosity(xs->zval,xx[0],xs->tau,xs->Q);

    ff[0]=xs->Lumi->L0;
    if(xx[0]<xs->tau/xs->zval)
        return 0;
    
    for(int i=0;i<8;i++)
        ff[i+1]=xs->Lumi->L[i];
    xs=0;
    return 0;
}




