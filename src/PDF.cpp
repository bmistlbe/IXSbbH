#include "../includes/PDF.h"

vector<double > beta(5,0);
vector<double > gammam(5,0);

vector<double> GridPoints;
double InterP=12;



double LagrangePolynomial(double x,int j,const vector<double> & xval)
{
    double fac=1;
    for(int i=0;i<xval.size();++i)
    {
        if(i==j)
            continue;
        fac*=(x-xval[i])/(xval[j]-xval[i]);
    }
    return fac;
}

int FindGridPoint(double x)
{
    int pt=GridPoints.size()-1;
    for(int i=0;i<GridPoints.size();++i)
        if(GridPoints[i]>x)
        {
            pt=i-1;
            break;
        }
    if(pt<InterP/2)
        pt=InterP/2;
    if(pt>GridPoints.size()-InterP/2-1)
        pt=GridPoints.size()-InterP/2-1;
    return pt;
}


double Interpolator(double xx ,const vector<double> & x, const vector<double> & v)
{
    double res=0;
    for(int i=0;i<x.size();++i)
        res+=LagrangePolynomial(xx,i,x)*v[i];
    return res;
}

double PDFF::GetX(int i,double a,double b,int n)
{
    return log((exp(a*b)*n)/(i - exp(a*b)*i + exp(a*b)*n))/a;
}

double PDFF::ComputeValue(double x,int flav)
{
    int pt=FindGridPoint(x);
    vector<double> fs(InterP,0);
    vector<double> xs(InterP,0);
    for(int i=0;i<InterP;++i)
    {
        fs[i]=GetValue(GridPoints[pt-InterP/2+i],gridmuf,flav);
        xs[i]=log(GridPoints[pt-InterP/2+i]);
    }
    return Interpolator(log(x),xs,fs);
}

void PDFF::CreateGrid(double muf)
{
    gridmode=false;
    gridmuf=muf;
    
    grid_small.assign(11,vector<double>(n,0));
    grid_mid.assign(11,vector<double>(n,0));
    grid_large.assign(11,vector<double>(n,0));
    
    double zz=0;
    for(int kk=-5;kk<=5;kk++)
    {
        if(kk==0)
            GridPoints=pdf->subgrid(muf*muf).get_pid(21).xs();
        else
            GridPoints=pdf->subgrid(muf*muf).get_pid(kk).xs();
        
        int i;
        for(i=0;i<n;++i)
        {
            zz=GetX(i,a1,b1,n);
            if(i!=0)
                grid_small[5+kk][i]=this->ComputeValue(zz,kk);
            
            zz=GetX(i,a2,b2,n);
            if(i!=0)
                grid_mid[5+kk][i]=this->ComputeValue(zz,kk);
            
            zz=GetX(i,a3,b3,n);
            if(i!=0)
                grid_large[5+kk][i]=this->ComputeValue(zz,kk);
        }
    }
    gridmode=true;
    return;
}


double PDFF::ReadFromGrid(double z,int flav)
{
    int bin;
    double a,b;
    vector<double> *grid;
    if(z>b2)
    {
        grid=&grid_large[5+flav];
        a=a3;
        b=b3;
    }
    else if(z>b1)
    {
        grid=&grid_mid[5+flav];
        a=a2;
        b=b2;
    }
    else
    {
        grid=&grid_small[5+flav];
        a=a1;
        b=b1;
    }
    
    bin=(1.0-exp(-a*z))/(1.0-exp(-a*b))*n;
    
    if(bin<InterpolationPower/2)
        bin=InterpolationPower/2;
    if(bin>n-InterpolationPower/2-2)
        bin=n-InterpolationPower/2-2;
    
    double res=0;
    vector<double> ppos(InterpolationPower,0);
    for(int i=0;i<InterpolationPower;++i)
        ppos[i]=GetX(bin-InterpolationPower/2+i,a,b,n);
    for(int i=0;i<InterpolationPower;++i)
        res+=LagrangePolynomial(z,i,ppos)*(*grid)[bin-InterpolationPower/2+i];
    return res;
    
}



void PDFF::InitiatePDF(string sname,int stnr)
{
    setname=sname;
    member=stnr;
    if(stnr==0)
        pdf=static_cast<LHAPDF::GridPDF*>( LHAPDF::mkPDF(setname.c_str()));
    else
        pdf= static_cast<LHAPDF::GridPDF*>(LHAPDF::mkPDF(setname.c_str(),stnr));
    
    return;
}

double PDFF::MyxfxQ2(double x,double Q2,int flavour)
{
    if(x<=0||x>=1)
        return 0;
    return pdf->xfxQ2(flavour, x,Q2);
}

double PDFF::GetValue(double x, double Q, int flavour)
{
    if(gridmode)
        return ReadFromGrid(x,flavour);
    else
        return MyxfxQ2(x,Q*Q,flavour);
}


double PDFF::GetAlphaFromPDF(double mu)
{
    return pdf->alphasQ2(mu*mu);
}

double PDFF::GetAlpha(double mu)
{
    int it=10000;
    double deltamu=(mu-mZ)/(double)it;
    double a=alphamZ;
    
    for(int i=0;i<it;++i)
    {
        double locmu=mZ+deltamu*i;
        a+=-a*a/Pi*(beta[0]+beta[1]*a/Pi+beta[2]*a/Pi*a/Pi+beta[3]*a*a*a/Pi/Pi/Pi)*(locmu*locmu-(locmu-deltamu)*(locmu-deltamu))/locmu/locmu;
    }
    return a;
}

double PDFF::GetAlpha(double mu,int o)
{
    int it=10000;
    double deltamu=(mu-mZ)/(double)it;
    double a=alphamZ;
    
    vector<double> bb=beta;
    if(o<3)
        bb[3]=0;
    if(o<2)
        bb[2]=0;
    if(o<1)
        bb[1]=0;
    
    for(int i=0;i<it;++i)
    {
        double locmu=mZ+deltamu*i;
        a+=-a*a/Pi*(bb[0]+bb[1]*a/Pi+bb[2]*a/Pi*a/Pi+bb[3]*a*a*a/Pi/Pi/Pi)*(locmu*locmu-(locmu-deltamu)*(locmu-deltamu))/locmu/locmu;
    }
    return a;
}


double PDFF::GetMT(double mu,int o)
{
    int it=10000;
    int nf=5;
    double locm=mT;
    double deltamu=(mu-mT)/(double)it;
    
    double gamma0=gammam[0];
    double gamma1=gammam[1];
    double gamma2=gammam[2];
    double gamma3=gammam[3];
    double gamma4=gammam[4];
    
    
    if(o<3)
        gamma3=0;
    if(o<2)
        gamma2=0;
    if(o<1)
        gamma1=0;
    
    for(int i=0;i<=it;++i)
    {
        double locmu=mT+deltamu*i;
        double as=GetAlpha(locmu,o)/Pi;
        double gamma=as*gamma0
            +as*as*gamma1
            +as*as*as*gamma2
            +as*as*as*as*gamma3
            +as*as*as*as*as*gamma4;
        locm+=-locm*gamma*(locmu*locmu-(locmu-deltamu)*(locmu-deltamu))/locmu/locmu;
    }
    return locm;
}



double PDFF::GetMT2(double mu,int o)
{
    int it=100000;
    int nf=5;
    double locm=mT;
    double deltamu=(mu-mT)/(double)it;
    
    double gamma0=gammam[0];
    double gamma1=gammam[1];
    double gamma2=gammam[2];
    double gamma3=gammam[3];
    double gamma4=gammam[4];
    
    vector<double> bb=beta;
    if(o<3)
        bb[3]=0;
    if(o<2)
        bb[2]=0;
    if(o<1)
        bb[1]=0;
    
    if(o<3)
        gamma3=0;
    if(o<2)
        gamma2=0;
    if(o<1)
        gamma1=0;
    
    double a=GetAlpha(mT,o);
    for(int i=0;i<=it;++i)
    {
        double locmu=mT+deltamu*i;
        a+=-a*a/Pi*(bb[0]+bb[1]*a/Pi+bb[2]*a/Pi*a/Pi+bb[3]*a*a*a/Pi/Pi/Pi)*(locmu*locmu-(locmu-deltamu)*(locmu-deltamu))/locmu/locmu;
        double as=a/Pi;
        double gamma=as*gamma0
        +as*as*gamma1
        +as*as*as*gamma2
        +as*as*as*as*gamma3
        +as*as*as*as*as*gamma4;
        locm+=-locm*gamma*(locmu*locmu-(locmu-deltamu)*(locmu-deltamu))/locmu/locmu;
    }
    return locm;
}
