#include "../includes/Luminosity.h"


void Luminosity::SetBottomLuminosity(double z, double x,double tau,double Q)
{
    L0=Lbab(1,x,tau,Q);
    L[0]=Lbab(z,x,tau,Q);
    return;
}

void Luminosity::SetLuminosity(double z, double x,double tau,double Q)
{
    L0=Lbab(1,x,tau,Q);
    L[0]=Lbab(z,x,tau,Q);
    L[1]=Lbg(z,x,tau,Q);
    L[2]=Lbq(z,x,tau,Q);
    L[3]=Lbqb(z,x,tau,Q);
    L[4]=Lbb(z,x,tau,Q);
    L[5]=Lgg(z,x,tau,Q);
    L[6]=Lqqbar(z,x,tau,Q);
    L[7]=Lqg(z,x,tau,Q);
    return;
}

double Luminosity::Lumi(double z, double x, double tau,double Q,int parton1,int parton2)
{
    return pdf.GetValue(x,Q,parton1)*pdf.GetValue(tau/x/z,Q,parton2)/x;
}


double Luminosity::Lbab(double z, double x, double tau,double Q)
{
    return 2*Lumi(z,x,tau,Q,-5,5);
}

double Luminosity::Lbb(double z, double x, double tau,double Q)
{
    return Lumi(z,x,tau,Q,5,5)+Lumi(z,x,tau,Q,-5,-5);
}

double Luminosity::Lbg(double z, double x, double tau,double Q)
{
    return 2*Lumi(z,x,tau,Q,5,0)+2*Lumi(z,x,tau,Q,-5,0);
}

double Luminosity::Lbq(double z, double x, double tau,double Q)
{
    double lumi=0;
    for(int i=1;i<=4;i++)
        lumi+=2*Lumi(z,x,tau,Q,i,5)+2*Lumi(z,x,tau,Q,-i,-5);
    return lumi;
}

double Luminosity::Lbqb(double z, double x, double tau,double Q)
{
    double lumi=0;
    for(int i=1;i<=4;i++)
        lumi+=2*Lumi(z,x,tau,Q,i,-5)+2*Lumi(z,x,tau,Q,-i,5);
    return lumi;
}

double Luminosity::Lgg(double z, double x, double tau,double Q)
{
    return Lumi(z,x,tau,Q,0,0);
}

double Luminosity::Lqg(double z, double x, double tau,double Q)
{
    double lumi=0;
    for(int i=-4;i<=4;i++)
    {
        if(i==0)
            continue;
        lumi+=2*Lumi(z,x,tau,Q,i,0);
    }
    return lumi;
}

double Luminosity::Lqqbar(double z, double x, double tau,double Q)
{
    double lumi=0;
    for(int i=1;i<=4;i++)
    {
        if(i==0)
            continue;
        lumi+=2*Lumi(z,x,tau,Q,-i,i);
    }
    return lumi;
}
