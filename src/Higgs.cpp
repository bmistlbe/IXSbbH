#include "../includes/Higgs.h"


void Higgs::Setmb(double Q)
{
    mb_input=Q;
    Lumi.pdf.SetMT(mb_input);
    mb0=Lumi.pdf.GetMT2(mur,0);
    mb1=Lumi.pdf.GetMT2(mur,1);
    mb2=Lumi.pdf.GetMT2(mur,2);
    mb3=Lumi.pdf.GetMT2(mur,3);
    return;
    
}

void Higgs::SetMur(double Q)
{
    mur=Q;
    ar=Lumi.pdf.GetAlpha(mur)/Pi;
    mb0=Lumi.pdf.GetMT2(mur,0);
    mb1=Lumi.pdf.GetMT2(mur,1);
    mb2=Lumi.pdf.GetMT2(mur,2);
    mb3=Lumi.pdf.GetMT2(mur,3);
    xs.ar=ar;

    return;
}

vector<double> Higgs::IntegrateLuminosity(double zz)
{
    Lumi.pdf.InitiatePDF(PDFSet,PDFMember);
    Lumi.pdf.SetMZ(mZ);
    SetMur(mur);
    xs.L=2*log(mh/muf);
    xs.tau=mh*mh/E/E;
    xs.Q=muf;
    xs.MCPrecision=MCPrecision;
    xs.MCVerbose=MCVerbose;
    xs.zval=zz;
    return xs.IntegrateLuminosity();
}

void Higgs::IntegrateCrossSection()
{
    Lumi.pdf.InitiatePDF(PDFSet,PDFMember);
    Lumi.pdf.SetMZ(mZ);
    Lumi.pdf.SetMT(mb_input);
    SetMur(mur);
    
    
    xs.L=2*log(mh/muf);
    xs.tau=mh*mh/E/E;
    xs.Q=muf;
    xs.MCPrecision=MCPrecision;
    xs.MCVerbose=MCVerbose;
    xs.Integrate();
    integrated=true;
    return;
}

void Higgs::ApplyPrefactors(vector<vector<double> > & result,vector<vector<double> > & error)
{
    if(!integrated)
        IntegrateCrossSection();

    ar=Lumi.pdf.GetAlpha(mur)/Pi;

    Born=conversion*Pi/6.0/mh/mh/vev/vev;
    vector<double> mbs={pow(mb0,2),pow(mb1,2),pow(mb2,2),pow(mb3,2)};
    double pref=Born;
    double LF=2*log(mh/muf);
    double LR=2*log(mh/mur);
    
    
    result=vector<vector<double> > (8,vector<double>(4,0));
    error=vector<vector<double> > (8,vector<double>(4,0));
    
    cout<<"Value of mb at order 0  mur="<<mur<<" is "<<mb0<<endl;
    cout<<"Value of mb at order 1  mur="<<mur<<" is "<<mb1<<endl;
    cout<<"Value of mb at order 2  mur="<<mur<<" is "<<mb2<<endl;
    cout<<"Value of mb at order 3  mur="<<mur<<" is "<<mb3<<endl;
    cout<<"Value of alpha at mur="<<mur<<" is "<<ar*Pi<<endl;
    cout<<"Value of alpha at mZ="<<mZ<<" is "<<Lumi.pdf.GetAlphaFromPDF(mZ)<<endl<<endl<<endl;
    
    
    vector<string> chans={"b ab","b g","b q","b qbar","b b","g g","q qbar","q g"};
    //Add all Factorisation Logs
    for(int i=0;i<xs.xs.size();++i)
        for(int j=0;j<xs.xs[i].size();++j)
            for(int k=0;k<xs.xs[i][j].size();++k)
            {
                result[j][i]+=pref*xs.xs[i][j][k]*pow(LF,k);
                error[j][i]+=pow(pref*xs.error[i][j][k]*pow(LF,k),2);
            }
    for(int i=0;i<error.size();++i)
        for(int j=0;j<error[i].size();++j)
            error[i][j]=sqrt(error[i][j]);
    
    
    //Display all partonic channels
    /*cout<<endl<<endl;
    cout<<"Before mur evolution: "<<endl;
    for(int i=0;i<result.size();i++)
        for(int j=0;j<result[i].size();j++)
            cout<<"channel "<<chans[i]<<" order "<<j<<" ="<<mbs[j]*pow(ar,j)*result[i][j]<<" +- "<<mbs[j]*pow(ar,j)*error[i][j]<<" +- "<<fabs(error[i][j]/result[i][j])*100<<" %"<<endl;
    //*/
    //Do mb evolution purely perturbatively
    if(perturbative_mass_evolution)
    {
        mbs={mb_input,mb_input,mb_input,mb_input,mb_input};
    }
    
    //Apply MuR evolution 
    for(int i=0;i<result.size();++i)
    {
        if(perturbative_mass_evolution)
        {
            result[i]=MassEvolution(mb_input,muf,result[i]);
            error[i]=MassEvolution(mb_input,muf,error[i]);
        }
        result[i]=MuREvolution(muf,mur,result[i]);
        error[i]=MuREvolution(muf,mur,error[i]);
    }
    
    //Display all partonic channels
    /*cout<<endl<<endl;
    cout<<"After mur evolution: "<<endl;
    for(int i=0;i<result.size();i++)
        for(int j=0;j<result[i].size();j++)
            cout<<"channel "<<chans[i]<<" order "<<j<<" ="<<mbs[j]*pow(ar,j)*result[i][j]<<" +- "<<mbs[j]*pow(ar,j)*error[i][j]<<" pb +- "<<fabs(error[i][j]/result[i][j])*100<<" %"<<endl;
    cout<<endl<<endl;
    //*/
    
    //Add powers of alpha_S
    for(int i=0;i<result.size();i++)
        for(int j=0;j<result[i].size();j++)
        {
            result[i][j]*=pow(ar,j);
            error[i][j]*=pow(ar,j);
        }
    
    
    //Add to lower order terms and sum over channels
    vector<double> xs_per_order(4,0),xs_tot(4,0);
    vector<double> error_per_order(4,0),error_tot(4,0);
    for(int i=0;i<xs_per_order.size();i++)
    {
        for(int j=0;j<result.size();++j)
        {
            xs_per_order[i]+=result[j][i];
            error_per_order[i]+=error[j][i]*error[j][i];
        }
        error_per_order[i]=sqrt(error_per_order[i]);
        cout<<"XS at N^"<<i<<"LO = "<<mbs[i]*xs_per_order[i]<<" +- "<<mbs[i]*error_per_order[i]<<" pb +- "<<fabs(error_per_order[i]/xs_per_order[i])*100<<" % "<<endl;
        
        xs_tot[i]=xs_per_order[i];
        error_tot[i]=error_per_order[i];
        if(i>0)
        {
            xs_tot[i]+=xs_tot[i-1];
            error_tot[i]=sqrt(error_tot[i-1]*error_tot[i-1]+error_tot[i]*error_tot[i]);
        }
    }
    cout<<endl;
    for(int i=0;i<xs_tot.size();i++)
    {
        xs_tot[i]*=mbs[i];
        error_tot[i]*=mbs[i];
    }
    
    //Display final result
    /*cout<<endl<<endl;
    for(int i=0;i<xs_per_order.size();i++)
        cout<<"XS through N^"<<i<<"LO = "<<xs_tot[i]<<" +- "<<error_tot[i]<<" pb +- "<<fabs(error_tot[i]/xs_tot[i])*100<<" % "<<endl;
    //*/

    for(int j=1;j<result[0].size();j++)
        for(int i=0;i<result.size();i++)
        {
            result[i][j]+=result[i][j-1];
            error[i][j]=sqrt(error[i][j-1]*error[i][j-1]+error[i][j]*error[i][j]);
        }
    for(int i=0;i<result.size();i++)
        for(int j=0;j<result[i].size();j++)
        {
            result[i][j]*=mbs[j];
            error[i][j]*=mbs[j];
        }

    return;
}

vector<double> MuREvolution(double muf,double mur,const vector<double> & xs)
{
    double Lfr=log(muf*muf/mur/mur);
    vector<double> evolved;
    evolved.push_back(xs[0]);
    evolved.push_back(-2*Lfr*gammam[0]*xs[0] + xs[1]);
    evolved.push_back(pow(Lfr,2)*beta[0]*gammam[0]*xs[0] + 2*pow(Lfr,2)*pow(gammam[0],2)*xs[0] - 2*Lfr*gammam[1]*xs[0] - Lfr*beta[0]*xs[1] - 2*Lfr*gammam[0]*xs[1] + xs[2]);
    evolved.push_back((-2*pow(Lfr,3)*pow(beta[0],2)*gammam[0]*xs[0])/3. + pow(Lfr,2)*beta[1]*gammam[0]*xs[0] - 2*pow(Lfr,3)*beta[0]*pow(gammam[0],2)*xs[0] - (4*pow(Lfr,3)*pow(gammam[0],3)*xs[0])/3. + 2*pow(Lfr,2)*beta[0]*gammam[1]*xs[0] + 4*pow(Lfr,2)*gammam[0]*gammam[1]*xs[0] - 2*Lfr*gammam[2]*xs[0] + pow(Lfr,2)*pow(beta[0],2)*xs[1] - Lfr*beta[1]*xs[1] + 3*pow(Lfr,2)*beta[0]*gammam[0]*xs[1] + 2*pow(Lfr,2)*pow(gammam[0],2)*xs[1] - 2*Lfr*gammam[1]*xs[1] - 2*Lfr*beta[0]*xs[2] - 2*Lfr*gammam[0]*xs[2] + xs[3]);
    return evolved;
}


vector<double> MassEvolution(double muf,double mur,const vector<double> & xs)
{
    double Lfr=log(muf*muf/mur/mur);
    vector<double> evolved;
    evolved.push_back(xs[0]);
    evolved.push_back(-2*Lfr*gammam[0]*xs[0] + xs[1]);
    evolved.push_back(pow(Lfr,2)*beta[0]*gammam[0]*xs[0] + 2*pow(Lfr,2)*pow(gammam[0],2)*xs[0] - 2*Lfr*gammam[1]*xs[0] - 2*Lfr*gammam[0]*xs[1] + xs[2]);
    evolved.push_back((-2*pow(Lfr,3)*pow(beta[0],2)*gammam[0]*xs[0])/3. + pow(Lfr,2)*beta[1]*gammam[0]*xs[0] - 2*pow(Lfr,3)*beta[0]*pow(gammam[0],2)*xs[0] - (4*pow(Lfr,3)*pow(gammam[0],3)*xs[0])/3. + 2*pow(Lfr,2)*beta[0]*gammam[1]*xs[0] + 4*pow(Lfr,2)*gammam[0]*gammam[1]*xs[0] - 2*Lfr*gammam[2]*xs[0] + pow(Lfr,2)*beta[0]*gammam[0]*xs[1] + 2*pow(Lfr,2)*pow(gammam[0],2)*xs[1] - 2*Lfr*gammam[1]*xs[1] - 2*Lfr*gammam[0]*xs[2] + xs[3]);
    return evolved;
}













