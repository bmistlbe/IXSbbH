//
//  main.cpp
//  IXS
//
//  Created by Bernhard Mistlberger on 20/10/17.
//  Copyright © 2017 Bernhard Mistlberger. All rights reserved.
//

#include "includes/Higgs.h"
#include "includes/RunCard.h"

int main(int argc, const char * argv[]) {

    RunCard card;
    if(argc>2)
    {
        cout<<"Usage: ./main (inputfile)"<<endl;
        return 0;
    }
    else if(argc==2)
    {
        string inputfile = argv[1];
        card.ReadCard(inputfile);
    }
    
    Higgs xs;
    xs.SetMh(card.mh);
    xs.SetMuf(card.muf);
    xs.SetMur(card.mur);
    xs.SetE(card.ECM);
    xs.Setmb(card.mb);
    xs.SetPDF(card.pdfmember,card.pdfstr);
    xs.SetVerbose(0);
    xs.SetPrecision(card.MCPrecision);
    xs.SetFONNLL(card.FONNLL);
    xs.SetMBOS(card.Mb_OS);
    
    cout<<"Running with the parameters: "<<endl;
    cout<<card.PrintParameters()<<endl<<endl;
    
    vector<vector<double> > res,err;
    xs.IntegrateCrossSection();
    xs.ApplyPrefactors(res,err);
    
    stringstream ss;
    vector<string> chans={"b ab","b g","b q","b qbar","b b","g g","q qbar","q g"};
    for(int j=0;j<chans.size();j++)
        for(int i=0;i<4;i++)
            ss<<"XS at N^"<<i<<"LO for the channel "<<chans[j]<<" \t= "<<res[j][i]<<" +- "<<err[j][i]<<" pb +- "<<fabs(err[j][i]/res[j][i])*100<<" % "<<endl;
    ss<<endl;
    
    vector<double> totxs(4,0);
    vector<double> toterr(4,0);
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<chans.size();j++)
        {
            totxs[i]+=res[j][i];
            toterr[i]+=err[j][i]*err[j][i];
        }
        toterr[i]=sqrt(toterr[i]);
    }
    for(int i=0;i<4;i++)
        ss<<"xs at N"<<i<<"LO \t="<<totxs[i]<<" +- "<<toterr[i]<<" pb +- "<<toterr[i]/totxs[i]*100<<" % "<<endl;
    cout<<ss.str()<<endl;
    
    ofstream out;
    out.open(card.outputfile);
    if(!out.is_open())
    {
        cout<<"Could not open outputfile"<<card.outputfile<<endl;
        return 0;
    }
    cout<<"Storing output in "<<card.outputfile<<endl;
    out<<card.PrintParameters()<<endl<<endl;
    out<<ss.str()<<endl;
    out.close();
    
    return 0;
}
