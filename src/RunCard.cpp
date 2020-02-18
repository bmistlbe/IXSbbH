

#include "../includes/RunCard.h"



void RunCard::ReadCard(string input)
{
    ifstream inputf;
    string buffer="";
    
    inputf.open(input.c_str());
    if(!inputf.is_open())
    {
        cout<<"ERROR: Could not find the RunCard file: "<<input<<endl;
        exit(0);
    }
    cout<<"Reading the RunCard: "<<input<<endl;
    string locstr="";
    while(getline(inputf,buffer))
    {
        str_replace(" ","",buffer);
        str_replace("\t","",buffer);
        string var1="",var2="";
        int check=0;
        
        for(int i=0;i<buffer.size();i++)
        {
            locstr=buffer[i];
            if(locstr.compare("=")==0)
            {
                check=1;
                continue;
            }
            else if(locstr.compare("#")==0)
                break;
            if(check==0)
                var1=var1+locstr;
            if(check==1)
                var2=var2+locstr;
        }
        buffer="";
        if(check==0)
            continue;

        if(var1.compare("mH")==0)
            mh=atof(var2.c_str());
        else if(var1.compare("ECM")==0)
            ECM=atof(var2.c_str());
        else if(var1.compare("muF")==0)
            muf=atof(var2.c_str());
        else if(var1.compare("muR")==0)
            mur=atof(var2.c_str());
        else if(var1.compare("mb")==0)
            mb=atof(var2.c_str());
        else if(var1.compare("mZ")==0)
            mZ=atof(var2.c_str());
        else if(var1.compare("PDFSet")==0)
            pdfstr=var2;
        else if(var1.compare("PDFMember")==0)
            pdfmember=atof(var2.c_str());
        else if(var1.compare("MCPrecision")==0)
            MCPrecision=atof(var2.c_str());
        else if(var1.compare("FONNLL")==0)
            FONNLL=atof(var2.c_str());
        else if(var1.compare("Mb_OS")==0)
            Mb_OS=atof(var2.c_str());

    }
    inputf.close();
    return;
}


string RunCard::PrintParameters()
{
    stringstream res;
    res<<"mH="<<mh<<endl;
    res<<"ECM="<<ECM<<endl;
    res<<"muF="<<muf<<endl;
    res<<"muR="<<mur<<endl;
    res<<"mb="<<mb<<endl;
    res<<"mZ="<<mZ<<endl;
    res<<"PDFSet="<<pdfstr<<endl;
    res<<"PDFMember="<<pdfmember<<endl;
    res<<"MCPrecision="<<MCPrecision<<endl;
    res<<"FONNLL="<<FONNLL<<endl;
    res<<"Mb_OS="<<Mb_OS<<endl;
    return res.str();
}
