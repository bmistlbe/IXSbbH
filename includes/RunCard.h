#ifndef Runcard_h
#define Runcard_h


#include "../includes/GlobalDefs.h"



class RunCard{
    
private:
public:
    
    double mh;
    double muf;
    double mur;
    double mb;
    double ECM;
    double mZ;
    double MCPrecision;
    string pdfstr;
    int pdfmember;
    string outputfile;
    int FONNLL;
    double Mb_OS;
    
    RunCard(){
        mh=125.09;
        muf=mh/4.0;
        mur=mh;
        pdfmember=0;
        mb=4.18;
        pdfstr="PDF4LHC15_nnlo_mc";
        ECM=13000;
        mZ=91.1876;
        outputfile="output.txt";
        MCPrecision=1e-3;
        FONNLL=0;
        Mb_OS=4.58;
    };
    
    void ReadCard(string input);
    string PrintParameters();
    
};


#endif
