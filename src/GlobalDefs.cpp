#include "../includes/GlobalDefs.h"


double Pi=3.1415926535897932;


vector<double> Add(const vector<double> & a,const vector<double> & b)
{
    if(a.size()!=b.size())
    {
        cout<<"Trying to add to vectors of unequal length!"<<endl;
        exit(0);
    }

    vector<double> res(a.size(),0);
    for(int i=0;i<a.size();++i)
        res[i]=a[i]+b[i];
    
    return res;
}

vector<double> Multiply(const vector<double> & a,const vector<double> & b)
{
    vector<double> res;
    double loc=0;
    double size=a.size();
    if(a.size()>b.size())
        size=b.size();
    
    for(int i=0;i<a.size();++i)
    {
        loc=0;
        for(int j=0;j<=i;j++)
        {
            loc+=a[j]*b[i-j];
        }
        res.push_back(loc);
    }
    return res;
}


void str_replace(string rep, string wit, string & in) {
    int pos=0;
    if(rep.compare(wit)==0)
    return;
    while (true) {
        pos = (int)in.find(rep);
        if (pos == -1) {
            break;
        } else {
            in.erase(pos, rep.length());
            in.insert(pos, wit);
        }
    }
}
