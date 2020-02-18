#include "../includes/CrossSection.h"


void CrossSection::SetZbCoefs()
{
    
    if(FONNLL==1)
    {
#include "XSZb_FONNLL.txt"
    }
    else
    {
#include "XSZb.txt"
    }
    return;
}

