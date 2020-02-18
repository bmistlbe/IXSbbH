#include "../includes/CrossSection.h"


void CrossSection::SetWCoefs()
{
    if(FONNLL==1)
    {
#include "XSW_FONNLL.txt"
    }
    else
    {
#include "XSW.txt"
    }
    return;
}

