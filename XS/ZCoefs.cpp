#include "../includes/CrossSection.h"


void CrossSection::SetZCoefs()
{
    if(FONNLL)
    {
#include "XSZ_FONNLL.txt"
    }
    else
    {
#include "XSZ.txt"
    }
    return;
}
