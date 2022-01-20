// Declare dynamic parameters
// --- Declare displacementValue for mesh moving
vector displacementValue;

if (startTime == 0)
{
    Info << "Creating velocity output file" << endl;    
}
else
{
    Info << "Preprocessing velocity output file" << endl;
    preProcessRaindrop(startTime);
}

Info << "Threshold = " << threshold << endl;

OFstream raindropOut ("raindrop", Foam::IOstream::ASCII, Foam::IOstream::currentVersion, Foam::IOstream::UNCOMPRESSED, true);
if (startTime == 0)
{
    if (Pstream::myProcNo() == 0)
    {
        raindropOut << "Time" << '\t' << "Ux\tUy\tUz\tUmag\tx\ty\tz\ta\tb\tc\tV" << endl;
    }
}   
#include "calculateRaindrop.H"
