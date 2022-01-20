scalar threshold = std::stod(getSettingValue ("system/controlDict", "threshold"));
//scalar startTime = std::stod(getSettingValue ("system/controlDict", "startTime"));
scalar startTime = runTime.startTime().value();
Info << "Simulation starts at " << startTime << "s" << endl;
List <vector> processorCenterList (Pstream::nProcs());
vector processorCenter = 0.5*(mesh.points()[0] + mesh.points()[mesh.points().size()-1]); 
processorCenterList[Pstream::myProcNo()] = processorCenter;
Pstream::gatherList(processorCenterList);
Pstream::scatterList(processorCenterList);

vector domainCenter (0., 0., 0.);
forAll(processorCenterList, processorInd)
{
    domainCenter += processorCenterList[processorInd]/Pstream::nProcs();
}

Info << "Domain center is " << domainCenter << endl;

vector initialRaindropCenter = getPhaseCenter(alpha1, threshold);
if (startTime !=0)
{
    Info << "The initial phase center is read from system/setFields" << endl;
    initialRaindropCenter = getInitialPhaseCenter();
    Info << "Initial phase center is " << initialRaindropCenter << endl;
}

