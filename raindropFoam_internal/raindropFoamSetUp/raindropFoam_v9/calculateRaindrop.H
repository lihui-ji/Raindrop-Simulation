{
//vector raindropCenter = getPhaseCenter (alpha1, threshold);
List<DynamicList<label>> phaseCellIndList = getPhaseCellIndList(alpha1, threshold);
volVectorField cellCenter = alpha1.mesh().C();
List<DynamicList<vector>> waterCoordProcessorList = getFieldValueList<vector>(cellCenter, phaseCellIndList);
scalarField cellVolume = alpha1.mesh().V();
List<DynamicList<scalar>> waterVolProcessorList = getFieldValueList<scalar>(cellVolume, phaseCellIndList);

DynamicList<vector> waterCoordList = mergeProcessorListValue (waterCoordProcessorList);
vector raindropShape = getRaindropShape (waterCoordList);
vector raindropCenter = getAverageListValue<vector>(waterCoordList);
List<label> raindropCenterGlobalInd = getGlobalIndex (mesh, raindropCenter);
DynamicList<scalar> waterVolList = mergeProcessorListValue (waterVolProcessorList);
scalar raindropVolume = getSumListValue<scalar>(waterVolList);

vector raindropCenterVelocity = getGlobalIndexValue<vector> (U, raindropCenterGlobalInd);
Info << "The raindrop volume is " << raindropVolume << endl
     << "The raindrop center velocity is " << raindropCenterVelocity << endl
     << "The raindrop center is " << raindropCenter << endl;

if (Pstream::myProcNo() == 0)
{
    raindropOut << runTime.timeName() << '\t' 
                << raindropCenterVelocity[0] << '\t' 
                << raindropCenterVelocity[1] << '\t'
                << raindropCenterVelocity[2] << '\t'
                << mag(raindropCenterVelocity) << '\t' 
                << raindropCenter[0] << '\t' 
                << raindropCenter[1] << '\t'
                << raindropCenter[2] << '\t'
                << raindropShape[0] << '\t' 
                << raindropShape[1] << '\t'
                << raindropShape[2] << '\t'
                << raindropVolume                                                         
                << endl;
}
displacementValue = raindropCenter - initialRaindropCenter;
}
