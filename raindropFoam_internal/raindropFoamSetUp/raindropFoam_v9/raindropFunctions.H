// All customized raindrop functions are defined in this file.

// Connect amoung processors and get indices of certain phase (value > threshold)
List<DynamicList<label>> getPhaseCellIndList (volScalarField &alpha1, scalar threshold);

// Get field value list based on indices list
template <typename T>
List<DynamicList<T>> getFieldValueList (GeometricField<T, fvPatchField, volMesh> &field, List<DynamicList<label>> &CellIndList);
template <typename T>
List<DynamicList<T>> getFieldValueList (Field<T> &field, List<DynamicList<label>> &cellIndList);

// Merge all processor values
template<typename T>
DynamicList<T> mergeProcessorListValue (List<DynamicList<T>> &fieldValueList);

// Get average of list values
template <typename T>
T getAverageListValue (DynamicList<T> &fullList);

// Get sum of list values
template <typename T>
T getSumListValue (DynamicList<T> &fullList);

// Get average of all processor values
template <typename T>
T getAverageProcessorListValue (List<DynamicList<T>> &fieldValueList);

// Get the phase center coordinates as vector
vector getPhaseCenter (volScalarField &alpha1, scalar threshold);

// Get the raindrop shape (max(x)-min(x), max(y)-min(y), max(z)-min(z)) from coordinates list
vector getRaindropShape (DynamicList<vector> &waterCoordList);

// Get processorInd and index of a certain coordinate
List<label> getGlobalIndex (fvMesh &mesh, vector &raindropCenter);

// Get field value of givin globalIndex
template <typename T>
T getGlobalIndexValue (GeometricField<T, fvPatchField, volMesh> &field, List<label> &globalIndex);

// Get domain size
vector getInitialPhaseCenter ();

// Get string of parameter value from given setting file and keyword
string getSettingValue (string settingFilename, string parameterName);

// Output string of parameter value and corresponding setting keyword to given setting file
void writeSettingValue (string settingFileName, string parameterName, string value);

// Preprocess raindrop file
void preProcessRaindrop(scalar startTime);
/*************************************************************************************************************************************/

List<DynamicList<label>> getPhaseCellIndList(volScalarField &alpha1, scalar threshold)
{
    label procNum = Pstream::nProcs();
    label procCurr = Pstream::myProcNo();
    List<DynamicList<label>> phaseCellIndList(procNum);
    for(label ind = 0; ind < alpha1.size(); ind++)
    {
    	if(alpha1[ind] >= threshold)
	{
	    phaseCellIndList[procCurr].append(ind);
	}
    }
    Pstream::gatherList(phaseCellIndList);
    Pstream::scatterList(phaseCellIndList);
    return phaseCellIndList;
}

template <typename T>
List<DynamicList<T>> getFieldValueList(GeometricField<T, fvPatchField, volMesh> &field, List<DynamicList<label>> &cellIndList)
{
    label procNum = Pstream::nProcs();
    label procCurr = Pstream::myProcNo();
    List<DynamicList<T>> fieldValueList(procNum);
    forAll(cellIndList[procCurr], cellIndInd)
    {    
        fieldValueList[procCurr].append(field[cellIndList[procCurr][cellIndInd]]);
        
    }
    Pstream::gatherList(fieldValueList);
    Pstream::scatterList(fieldValueList);
    return fieldValueList;
}

template <typename T>
List<DynamicList<T>> getFieldValueList(Field<T> &field, List<DynamicList<label>> &cellIndList)
{
    label procNum = Pstream::nProcs();
    label procCurr = Pstream::myProcNo();
    List<DynamicList<T>> fieldValueList(procNum);
    forAll(cellIndList[procCurr], cellIndInd)
    {    
        fieldValueList[procCurr].append(field[cellIndList[procCurr][cellIndInd]]);
        
    }
    Pstream::gatherList(fieldValueList);
    Pstream::scatterList(fieldValueList);
    return fieldValueList;
}

template<typename T>
DynamicList<T> mergeProcessorListValue (List<DynamicList<T>> &fieldValueList)
{
    DynamicList<T> fullList;
    forAll(fieldValueList, processorInd)
    {
    	fullList.append(fieldValueList[processorInd]);
    }
    return fullList;    
}

template <typename T>
T getAverageListValue (DynamicList<T> &fullList)
{
    T averageProcessorListValue = fullList[0]-fullList[0];
    label fullListSize = fullList.size();
    forAll(fullList, fullListInd)
    {
    	averageProcessorListValue += fullList[fullListInd]/fullListSize;
    }
    return averageProcessorListValue;
}

template <typename T>
T getSumListValue (DynamicList<T> &fullList)
{
    T sumProcessorListValue = fullList[0]-fullList[0];
    forAll(fullList, fullListInd)
    {
    	sumProcessorListValue += fullList[fullListInd];
    }
    return sumProcessorListValue;
}

template <typename T>
T getAverageProcessorListValue(List<DynamicList<T>> &fieldValueList)
{
    DynamicList<T> fullList;
    forAll(fieldValueList, processorInd)
    {
    	fullList.append(fieldValueList[processorInd]);
    }
    T averageProcessorListValue = fullList[0]-fullList[0];
    label fullListSize = fullList.size();
    forAll(fullList, fullListInd)
    {
    	averageProcessorListValue += fullList[fullListInd]/fullListSize;
    }
    return averageProcessorListValue;
}

vector getPhaseCenter(volScalarField &alpha1, scalar threshold)
{
    List<DynamicList<label>> phaseCellIndList = getPhaseCellIndList(alpha1, threshold);
    volVectorField cellCenter = alpha1.mesh().C();
    List<DynamicList<vector>> fieldValueList = getFieldValueList<vector>(cellCenter, phaseCellIndList);
    vector raindropCenter = getAverageProcessorListValue<vector>(fieldValueList);
    return raindropCenter;
}

vector getRaindropShape (DynamicList<vector> &waterCoordList)
{
    vector coordMax (waterCoordList[0]);
    vector coordMin (waterCoordList[0]);
    forAll (waterCoordList, coordInd)
    {
        for (label xyzInd = 0; xyzInd < 3; xyzInd++)
        {
            coordMax[xyzInd] = max (coordMax[xyzInd],waterCoordList[coordInd][xyzInd]);
            coordMin[xyzInd] = min (coordMin[xyzInd],waterCoordList[coordInd][xyzInd]);
        }
    }
    return coordMax - coordMin;
}

List<label> getGlobalIndex (fvMesh &mesh, vector &raindropCenter)
{
    label procNum = Pstream::nProcs();
    label procCurr = Pstream::myProcNo();
    List<label> indexList(procNum);
    List<scalar> normDistance(procNum);
    label searchedInd = mesh.findCell(raindropCenter);
    indexList[procCurr] = searchedInd;
    normDistance[procCurr] = searchedInd==-1 ? -1 : mag(mesh.C()[searchedInd]-raindropCenter);
    scalar minNormDistance = -1;
    Pstream::gatherList(indexList);
    Pstream::scatterList(indexList);
    Pstream::gatherList(normDistance);
    Pstream::scatterList(normDistance);
    List<label> globalIndex(2);
    for (label processorInd = 0; processorInd < procNum; processorInd++) 
    {
    	if (normDistance[processorInd] > 0)
	{	
	    if (minNormDistance == -1 || minNormDistance > normDistance[processorInd])
	    {
	    	minNormDistance = normDistance[processorInd];
		globalIndex[0] = processorInd;
		globalIndex[1] = indexList[processorInd];
	    }  
	}
    }
    return globalIndex;
}

vector getInitialPhaseCenter()
{
    vector phaseCenter;
    std::string line;
    std::ifstream myfile ("system/setFieldsDict");
    if (myfile.is_open())
    {
        while ( std::getline (myfile,line) )
        {
            if (line.find("centre") != std::string::npos)
            {
                line = line.substr(line.find('(')+1, line.find(')')-line.find('(')-1);
                phaseCenter[0] = std::stod(line.substr(0, line.find(' ')));
                line = line.substr(line.find(' ')+1, line.size()-line.find(' '));
                phaseCenter[1] = std::stod (line.substr(0, line.find(' ')));
                line = line.substr(line.find(' ')+1, line.size()-line.find(' '));
                phaseCenter[2] = std::stod(line);
                return phaseCenter;
            }
        }
        myfile.close();
    }
    else
    {
        cout << "Can not open file.";
    }
    cout << "Can not find phase centre " << endl;
    return phaseCenter;
}    

template <typename T>
T getGlobalIndexValue (GeometricField<T, fvPatchField, volMesh> &field, List<label> &globalIndex)
{
    List<T> processorValueList (Pstream::nProcs());
    if (globalIndex[0] == Pstream::myProcNo())
    {
    	processorValueList[globalIndex[0]] = field[globalIndex[1]];
    }
    Pstream::gatherList (processorValueList);
    Pstream::scatterList (processorValueList);
    return processorValueList[globalIndex[0]];
}

string getSettingValue (string settingFileName, string parameterName)
{
    std::string line;
    std::ifstream myfile (settingFileName);
    if (myfile.is_open())
    {
        while ( std::getline (myfile,line) )
        {
            if (line.find(parameterName)==0)
            {
            size_t numPos = 0;
            for (size_t linei = 0; linei < line.size(); linei++)
            {
                if ((line[linei] >= '0' && line[linei] <= '9') || line[linei] == '-')
                {
                    numPos = linei;
                    break;
                }
            }
            return line.substr(numPos, line.find(';')-numPos);
            }
        }
        myfile.close();
    }
    else
    {
        FatalErrorInFunction << "Can not open file." << abort(FatalError);                       
    }
    FatalErrorInFunction << "Can not find parameter " << parameterName << abort(FatalError); 
    return line;
}          

void writeSettingValue (string settingFileName, string parameterName, string value)
{
    std::string line;
    std::string temp;
    std::ifstream myfile (settingFileName);
    bool if_found = false;
    if (myfile.is_open())
    {
        while ( std::getline (myfile,line) )
        {
            if (line.find(parameterName)!=0)
            {
                temp.append(line.append("\n"));
            }
            else
            {
                temp.append(parameterName.append("\t").append(value).append(";\n"));
                if_found = true;
            }
        }
        if (if_found == false)
        {
            temp.append("\n");
            temp.append(parameterName.append("\t").append(value).append(";\n"));            
        }
        myfile.close();
    }
    else
    {
        FatalErrorInFunction << "Can not open file." << abort(FatalError);                       
    }

    std::ofstream out(settingFileName);
    out << temp;
    out.close();
    
}  

void preProcessRaindrop(scalar startTime)
{
    if (Pstream::myProcNo() !=0 )
    {
        return;
    }
    std::string line;
    // open input file
    std::ifstream in("raindrop");
    if( !in.is_open())
    {
          FatalErrorInFunction << "Fail to open raindrop file.\n" << abort(FatalError);
    }
    
    std::string temp;
    // write the header
    std::getline(in,line);
    temp.append(line.append("\n"));
    // loop to read/write the file.
    while( std::getline(in,line) )
    {        
        if (line[0] >= '0' && line[0] <= '9')
        {
            if (std::stod(line.substr(0, line.find('\t')+1)) >= startTime)
            {
                continue;
            }
            temp.append(line.append("\n"));
        }
    }
    in.close();
    std::ofstream out("raindrop");
    out << temp;
    out.close();
}