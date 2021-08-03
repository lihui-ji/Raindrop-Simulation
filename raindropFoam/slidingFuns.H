// The prototypes of all functions defined in this file.

// Convert local index to global index, return a vector of (xi_global, yi_global, zi_global)
Vector<label> xyzIndLocalToGlobal (label npi, Vector<label> xyzIndLocal, Vector<label> dimension_local, Vector<label> dimension_processor);

// Convert global index to local index, return a list (processor index, xi_global, yi_global, zi_global)
List<label> pxyzIndGlobalToLocal  (Vector<label> xyzIndGlobal, Vector<label> dimension_local, Vector<label> dimension_processor);

// Connect amoung processors and get indices of certain phase (value > threshold)
List<DynamicList<Vector<label>>> getPhaseCellIndLst(volScalarField alpha1, Vector<label> dimension_local, Vector<label> dimension_processor, label num_proc, label proc_curr, scalar threshold);

// Calculate Index center based on stored indice from certain criteria
Vector<label> getIndCenter(List<DynamicList<Vector<label>>> phase_cell_ind_lst);

// Calculate and print phase information (center, volume, diameter)
void phaseInfo (Vector<label> ind_phase_center, List<DynamicList<Vector<label>>> phase_cell_ind_lst, List<List<vector>> cell_center_lst, List<List<scalar>> cell_volume_lst, Vector<label> dimension_local, Vector<label> dimension_processor, scalar threshold);

//print all fields information
void printAllField(fvMesh& mesh);

// Return current mesh size dimensions (xidim, yidim, zidim) with input as points (the points number is one higher than dimension number)
template<typename T>
Vector<label> findMeshDimensionPoints (T& coordinates_points);

// Return current mesh size dimensions (xidim, yidim, zidim) with input as points (the points number is one higher than dimension number)
template<typename T>
Vector<label> findMeshDimensionPoints (T& coordinates_points);

// Move all T fields in mesh (scalar, vector, sphericalTensor, symmTensor, tensor). Currently only scalar and vector fields are involved.
template<typename T>
void moveTypeFields(fvMesh& mesh, label num_proc, label proc_curr, Vector<label> xyz_ind_distance, Vector<label> dimension_global, Vector<label> dimension_local, Vector<label> dimension_processor);

// Move the given T field (scalar, vector, sphericalTensor, symmTensor, tensor)
template<typename T>
void moveTypeField(GeometricField<T, fvPatchField, volMesh> & fld_tmp, label num_proc, label proc_curr, Vector<label> xyz_ind_distance, Vector<label> dimension_global, Vector<label> dimension_local, Vector<label> dimension_processor);  

// Loop over all fields to correct boundary conditions
void correctAllBoundaryConditions(fvMesh& mesh);

// Move all fields given key_word
void moveAllField(fvMesh& mesh, label dist_zi, Vector<label> dimension, std::string key_word);

// Obtain the boundary face indice of each processor, ordered with left, right, front, back, bot, top
List<label> get_patchinds (label npi, Vector<label> dimension_processor, fvMesh& mesh);

List<DynamicList<Vector<label>>> getPhaseCellIndLst(volScalarField alpha1, Vector<label> dimension_local, Vector<label> dimension_processor, label num_proc, label proc_curr, scalar threshold){
    List<DynamicList<Vector<label>>> phase_cell_ind_lst(num_proc);
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    for (label xi = 0; xi < xidim_local; xi++){
        for (label yi = 0; yi < yidim_local; yi++){
            for (label zi = 0; zi < zidim_local; zi++){
                label ind_curr = xi + yi * xidim_local + zi * xidim_local * yidim_local;
                if (alpha1[ind_curr] > threshold){
                    Vector<label> xyz_ind_curr_local (xi, yi, zi);
                    Vector<label> xyz_ind_curr_global = xyzIndLocalToGlobal (proc_curr, xyz_ind_curr_local, dimension_local, dimension_processor);
                    phase_cell_ind_lst[proc_curr].append(xyz_ind_curr_global);
                }
            }
        }
    }
    Pstream::gatherList(phase_cell_ind_lst);
    Pstream::scatterList(phase_cell_ind_lst);    
    return phase_cell_ind_lst;
}

Vector<label> getIndCenter(List<DynamicList<Vector<label>>> phase_cell_ind_lst){
    Vector<label> ind_phase_center(0, 0, 0);
    label num_phase = 0;
    forAll(phase_cell_ind_lst, proc_i){
        forAll(phase_cell_ind_lst[proc_i], i){
            ind_phase_center += phase_cell_ind_lst[proc_i][i];
            num_phase++;
        }    
    }
    if (num_phase == 0){
        FatalErrorInFunction
        << "No cell (phase value > threshold) exists"
        << abort(FatalError);
    } 
    ind_phase_center/=num_phase;
    return ind_phase_center;
}

void phaseInfo (Vector<label> ind_phase_center, List<DynamicList<Vector<label>>> phase_cell_ind_lst, List<List<vector>> cell_center_lst, List<List<scalar>> cell_volume_lst, Vector<label> dimension_local, Vector<label> dimension_processor, scalar threshold){
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    Info << endl;
    Info << "Raindrop center grid is " << ind_phase_center << endl;
    List<label> info_phase_center_ans = pxyzIndGlobalToLocal (ind_phase_center, dimension_local, dimension_processor);
    label pi_ans = info_phase_center_ans[0], xi_ans = info_phase_center_ans[1], yi_ans = info_phase_center_ans[2], zi_ans = info_phase_center_ans[3];
    Info << "Raindrop center (m) is " << cell_center_lst[pi_ans][xi_ans + yi_ans * xidim_local + zi_ans * xidim_local * yidim_local] << endl;
    //Calculate the volume of water
    scalar volume_phase = 0;
    forAll(phase_cell_ind_lst, proc_i){
        forAll(phase_cell_ind_lst[proc_i], i){
            volume_phase += cell_volume_lst[proc_i][i];
        }    
    }
    scalar diameter_phase = Foam::cbrt(volume_phase*3/4/M_PI);
    Info << "Water volume (water fraction > " << threshold << ") is " << volume_phase << endl;
    Info << "Equivalent raindrop diameter is " << diameter_phase <<" m" << endl;
}

void moveInfo (Vector<label> xyz_ind_distance){
    Info << endl;
    Info << "Sliding Operation: ";
    if (xyz_ind_distance[0] < 0){
        Info << "move properties right " << -xyz_ind_distance[0] << " cell, ";
    }
    else if (xyz_ind_distance[0] > 0){
        Info << "move properties left " << xyz_ind_distance[0] << " cell, ";
    }
    else{
        Info << "no operation in x direction, ";
    }

    if (xyz_ind_distance[1] < 0){
        Info << "move properties back " << -xyz_ind_distance[1] << " cell, ";
    }
    else if (xyz_ind_distance[1] > 0){
        Info << "move properties front " << xyz_ind_distance[1] << " cell, ";
    }
    else{
        Info << "no operation in y direction, ";
    }

    if (xyz_ind_distance[2] < 0){
        Info << "move properties up " << -xyz_ind_distance[2] << " cell. ";
    }
    else if (xyz_ind_distance[2] > 0){
        Info << "move properties down " << xyz_ind_distance[2] << " cell. ";
    }
    else{
        Info << "no operation in z direction. ";
    }    
    Info << endl;
}

void printAllField(fvMesh& mesh){
    {typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;
    HashTable<fieldType*>
    flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        Info << "field type: " << iter()->type() << "   " << "field name: " << iter()->name() << endl;}
    }   
    {typedef GeometricField<vector, fvPatchField, volMesh> fieldType;
    HashTable<fieldType*>
    flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        Info << "field type: " << iter()->type() << "   " << "field name: " << iter()->name() << endl;}
    }

    {typedef GeometricField<scalar, fvsPatchField, surfaceMesh > fieldType;
    HashTable<fieldType*>
    flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        Info << "field type: " << iter()->type() << "   " << "field name: " << iter()->name() << endl;}
    }   

    {typedef GeometricField<vector, fvsPatchField, surfaceMesh > fieldType;
    HashTable<fieldType*>
    flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        Info << "field type: " << iter()->type() << "   " << "field name: " << iter()->name() << endl;}
    }   

    // The following fields are not in interFoam
    {typedef GeometricField<sphericalTensor, fvPatchField, volMesh> fieldType;
    HashTable<fieldType*>
    flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        Info << "field type: " << iter()->type() << "   " << "field name: " << iter()->name() << endl;}
    }
    {typedef GeometricField<symmTensor, fvPatchField, volMesh> fieldType;
    HashTable<fieldType*>
    flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        Info << "field type: " << iter()->type() << "   " << "field name: " << iter()->name() << endl;}
    }
    {typedef GeometricField<tensor, fvPatchField, volMesh> fieldType;
    HashTable<fieldType*>
    flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        Info << "field type: " << iter()->type() << "   " << "field name: " << iter()->name() << endl;}
    }
}

template<typename T>
void moveTypeFields(fvMesh& mesh, label num_proc, label proc_curr, Vector<label> xyz_ind_distance, Vector<label> dimension_global, Vector<label> dimension_local, Vector<label> dimension_processor){
    typedef GeometricField<T, fvPatchField, volMesh> fieldType;
    HashTable<fieldType*> flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
    forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {
        fieldType& fld_tmp = *iter();
        moveTypeField<T>(fld_tmp, num_proc, proc_curr, xyz_ind_distance, dimension_global, dimension_local, dimension_processor); 
    }
}

template<typename T>
void moveTypeField(GeometricField<T, fvPatchField, volMesh> & fld_tmp, label num_proc, label proc_curr, Vector<label> xyz_ind_distance, Vector<label> dimension_global, Vector<label> dimension_local, Vector<label> dimension_processor){    
    List<List<T>> field_lst(num_proc);
    field_lst[proc_curr] = fld_tmp;
    Pstream::gatherList(field_lst);
    Pstream::scatterList(field_lst);

    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    for (label xi_local = 0; xi_local < dimension_local[0]; xi_local++){
        for (label yi_local = 0; yi_local < dimension_local[1]; yi_local++){
            for (label zi_local = 0; zi_local < dimension_local[2]; zi_local++){
                label ind_curr_local = xi_local + yi_local * xidim_local + zi_local * xidim_local * yidim_local;
                Vector<label> xyz_ind_curr_local(xi_local, yi_local, zi_local);
                Vector<label> xyz_ind_curr_global = xyzIndLocalToGlobal (proc_curr, xyz_ind_curr_local, dimension_local, dimension_processor);
                Vector<label> xyz_ind_target_global = xyz_ind_curr_global + xyz_ind_distance;
                if (xyz_ind_target_global[0] < 0){xyz_ind_target_global[0] = 0;}
                if (xyz_ind_target_global[0] > dimension_global[0]-1){xyz_ind_target_global[0] = dimension_global[0]-1;}
                if (xyz_ind_target_global[1] < 0){xyz_ind_target_global[1] = 0;}
                if (xyz_ind_target_global[1] > dimension_global[1]-1){xyz_ind_target_global[1] = dimension_global[1]-1;}
                if (xyz_ind_target_global[2] < 0){xyz_ind_target_global[2] = 0;}
                if (xyz_ind_target_global[2] > dimension_global[2]-1){xyz_ind_target_global[2] = dimension_global[2]-1;}
                List<label> pxyz_ind_target_local = pxyzIndGlobalToLocal (xyz_ind_target_global, dimension_local, dimension_processor);
                label ind_target_local = pxyz_ind_target_local[1] + pxyz_ind_target_local[2] * xidim_local + pxyz_ind_target_local[3] * xidim_local * yidim_local;
                fld_tmp[ind_curr_local] = field_lst[pxyz_ind_target_local[0]][ind_target_local];
            }
        }
    } 
    
}

void correctAllBoundaryConditions(fvMesh& mesh){

        {typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;
        HashTable<fieldType*>
        flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
        forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {iter()->correctBoundaryConditions();}}   

        {typedef GeometricField<vector, fvPatchField, volMesh> fieldType;
        HashTable<fieldType*>
        flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
        forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {iter()->correctBoundaryConditions();}}

        {typedef GeometricField<sphericalTensor, fvPatchField, volMesh> fieldType;
        HashTable<fieldType*>
        flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
        forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {iter()->correctBoundaryConditions();}}

        {typedef GeometricField<symmTensor, fvPatchField, volMesh> fieldType;
        HashTable<fieldType*>
        flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
        forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {iter()->correctBoundaryConditions();}}

        {typedef GeometricField<tensor, fvPatchField, volMesh> fieldType;
        HashTable<fieldType*>
        flds_tmp(mesh.objectRegistry::lookupClass<fieldType>());
        forAllIter(typename HashTable<fieldType*>, flds_tmp, iter) {iter()->correctBoundaryConditions();}}              
}

// This function has input as points (the points number is one higher than dimension number)
template<typename T>
Vector<label> findMeshDimensionPoints (T& coordinates_points){
    label nx = 0, ny = 0, nz = 0;
    for (label i = 0; i < coordinates_points.size(); i++){
        if (nx == 0 && (i+1 == coordinates_points.size() || coordinates_points[i][0] >= coordinates_points[i+1][0])){
            nx = i;
        }
        if (ny == 0 && (i+1 == coordinates_points.size() || coordinates_points[i][2] != coordinates_points[i+1][2])){
            ny = (i + 1) /(nx + 1) - 1;
        }
    }
    nz = coordinates_points.size() /(nx + 1) /(ny + 1) - 1;
    Vector<label> dimension(nx, ny, nz);
    return dimension;
}

// This function has input as cell center coordinates (the cell number is equal to dimension number)
template<typename T>
Vector<label> findMeshDimensionCells (T& coordinates_cells){
    label nx = 0, ny = 0, nz = 0;
    for (label i = 0; i < coordinates_cells.size(); i++){
        if (nx == 0 && (i+1 == coordinates_cells.size() || coordinates_cells[i][0] >= coordinates_cells[i+1][0])){
            nx = i + 1;
        }
        if (ny == 0 && (i+1 == coordinates_cells.size() || coordinates_cells[i][2] != coordinates_cells[i+1][2])){
            ny = (i + 1) /nx;
        }
    }
    nz = coordinates_cells.size() /nx /ny;
    Vector<label> dimension(nx, ny, nz);
    return dimension;
}

// This function is to convert processor local index to gloabl index
Vector<label> xyzIndLocalToGlobal (label npi, Vector<label> xyzIndLocal, Vector<label> dimension_local, Vector<label> dimension_processor){
    label xi_local = xyzIndLocal[0], yi_local = xyzIndLocal[1], zi_local = xyzIndLocal[2];
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    label npxidim = dimension_processor[0], npyidim = dimension_processor[1], npzidim = dimension_processor[2];
    label npxi = (npi%(npxidim*npyidim))%npxidim, npyi = (npi%(npxidim*npyidim))/npxidim, npzi = npi/(npxidim*npyidim);
    label xi_global = xi_local + xidim_local * npxi, yi_global = yi_local + yidim_local * npyi, zi_global = zi_local + zidim_local * npzi;
    Vector<label> xyzIndGlobal(xi_global, yi_global, zi_global);
    return xyzIndGlobal;
}

// This function is to convert global index to processor local index
List<label> pxyzIndGlobalToLocal  (Vector<label> xyzIndGlobal, Vector<label> dimension_local, Vector<label> dimension_processor){
    label xi_global = xyzIndGlobal[0], yi_global = xyzIndGlobal[1], zi_global = xyzIndGlobal[2];
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    label npxidim = dimension_processor[0], npyidim = dimension_processor[1], npzidim = dimension_processor[2];
    label npxi = xi_global/xidim_local, npyi = yi_global/yidim_local, npzi = zi_global/zidim_local;
    label npi = npxi + npyi * npxidim + npzi * npxidim * npyidim;
    label xi_local = xi_global%xidim_local, yi_local = yi_global%yidim_local, zi_local = zi_global%zidim_local;
    List<label> pxyzIndLocal(4);
    pxyzIndLocal[0] = npi; pxyzIndLocal[1] = xi_local; pxyzIndLocal[2] = yi_local; pxyzIndLocal[3] = zi_local;
    return pxyzIndLocal;
}

// This function is to obtain the boundary face indice of each processor, ordered with left, right, front, back, bot, top
List<label> get_patchinds (label npi, Vector<label> dimension_processor, fvMesh& mesh){
    label npxidim = dimension_processor[0], npyidim = dimension_processor[1], npzidim = dimension_processor[2];
    label npxi = (npi%(npxidim*npyidim))%npxidim, npyi = (npi%(npxidim*npyidim))/npxidim, npzi = npi/(npxidim*npyidim);
    
    List<label> patchinds(6);

    if (npxi == 0){
        patchinds[0] = mesh.boundaryMesh().findPatchID("left" );
    }
    else {
        std::string a = "procBoundary";
        a.append(std::to_string(npi)).append("to").append(std::to_string(npi-1));
        patchinds[0] = mesh.boundaryMesh().findPatchID(a);
    }
    if (npxi == npxidim-1){
        patchinds[1] = mesh.boundaryMesh().findPatchID("right" );
    }
    else {
        std::string a = "procBoundary";
        a.append(std::to_string(npi)).append("to").append(std::to_string(npi+1));
        patchinds[1] = mesh.boundaryMesh().findPatchID(a);
    }
    
    if (npyi == 0){
        patchinds[2] = mesh.boundaryMesh().findPatchID("front" );
    }
    else {
        std::string a = "procBoundary";
        a.append(std::to_string(npi)).append("to").append(std::to_string(npi-npxidim));
        patchinds[2] = mesh.boundaryMesh().findPatchID(a);
    }
    if (npyi == npyidim-1){
        patchinds[3] = mesh.boundaryMesh().findPatchID("back" );
    }
    else {
        std::string a = "procBoundary";
        a.append(std::to_string(npi)).append("to").append(std::to_string(npi+npxidim));
        patchinds[3] = mesh.boundaryMesh().findPatchID(a);
    }    
    
    if (npzi == 0){
        patchinds[4] = mesh.boundaryMesh().findPatchID("bot" );
    }
    else {
        std::string a = "procBoundary";
        a.append(std::to_string(npi)).append("to").append(std::to_string(npi-npxidim*npyidim));
        patchinds[4] = mesh.boundaryMesh().findPatchID(a);
    }
    if (npzi == npzidim-1){
        patchinds[5] = mesh.boundaryMesh().findPatchID("top" );
    }
    else {
        std::string a = "procBoundary";
        a.append(std::to_string(npi)).append("to").append(std::to_string(npi+npxidim*npyidim));
        patchinds[5] = mesh.boundaryMesh().findPatchID(a);
    }   
    for (label indind = 1; indind < 6; indind++){
        if (patchinds[indind]==-1){
            Pout << patchinds;
            FatalErrorInFunction
            << "The boundary processing is incorrect (should be left, right, front, back, bot, top)"
            << abort(FatalError);
        }
    }
    return patchinds;
}