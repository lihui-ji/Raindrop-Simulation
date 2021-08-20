// The prototypes of all functions defined in this file.

// Convert local index to global index, return a vector of (xi_global, yi_global, zi_global)
Vector<label> xyzIndLocalToGlobal (label npi, Vector<label> xyzIndLocal, Vector<label> dimension_local, Vector<label> dimension_processor);

// Convert global index to local index, return a list (processor index, xi_global, yi_global, zi_global)
List<label> pxyzIndGlobalToLocal  (Vector<label> xyzIndGlobal, Vector<label> dimension_local, Vector<label> dimension_processor);

Vector<label> pxyzTargetIndCellToIndFace (List<label> pxyz_ind_local, Vector<label> dimension_local, List<Vector<label>> cell_face_dict, List<List<List<label>>> boundary_order_dict_lst);

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

void moveSurfaceScalarField(surfaceScalarField &fld_tmp, label num_proc, label proc_curr, Vector<label> xyz_ind_distance, Vector<label> dimension_global, Vector<label> dimension_local, Vector<label> dimension_processor, List<label> patchinds, Vector<label> xyz_onesidebdr_face_num, label internal_face_num, List<Vector<label>> cell_face_dict, List<List<List<label>>> boundary_face_order_dict_lst){
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    label x_onesidebdr_face_num = xyz_onesidebdr_face_num[0], y_onesidebdr_face_num = xyz_onesidebdr_face_num[1], z_onesidebdr_face_num = xyz_onesidebdr_face_num[2];

    DynamicList<scalar> faces_value;
    DynamicList<label> empty_patch_inds;
    faces_value.append(fld_tmp.internalField());
    forAll(patchinds, ind_tmp){
        List<scalar> faces_value_curr = fld_tmp.boundaryField()[patchinds[ind_tmp]];
        // this is to solve empty boundary conditions
        if (faces_value_curr.size()==0){
            empty_patch_inds.append(ind_tmp);
            Info << "Empty boundary faces in " << fld_tmp.name() << endl;
            List<scalar> zero_empty(xyz_onesidebdr_face_num[ind_tmp/2]);
            forAll(zero_empty, zero_ind){
                zero_empty[zero_ind] = 0;
            }
            faces_value_curr = zero_empty;
        }
        faces_value.append(faces_value_curr);
    }
    if (empty_patch_inds.size()!=2 && empty_patch_inds.size()!=0){
            FatalErrorInFunction
            << "The empty boundary faces number is incorrect (0 or 2)"
            << abort(FatalError);
    }

    List<List<scalar> > faces_value_lst(num_proc);
    faces_value_lst[proc_curr] = faces_value;
    Pstream::gatherList(faces_value_lst);
    Pstream::scatterList(faces_value_lst);

    for (label xi_local = -1; xi_local < dimension_local[0]; xi_local++){
        for (label yi_local = -1; yi_local < dimension_local[1]; yi_local++){
            for (label zi_local = -1; zi_local < dimension_local[2]; zi_local++){
                List<label> pxyz_ind_curr_local(4);
                pxyz_ind_curr_local[0] = proc_curr; pxyz_ind_curr_local[1] = xi_local; pxyz_ind_curr_local[2] = yi_local; pxyz_ind_curr_local[3] = zi_local;
                Vector<label> xyz_ind_curr_local (xi_local, yi_local, zi_local);
                Vector<label> faces_ind_curr_local = pxyzTargetIndCellToIndFace (pxyz_ind_curr_local, dimension_local, cell_face_dict, boundary_face_order_dict_lst);
                Vector<label> xyz_ind_curr_global = xyzIndLocalToGlobal (proc_curr, xyz_ind_curr_local, dimension_local, dimension_processor);
                Vector<label> xyz_ind_target_global = xyz_ind_curr_global + xyz_ind_distance;
                // The face boundary assignment is different from cell center boundary assignment
                if (xyz_ind_target_global[0] < 0){xyz_ind_target_global[0] = -1;}
                if (xyz_ind_target_global[0] > dimension_global[0]-1){xyz_ind_target_global[0] = dimension_global[0]-1;}
                if (xyz_ind_target_global[1] < 0){xyz_ind_target_global[1] = -1;}
                if (xyz_ind_target_global[1] > dimension_global[1]-1){xyz_ind_target_global[1] = dimension_global[1]-1;}
                if (xyz_ind_target_global[2] < 0){xyz_ind_target_global[2] = -1;}
                if (xyz_ind_target_global[2] > dimension_global[2]-1){xyz_ind_target_global[2] = dimension_global[2]-1;}
                List<label> pxyz_ind_target_local = pxyzIndGlobalToLocal (xyz_ind_target_global, dimension_local, dimension_processor);
                Vector<label> xyz_ind_target_local (pxyz_ind_target_local[1], pxyz_ind_target_local[2], pxyz_ind_target_local[3]);
                Vector<label> faces_ind_target_local = pxyzTargetIndCellToIndFace (pxyz_ind_target_local, dimension_local, cell_face_dict, boundary_face_order_dict_lst);

                Vector<scalar> faces_value_target;
                forAll(faces_value_target, faces_value_target_ind_tmp){
                    if (xyz_ind_target_local[faces_value_target_ind_tmp] == -1){
                        faces_value_target[faces_value_target_ind_tmp] = -faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[faces_value_target_ind_tmp]];
                    }
                    else{
                        faces_value_target[faces_value_target_ind_tmp] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[faces_value_target_ind_tmp]];
                    }
                }
                
                //Pout << faces_ind_curr_local << faces_ind_target_local << endl;
                List<scalar> faces_value_new(faces_value);
                for (label ind_tmp = 0; ind_tmp < 3; ind_tmp++){
                    if (xi_local == -1 && ind_tmp != 0){continue;}
                    if (yi_local == -1 && ind_tmp != 1){continue;}
                    if (zi_local == -1 && ind_tmp != 2){continue;}
                    label face_ind_curr_local = faces_ind_curr_local[ind_tmp];
                    if (face_ind_curr_local < internal_face_num){
                        fld_tmp[face_ind_curr_local] = faces_value_target[ind_tmp];
                    }
                    else if (face_ind_curr_local < internal_face_num + x_onesidebdr_face_num){
                        if (empty_patch_inds[0] == 0 || empty_patch_inds[1] == 0){continue;}                        
                        fld_tmp.boundaryFieldRef()[patchinds[0]][face_ind_curr_local - internal_face_num] = -faces_value_target[ind_tmp];
                    }
                    else if (face_ind_curr_local < internal_face_num + x_onesidebdr_face_num * 2){
                        if (empty_patch_inds[0] == 1 || empty_patch_inds[1] == 1){continue;}  
                        fld_tmp.boundaryFieldRef()[patchinds[1]][face_ind_curr_local - internal_face_num - x_onesidebdr_face_num] = faces_value_target[ind_tmp];
                    }
                    else if (face_ind_curr_local < (internal_face_num + x_onesidebdr_face_num * 2 + y_onesidebdr_face_num)){
                        if (empty_patch_inds[0] == 2 || empty_patch_inds[1] == 2){continue;}  
                        fld_tmp.boundaryFieldRef()[patchinds[2]][face_ind_curr_local - internal_face_num - x_onesidebdr_face_num * 2] = -faces_value_target[ind_tmp];
                    }
                    else if (face_ind_curr_local < (internal_face_num + x_onesidebdr_face_num * 2 + y_onesidebdr_face_num * 2)){
                        if (empty_patch_inds[0] == 3 || empty_patch_inds[1] == 3){continue;}  
                        fld_tmp.boundaryFieldRef()[patchinds[3]][face_ind_curr_local - internal_face_num - x_onesidebdr_face_num * 2 - y_onesidebdr_face_num] = faces_value_target[ind_tmp];
                    }
                    else if (face_ind_curr_local < (internal_face_num + x_onesidebdr_face_num * 2 + y_onesidebdr_face_num * 2 + z_onesidebdr_face_num)){
                        if (empty_patch_inds[0] == 4 || empty_patch_inds[1] == 4){continue;}  
                        fld_tmp.boundaryFieldRef()[patchinds[4]][face_ind_curr_local - internal_face_num - x_onesidebdr_face_num * 2 - y_onesidebdr_face_num * 2] = -faces_value_target[ind_tmp];
                    }
                    else{
                        if (empty_patch_inds[0] == 5 || empty_patch_inds[1] == 5){continue;}  
                        fld_tmp.boundaryFieldRef()[patchinds[5]][face_ind_curr_local - internal_face_num - x_onesidebdr_face_num * 2 - y_onesidebdr_face_num * 2 - z_onesidebdr_face_num] = faces_value_target[ind_tmp];
                    }
                }
                // faces_value_new[faces_ind_curr_local[0]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[0]];
                // faces_value_new[faces_ind_curr_local[1]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[1]];
                // faces_value_new[faces_ind_curr_local[2]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[2]];
                
                // Pout << "xyz_ind_curr_local" << xyz_ind_curr_local << "xyz_ind_curr_global" << xyz_ind_curr_global << "xyz_ind_target_global" << xyz_ind_target_global << "pxyz_ind_target_local" << pxyz_ind_target_local << endl
                //      << faces_value[faces_ind_curr_local[0]] << faces_value_new[faces_ind_curr_local[0]] << endl
                //      << faces_value[faces_ind_curr_local[1]] << faces_value_new[faces_ind_curr_local[1]] << endl
                //      << faces_value[faces_ind_curr_local[2]] << faces_value_new[faces_ind_curr_local[2]] << endl;
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
    label npxidim = dimension_processor[0], npyidim = dimension_processor[1], npzidim = dimension_processor[2];
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    label xi_global = xyzIndGlobal[0], yi_global = xyzIndGlobal[1], zi_global = xyzIndGlobal[2];
    label npxi, npyi, npzi, xi_local, yi_local, zi_local;

    if (xi_global < 0){npxi = 0; xi_local = -1;}
    //else if (xi_global == xidim_local * npxidim){npxi = npxidim-1; xi_local = xidim_local+1;}
    else{npxi = xi_global/xidim_local; xi_local = xi_global%xidim_local;}

    if (yi_global < 0){npyi = 0; yi_local = -1;}
    //else if (yi_global == yidim_local * npyidim){npyi = npyidim-1; yi_local = yidim_local+1;}
    else{npyi = yi_global/yidim_local; yi_local = yi_global%yidim_local;}

    if (zi_global < 0){npzi = 0; zi_local = -1;}
    //else if (zi_global == zidim_local * npzidim){npzi = npzidim-1; zi_local = zidim_local+1;}
    else{npzi = zi_global/zidim_local; zi_local = zi_global%zidim_local;}

    label npi = npxi + npyi * npxidim + npzi * npxidim * npyidim;
    List<label> pxyzIndLocal(4);
    pxyzIndLocal[0] = npi; pxyzIndLocal[1] = xi_local; pxyzIndLocal[2] = yi_local; pxyzIndLocal[3] = zi_local;
    return pxyzIndLocal;
}

Vector<label> pxyzTargetIndCellToIndFace (List<label> pxyz_ind_local, Vector<label> dimension_local, List<Vector<label>> cell_face_dict, List<List<List<label>>> boundary_order_dict_lst){
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    label npi = pxyz_ind_local[0], xi_local = pxyz_ind_local[1], yi_local = pxyz_ind_local[2], zi_local = pxyz_ind_local[3];
    label internal_face_num = yidim_local*zidim_local*(xidim_local-1) + xidim_local*zidim_local*(yidim_local-1) + xidim_local*yidim_local*(zidim_local-1);
    label x_onesidebdr_face_num = yidim_local*zidim_local,  y_onesidebdr_face_num = xidim_local*zidim_local, z_onesidebdr_face_num = xidim_local*yidim_local;
    label face0_ind_local, face1_ind_local, face2_ind_local;
    Vector<label> initial_xyz_ind_local (xi_local, yi_local, zi_local);
    forAll(initial_xyz_ind_local, ind_tmp){
        initial_xyz_ind_local[ind_tmp] = (initial_xyz_ind_local[ind_tmp] != -1) * initial_xyz_ind_local[ind_tmp];
    }
    label initial_ind_local = initial_xyz_ind_local[0] + initial_xyz_ind_local[1] * xidim_local + initial_xyz_ind_local[2] * xidim_local * yidim_local;
    Vector<label> faces_ind_local = cell_face_dict[initial_ind_local];
    if (xi_local == -1){
        faces_ind_local[0] = boundary_order_dict_lst[npi][0][(yi_local != -1) * yi_local + (zi_local != -1) * zi_local * yidim_local];
    }
    if (xi_local == xidim_local-1){
        faces_ind_local[0] = boundary_order_dict_lst[npi][1][(yi_local != -1) * yi_local + (zi_local != -1) * zi_local * yidim_local];
    }

    if (yi_local == -1){
        faces_ind_local[1] = boundary_order_dict_lst[npi][2][(xi_local != -1) * xi_local + (zi_local != -1) * zi_local * xidim_local];
    }
    if (yi_local == yidim_local-1){
        faces_ind_local[1] = boundary_order_dict_lst[npi][3][(xi_local != -1) * xi_local + (zi_local != -1) * zi_local * xidim_local];
    }
    if (zi_local == -1){
        faces_ind_local[2] = boundary_order_dict_lst[npi][4][(xi_local != -1) * xi_local + (yi_local != -1) * yi_local * xidim_local];
    }
    if (zi_local == zidim_local-1){
        faces_ind_local[2] = boundary_order_dict_lst[npi][5][(xi_local != -1) * xi_local + (yi_local != -1) * yi_local * xidim_local];
    }
    return faces_ind_local;
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

// This function is to return boundary patch face order
List<label> get_patch_order_dict(List<vector> faces_value_curr){
    
    label size = faces_value_curr.size();
    List<label> patch_order_dict(size);
    for (label ind_tmp = 0; ind_tmp < size; ind_tmp++){
        patch_order_dict[ind_tmp] = ind_tmp;       
    }

    for (int i =0; i < size - 1; i++){
        for (int j = 0; j < size - 1; j++){
            label if_correct;
            for (label xyzind_tmp = 2; xyzind_tmp >= 0; xyzind_tmp--){
                if (faces_value_curr[patch_order_dict[j]][xyzind_tmp] > faces_value_curr[patch_order_dict[j+1]][xyzind_tmp]){
                    if_correct = 0;
                    break;
                }
                else if (faces_value_curr[patch_order_dict[j]][xyzind_tmp] < faces_value_curr[patch_order_dict[j+1]][xyzind_tmp]){
                    if_correct = 1;
                    break;
                }
                else{if_correct = -1;}
            }
            if (if_correct == 0){
                label tmp = patch_order_dict[j];
                patch_order_dict[j] = patch_order_dict[j+1];
                patch_order_dict[j+1] = tmp;
            }
        }
    }
// check if the order is correct
    List<vector> faces_value_ordered(size);
    for (int i =0; i < size - 1; i++){
        faces_value_ordered[i] = faces_value_curr[patch_order_dict[i]];
    }

    for (int i =0; i < size - 1; i++){
        label if_correct;
        vector diff = faces_value_curr[patch_order_dict[i+1]] - faces_value_curr[patch_order_dict[i]];
        for (label xyzind = 2; xyzind >= 0; xyzind--){
            if (diff[xyzind] > 0){
                break;
            }
            if (diff[xyzind] < 0){

                Pout << "!!!" << faces_value_curr[patch_order_dict[i+1]] << faces_value_curr[patch_order_dict[i]] << endl;
                Pout << "!!!" << faces_value_ordered << endl;
                FatalErrorInFunction
                << "The order of patch face is incorrect"
                << abort(FatalError);
            }
        }
    }
    //Info << patch_order_dict << endl;
    return patch_order_dict;
}

// Get boundary face indices
List<List<label>> get_boundary_face_order_dict(fvMesh& mesh, List<label> patchinds, label internal_face_num, Vector<label> xyz_onesidebdr_face_num){
    label x_onesidebdr_face_num = xyz_onesidebdr_face_num[0], y_onesidebdr_face_num = xyz_onesidebdr_face_num[1], z_onesidebdr_face_num = xyz_onesidebdr_face_num[2];
    List<List<label>> boundary_face_order_dict(6);
    forAll(patchinds, ind_tmp){
        List<vector> faces_value_curr = mesh.Cf().boundaryField()[patchinds[ind_tmp]];
        // this is to solve empty boundary conditions
        if (faces_value_curr.size()==0){
            //Info << "Empty boundary faces in " << mesh.Cf().name() << endl;
            List<vector> zero_empty(xyz_onesidebdr_face_num[ind_tmp/2], vector(0, 0, 0));
            faces_value_curr = zero_empty;
        }
        boundary_face_order_dict[ind_tmp] = get_patch_order_dict(faces_value_curr);
        //Info << mesh.Cf().boundaryField()[patchinds[ind_tmp]].size() << " ";
    }    

    for (label ind_tmp = 0; ind_tmp < 6; ind_tmp++){
        label boundary_ind_base = 0;
        if (ind_tmp == 1){boundary_ind_base = x_onesidebdr_face_num;}
        if (ind_tmp == 2){boundary_ind_base = x_onesidebdr_face_num * 2;}
        if (ind_tmp == 3){boundary_ind_base = x_onesidebdr_face_num * 2 + y_onesidebdr_face_num;}
        if (ind_tmp == 4){boundary_ind_base = x_onesidebdr_face_num * 2 + y_onesidebdr_face_num * 2;}
        if (ind_tmp == 5){boundary_ind_base = x_onesidebdr_face_num * 2 + y_onesidebdr_face_num * 2 + z_onesidebdr_face_num;}
        forAll(boundary_face_order_dict[ind_tmp], ind_patch_tmp){
            boundary_face_order_dict[ind_tmp][ind_patch_tmp] += (internal_face_num + boundary_ind_base);
        }   
    }
    return boundary_face_order_dict; 
}


List<Vector<label>> get_cell_face_dict(Vector<label> dimension_local, label internal_face_num, Vector<label> xyz_onesidebdr_face_num){
    label x_onesidebdr_face_num = xyz_onesidebdr_face_num[0], y_onesidebdr_face_num = xyz_onesidebdr_face_num[1], z_onesidebdr_face_num = xyz_onesidebdr_face_num[2];
    label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
    label zi_face_num = xidim_local*yidim_local+(xidim_local-1)*yidim_local+(yidim_local-1)*xidim_local;
    List<Vector<label>> cell_face_dict(xidim_local*yidim_local*zidim_local);
    for (label zi_local = 0; zi_local < zidim_local; zi_local++){
        for (label yi_local = 0; yi_local < yidim_local; yi_local++){
            for (label xi_local = 0; xi_local < xidim_local; xi_local++){
                label face_curr = zi_local * zi_face_num;
                label cell_curr = zi_local * xidim_local * yidim_local + yi_local * xidim_local + xi_local;
                if (zi_local == zidim_local-1){
                    face_curr += yi_local * (xidim_local * 2-1);
                    if (yi_local == yidim_local-1){
                        face_curr += xi_local;
                    }
                    else{
                        face_curr += xi_local*2;
                    }
                }
                else{
                    face_curr += yi_local * (xidim_local * 3-1);
                    if (yi_local == yidim_local-1){
                        face_curr += xi_local*2;
                    }
                    else{
                        face_curr += xi_local*3;
                    }
                }
                Vector<label> face_inds_curr(-1, -1, -1);
                if (xi_local == xidim_local-1){
                    face_inds_curr[0] = internal_face_num +  x_onesidebdr_face_num + zi_local*yidim_local + yi_local;
                }
                if (yi_local == yidim_local-1){
                    face_inds_curr[1] = internal_face_num +  x_onesidebdr_face_num*2 + y_onesidebdr_face_num + zi_local*xidim_local + xi_local;
                }
                if (zi_local == zidim_local-1){
                    face_inds_curr[2] = internal_face_num +  x_onesidebdr_face_num*2 + y_onesidebdr_face_num*2 + z_onesidebdr_face_num + yi_local*xidim_local + xi_local;
                }
                forAll(face_inds_curr, ind_tmp){
                    if(face_inds_curr[ind_tmp] == -1){
                        face_inds_curr[ind_tmp] = face_curr;
                        face_curr++;
                    }
                }
                cell_face_dict[cell_curr] = face_inds_curr;
            }
        }
    }
    return cell_face_dict;
}
