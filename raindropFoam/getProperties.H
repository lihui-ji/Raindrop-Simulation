label num_proc = Pstream::nProcs();
label proc_curr = Pstream::myProcNo();

List<List<scalar>> cell_volume_lst (num_proc);
cell_volume_lst[Pstream::myProcNo()] = mesh.V();
Pstream::gatherList(cell_volume_lst);
Pstream::scatterList(cell_volume_lst);

List<List<vector>> cell_center_lst (num_proc);
cell_center_lst[Pstream::myProcNo()] = mesh.C();
Pstream::gatherList(cell_center_lst);
Pstream::scatterList(cell_center_lst);

List<vector> dimension_lst (num_proc);
Vector<label> dimension_local = findMeshDimensionPoints (mesh.points());
label xidim_local = dimension_local[0], yidim_local = dimension_local[1], zidim_local = dimension_local[2];
dimension_lst[Pstream::myProcNo()] = dimension_local;
Pstream::gatherList(dimension_lst);
Pstream::scatterList(dimension_lst);
forAll(dimension_lst, i){
    if (xidim_local != dimension_lst[i][0] || yidim_local != dimension_lst[i][1] || zidim_local != dimension_lst[i][2]){
        FatalErrorInFunction
        << "The decomposition is uneven."
        << abort(FatalError);
    }
}
label cell_num_local = xidim_local * yidim_local * zidim_local;
label cell_num_global = cell_num_local * num_proc;

List<vector> domain_center_local_lst (num_proc);
vector domain_center_local = 0.5*(mesh.points()[0] + mesh.points()[cell_num_local-1]);
domain_center_local_lst[Pstream::myProcNo()] = domain_center_local;
Pstream::gatherList(domain_center_local_lst);
Pstream::scatterList(domain_center_local_lst);
vector domain_center_global (0., 0., 0.);
forAll(domain_center_local_lst, i){
    domain_center_global += domain_center_local_lst[i];
}

domain_center_global /= domain_center_local_lst.size();
Vector <label> dimension_processor = findMeshDimensionCells (domain_center_local_lst);
label npxidim = dimension_processor[0], npyidim = dimension_processor[1], npzidim = dimension_processor[2];
label xidim_global = xidim_local * npxidim, yidim_global = yidim_local * npyidim, zidim_global = zidim_local * npzidim;
label xi_center_global = (xidim_global-1)/2, yi_center_global = (yidim_global-1)/2, zi_center_global = (zidim_global-1)/2;
label ximax_global = xidim_global + 1, yimax_global= yidim_global + 1, zimax_global = zidim_global + 1;
Vector<label> dimension_global (xidim_global, yidim_global, zidim_global);
Vector<label> ind_center_global (xi_center_global, yi_center_global, zi_center_global);
Vector<label> ind_max_global (ximax_global, yimax_global, zimax_global);

Info << "Domain physical center is " << domain_center_global << " m" << endl;
Info << "Domain mesh dimension is " << dimension_global << endl;
Info << "The global mesh grid center is " << ind_center_global << endl;
Info << "Processor layout is " << dimension_processor <<endl;
Info << "For each processor, the mesh dimension is " << dimension_local << endl;
Info << endl;
printAllField(mesh);
Info << endl;

List<label> patchinds = get_patchinds (proc_curr, dimension_processor, mesh);
List<List<label> > patchinds_lst(num_proc);
patchinds_lst[proc_curr] = patchinds;
Pstream::gatherList(patchinds_lst);
Pstream::scatterList(patchinds_lst);
Info << "The boundary indices in each processor:" << endl << patchinds_lst << endl;

label internal_face_num = yidim_local*zidim_local*(xidim_local-1) + xidim_local*zidim_local*(yidim_local-1) + xidim_local*yidim_local*(zidim_local-1);
label x_onesidebdr_face_num = yidim_local*zidim_local;
label y_onesidebdr_face_num = xidim_local*zidim_local;
label z_onesidebdr_face_num = xidim_local*yidim_local;
Vector<label> xyz_onesidebdr_face_num(x_onesidebdr_face_num, y_onesidebdr_face_num, z_onesidebdr_face_num);

List<List<List<label>>> boundary_face_order_dict_lst (num_proc);
boundary_face_order_dict_lst[proc_curr] = get_boundary_face_order_dict(mesh, patchinds, internal_face_num, xyz_onesidebdr_face_num);
Pstream::gatherList(boundary_face_order_dict_lst);
Pstream::scatterList(boundary_face_order_dict_lst);

List<Vector<label>> cell_face_dict = get_cell_face_dict(dimension_local, internal_face_num, xyz_onesidebdr_face_num);


Vector<label> xyz_ind_distance(3, 0, -2);
phi = fvc::interpolate(U) & mesh.Sf();
// store phi value before shift
DynamicList<scalar> faces_value_before;
faces_value_before.append(phi.internalField());
forAll(patchinds, ind_tmp){
    List<scalar> faces_value_curr = phi.boundaryField()[patchinds[ind_tmp]];
    if (faces_value_curr.size()==0){
        Info << "Empty boundary faces" << endl;
            List<scalar> zero_empty(xyz_onesidebdr_face_num[ind_tmp/2]);
            forAll(zero_empty, zero_ind){
                zero_empty[zero_ind] = 0;
            }
             faces_value_curr = zero_empty;
    }
    faces_value_before.append(faces_value_curr);
}

List<List<scalar> > faces_value_before_lst(num_proc);
faces_value_before_lst[proc_curr] = faces_value_before;
Pstream::gatherList(faces_value_before_lst);
Pstream::scatterList(faces_value_before_lst);


moveSurfaceScalarField(phi, num_proc, proc_curr, xyz_ind_distance, dimension_global, dimension_local, dimension_processor, patchinds, xyz_onesidebdr_face_num, internal_face_num, cell_face_dict, boundary_face_order_dict_lst);
// store phi value after shift
DynamicList<scalar> faces_value_after;
faces_value_after.append(phi.internalField());
forAll(patchinds, ind_tmp){
    List<scalar> faces_value_curr = phi.boundaryField()[patchinds[ind_tmp]];
    // this is to solve empty boundary conditions
    if (faces_value_curr.size()==0){
        Info << "Empty boundary faces" << endl;
            List<scalar> zero_empty(xyz_onesidebdr_face_num[ind_tmp/2]);
            forAll(zero_empty, zero_ind){
                zero_empty[zero_ind] = 0;
            }
             faces_value_curr = zero_empty;
    }
    faces_value_after.append(faces_value_curr);
}

List<List<scalar> > faces_value_after_lst(num_proc);
faces_value_after_lst[proc_curr] = faces_value_after;
Pstream::gatherList(faces_value_after_lst);
Pstream::scatterList(faces_value_after_lst);



// store face center value before shift
DynamicList<vector> faces_value;
faces_value.append(mesh.Cf());
forAll(patchinds, ind_tmp){
    List<vector> faces_value_curr = mesh.Cf().boundaryField()[patchinds[ind_tmp]];
    // this is to solve empty boundary conditions
    if (faces_value_curr.size()==0){
        Info << "Empty boundary faces" << endl;
            List<vector> zero_empty(xyz_onesidebdr_face_num[ind_tmp/2], vector(0, 0, 0));
            faces_value_curr = zero_empty;
    }
    faces_value.append(faces_value_curr);
}

List<List<vector> > faces_value_lst(num_proc);
faces_value_lst[proc_curr] = faces_value;
Pstream::gatherList(faces_value_lst);
Pstream::scatterList(faces_value_lst);


List<vector> faces_value_new(faces_value);
//Let's say, move face properties left 1 cell
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
            //Pout << faces_ind_curr_local << faces_ind_target_local << endl;
            
            faces_value_new[faces_ind_curr_local[0]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[0]];
            faces_value_new[faces_ind_curr_local[1]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[1]];
            faces_value_new[faces_ind_curr_local[2]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[2]];
            
            // Pout << "xyz_ind_curr_local" << xyz_ind_curr_local << "xyz_ind_curr_global" << xyz_ind_curr_global << "xyz_ind_target_global" << xyz_ind_target_global << "pxyz_ind_target_local" << pxyz_ind_target_local << endl
            //      << faces_value[faces_ind_curr_local[0]] << faces_value_new[faces_ind_curr_local[0]] << endl
            //      << faces_value[faces_ind_curr_local[1]] << faces_value_new[faces_ind_curr_local[1]] << endl
            //      << faces_value[faces_ind_curr_local[2]] << faces_value_new[faces_ind_curr_local[2]] << endl;
        }
    }
}

List<List<vector> > faces_value_new_lst(num_proc);
faces_value_new_lst[proc_curr] = faces_value_new;
Pstream::gatherList(faces_value_new_lst);
Pstream::scatterList(faces_value_new_lst);

fileName outputFile(runTime.path()/runTime.timeName()/"phi_comparsion");
OFstream os(outputFile);
forAll (faces_value, faces_value_ind){
    if (faces_value_ind==0){os << endl<< "Internal:" << endl;}
    if (faces_value_ind==internal_face_num){os << endl << "Left:" << endl;}
    if (faces_value_ind==internal_face_num + x_onesidebdr_face_num){os << endl << "Right:" << endl;}
    if (faces_value_ind==internal_face_num + x_onesidebdr_face_num * 2){os << endl<< "Front:" << endl;}
    if (faces_value_ind==internal_face_num + x_onesidebdr_face_num * 2 + y_onesidebdr_face_num){os << endl<< "Back:" << endl;}
    if (faces_value_ind==internal_face_num + x_onesidebdr_face_num * 2 + y_onesidebdr_face_num * 2){os << endl<< "Bottom:" << endl;}
    if (faces_value_ind==internal_face_num + x_onesidebdr_face_num * 2 + y_onesidebdr_face_num * 2 + z_onesidebdr_face_num){os << endl<< "Top:" << endl;}
    os << faces_value_before[faces_value_ind] << faces_value[faces_value_ind] << faces_value_after[faces_value_ind] << endl;
}

phi.write();

for (label xi_local = -1; xi_local < dimension_local[0]; xi_local++){
    for (label yi_local = -1; yi_local < dimension_local[1]; yi_local++){
        for (label zi_local = -1; zi_local < dimension_local[2]; zi_local++){
        //     if ((xi_local == -1 && yi_local == -1) || (yi_local == -1 && zi_local == -1) || (xi_local == -1 && zi_local == -1)){
        //         continue;
        //     }
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
            //Pout << faces_ind_curr_local << faces_ind_target_local << endl;
            bool if_pass = 1;
            for (label xyz_ind_tmp = 0; xyz_ind_tmp < 3; xyz_ind_tmp++){
                if (xi_local == -1 && xyz_ind_tmp != 0){continue;}
                if (yi_local == -1 && xyz_ind_tmp != 1){continue;}
                if (zi_local == -1 && xyz_ind_tmp != 2){continue;}

                if ((xyz_ind_curr_local[xyz_ind_tmp] == -1 && xyz_ind_target_local[xyz_ind_tmp] != -1) || (xyz_ind_curr_local[xyz_ind_tmp] != -1 && xyz_ind_target_local[xyz_ind_tmp] == -1)){
                    if_pass = (faces_value_after_lst[pxyz_ind_curr_local[0]][faces_ind_curr_local[xyz_ind_tmp]] == -faces_value_before_lst[pxyz_ind_target_local[0]][faces_ind_target_local[xyz_ind_tmp]]);
                }
                else{
                    if_pass = (faces_value_after_lst[pxyz_ind_curr_local[0]][faces_ind_curr_local[xyz_ind_tmp]] == faces_value_before_lst[pxyz_ind_target_local[0]][faces_ind_target_local[xyz_ind_tmp]]);
                }
                if (faces_value_new_lst[pxyz_ind_curr_local[0]][faces_ind_curr_local[xyz_ind_tmp]] != faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[xyz_ind_tmp]]){
                    if_pass = 0;
                }

                if (!if_pass){
                    FatalErrorInFunction
                    << xyz_ind_tmp << " "
                    << xi_local << " " << yi_local << " " << zi_local
                    << pxyz_ind_curr_local << pxyz_ind_target_local
                    << faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[xyz_ind_tmp]] << faces_value_new_lst[pxyz_ind_curr_local[0]][faces_ind_curr_local[xyz_ind_tmp]]
                    << faces_value_after_lst[pxyz_ind_curr_local[0]][faces_ind_curr_local[xyz_ind_tmp]] << " " << faces_value_before_lst[pxyz_ind_target_local[0]][faces_ind_target_local[xyz_ind_tmp]] 
                    << "The face coordinates move check failed"
                    << abort(FatalError);
                }
            }

            // Pout << "xyz_ind_curr_local" << xyz_ind_curr_local << "xyz_ind_curr_global" << xyz_ind_curr_global << "xyz_ind_target_global" << xyz_ind_target_global << "pxyz_ind_target_local" << pxyz_ind_target_local << endl
            //      << faces_value[faces_ind_curr_local[0]] << faces_value_new[faces_ind_curr_local[0]] << endl
            //      << faces_value[faces_ind_curr_local[1]] << faces_value_new[faces_ind_curr_local[1]] << endl
            //      << faces_value[faces_ind_curr_local[2]] << faces_value_new[faces_ind_curr_local[2]] << endl;
        }
    }
}





/*
List<vector> faces_value_new(faces_value);
//Let's say, move face properties left 1 cell
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
            //Pout << faces_ind_curr_local << faces_ind_target_local << endl;
            
            faces_value_new[faces_ind_curr_local[0]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[0]];
            faces_value_new[faces_ind_curr_local[1]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[1]];
            faces_value_new[faces_ind_curr_local[2]] = faces_value_lst[pxyz_ind_target_local[0]][faces_ind_target_local[2]];
            
            // Pout << "xyz_ind_curr_local" << xyz_ind_curr_local << "xyz_ind_curr_global" << xyz_ind_curr_global << "xyz_ind_target_global" << xyz_ind_target_global << "pxyz_ind_target_local" << pxyz_ind_target_local << endl
            //      << faces_value[faces_ind_curr_local[0]] << faces_value_new[faces_ind_curr_local[0]] << endl
            //      << faces_value[faces_ind_curr_local[1]] << faces_value_new[faces_ind_curr_local[1]] << endl
            //      << faces_value[faces_ind_curr_local[2]] << faces_value_new[faces_ind_curr_local[2]] << endl;
        }
    }
}
*/
//Let's do some validation work
vector proc_size = mesh.points()[(xidim_local+1)*(yidim_local+1)*(zidim_local+1)-1] - mesh.points()[0];
vector domain_size (dimension_processor[0] * proc_size[0], dimension_processor[1] * proc_size[1], dimension_processor[2] * proc_size[2]);
vector cell_size (proc_size[0]/xidim_local, proc_size[1]/yidim_local, proc_size[2]/zidim_local);

Info << "domain size: " << domain_size << endl;
Info << "cell size: " << cell_size << endl;

forAll(faces_value, ind_tmp){
    if (faces_value[ind_tmp] == vector(0, 0, 0)){
        continue;
    }
    vector xyz_distance (xyz_ind_distance[0] * cell_size[0], xyz_ind_distance[1] * cell_size[1], xyz_ind_distance[2] * cell_size[2]);
    vector faces_value_correct = faces_value[ind_tmp] + xyz_distance;
    for (label xyzind_tmp = 0; xyzind_tmp < 3; xyzind_tmp++){
        if (faces_value_correct[xyzind_tmp] < 0){
            double quo = (-faces_value_correct[xyzind_tmp])/cell_size[xyzind_tmp];
            if (quo - floor (quo) > 0.1){
                faces_value_correct[xyzind_tmp] = cell_size[xyzind_tmp]/2;
            }
            else{
                faces_value_correct[xyzind_tmp] = 0;
            }
        }
        if (faces_value_correct[xyzind_tmp] > domain_size[xyzind_tmp]){
            double quo = (faces_value_correct[xyzind_tmp] - domain_size[xyzind_tmp])/cell_size[xyzind_tmp];
            if (quo - floor (quo) > 0.1){
                faces_value_correct[xyzind_tmp] = domain_size[xyzind_tmp] - cell_size[xyzind_tmp]/2;
            }
            else{
                faces_value_correct[xyzind_tmp] = domain_size[xyzind_tmp];
            }
        }        
    }    
    
    if (faces_value_new[ind_tmp] != faces_value_correct){
            FatalErrorInFunction
            << "current face value: " << faces_value[ind_tmp] 
            << "calculated value: " << faces_value_new[ind_tmp] 
            << "correct face value: " << faces_value_correct
            << "The face coordinates move check failed"
            << abort(FatalError);
    }
    //Pout << faces_value[ind_tmp] << faces_value_new[ind_tmp] << faces_value_correct << (faces_value_new[ind_tmp] == faces_value_correct) << endl;
}

/*
DynamicList<scalar> phis_value;
phis_value.append(phi.internalField());
forAll(patchinds, ind_tmp){
    phis_value.append(phi.boundaryField()[patchinds[ind_tmp]]);
}
Info << phis_value << endl;
*/

