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

DynamicList<vector> faces_value;
faces_value.append(mesh.Cf());
forAll(patchinds, ind_tmp){
    List<vector> faces_value_curr = mesh.Cf().boundaryField()[patchinds[ind_tmp]];
    // this is to solve empty boundary conditions
    if (faces_value_curr.size()!=0){
        Info << "Empty boundary faces" << endl;
        faces_value.append(faces_value_curr);
    }
    else{
        List<vector> zero_empty(xyz_onesidebdr_face_num[ind_tmp/2], vector(0, 0, 0));
        faces_value.append(zero_empty);
    }

    Info << mesh.Cf().boundaryField()[patchinds[ind_tmp]].size() << " ";
}
Info << mesh.Cf() << endl;
Info << faces_value << endl;


label zi_face_num = xidim_local*yidim_local+(xidim_local-1)*yidim_local+(yidim_local-1)*xidim_local;
List<List<label>> cell_face_dict(xidim_local*yidim_local*zidim_local);
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
            List<label> face_inds_curr(3, -1);
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
Info << cell_face_dict << endl;

forAll(cell_face_dict, cell_ind_tmp){
    Pout << faces_value[cell_face_dict[cell_ind_tmp][0]] << " " << faces_value[cell_face_dict[cell_ind_tmp][1]] << " " << faces_value[cell_face_dict[cell_ind_tmp][2]] << endl;
}



DynamicList<scalar> phis_value;
phis_value.append(phi.internalField());
forAll(patchinds, ind_tmp){
    phis_value.append(phi.boundaryField()[patchinds[ind_tmp]]);
}
Info << phis_value << endl;


