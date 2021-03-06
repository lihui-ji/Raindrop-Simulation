scalar threshold = 0.5;
List<DynamicList<Vector<label>>> phase_cell_ind_lst = 
                    getPhaseCellIndLst(alpha1, dimension_local, dimension_processor, num_proc, proc_curr, threshold);
Vector<label> ind_phase_center = getIndCenter(phase_cell_ind_lst);
phaseInfo (ind_phase_center, phase_cell_ind_lst, cell_center_lst, cell_volume_lst, dimension_local, dimension_processor, threshold);
Vector<label> xyz_ind_distance = ind_phase_center - ind_center_global;
moveInfo (xyz_ind_distance);
if (xyz_ind_distance[0]!=0 || xyz_ind_distance[1]!=0 || xyz_ind_distance[2]!=0){
    moveTypeFields<scalar>(mesh, num_proc, proc_curr, xyz_ind_distance, dimension_global, dimension_local, dimension_processor);
    moveTypeFields<vector>(mesh, num_proc, proc_curr, xyz_ind_distance, dimension_global, dimension_local, dimension_processor);
    correctAllBoundaryConditions(mesh);
    moveSurfaceScalarField(phi, num_proc, proc_curr, xyz_ind_distance, dimension_global, dimension_local, dimension_processor, patchinds, xyz_onesidebdr_face_num, internal_face_num, cell_face_dict, boundary_face_order_dict_lst);
    moveSurfaceScalarField(rhoPhi, num_proc, proc_curr, xyz_ind_distance, dimension_global, dimension_local, dimension_processor, patchinds, xyz_onesidebdr_face_num, internal_face_num, cell_face_dict, boundary_face_order_dict_lst);
    moveSurfaceScalarField(alphaPhi10, num_proc, proc_curr, xyz_ind_distance, dimension_global, dimension_local, dimension_processor, patchinds, xyz_onesidebdr_face_num, internal_face_num, cell_face_dict, boundary_face_order_dict_lst);
}
//phi = fvc::interpolate(U) & mesh.Sf();


