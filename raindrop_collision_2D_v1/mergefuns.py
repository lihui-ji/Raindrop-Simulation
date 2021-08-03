def mergeInit(d2_drop_path,d1_drop_path,m,n,m_new,n_new,threShold,p_rgh,hor_adjust):

    from mergefuns import downScale
    from mergefuns import pieceInfo
    from mergefuns import readFoamScalarOutput
    from mergefuns import mergeTwoPieces
    from mergefuns import readFoamVectorOutput
    from mergefuns import writeScalar
    from mergefuns import writeVector
    
    # Calculate transfer arrays
    aw_d1=readFoamScalarOutput(d1_drop_path+'/alpha.water')
    transfer_1=pieceInfo(aw_d1,threShold,m)
    aw_d2_raw=readFoamScalarOutput(d2_drop_path+'/alpha.water')
    aw_d2=downScale(aw_d2_raw,m,n,m_new,n_new,0)
    transfer_2=pieceInfo(aw_d2,threShold,m_new)

    # Merge alpha.water file
    aw_d1=readFoamScalarOutput(d1_drop_path+'/alpha.water')
    aw_d2_raw=readFoamScalarOutput(d2_drop_path+'/alpha.water')
    aw_d2=downScale(aw_d2_raw,m,n,m_new,n_new,0)
    merged_aw=mergeTwoPieces(aw_d2,aw_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)

    # Merge p_rgh file
    p_rgh_d1=readFoamScalarOutput(d1_drop_path+'/p_rgh')
    p_rgh_d2_raw=readFoamScalarOutput(d2_drop_path+'/p_rgh')
    p_rgh_d2=downScale(p_rgh_d2_raw,m,n,m_new,n_new,p_rgh)
    merged_p_rgh=mergeTwoPieces(p_rgh_d2,p_rgh_d1,m_new,n_new,transfer_2,transfer_1,p_rgh,hor_adjust)

    # Merge U
    [U0_d1,U1_d1,U2_d1]=readFoamVectorOutput(d1_drop_path+'/U')
    [U0_d2_raw,U1_d2_raw,U2_d2_raw]=readFoamVectorOutput(d2_drop_path+'/U')
    U0_d2=downScale(U0_d2_raw,m,n,m_new,n_new,0)
    U1_d2=downScale(U1_d2_raw,m,n,m_new,n_new,0)
    U2_d2=downScale(U2_d2_raw,m,n,m_new,n_new,0)
    merged_U0=mergeTwoPieces(U0_d2,U0_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)
    merged_U1=mergeTwoPieces(U1_d2,U1_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)
    merged_U2=mergeTwoPieces(U2_d2,U1_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)
    merged_U=[merged_U0,merged_U1,merged_U2]

    # Write files
    writeScalar('alpha.water',merged_aw,'0',m_new,n_new)
    writeScalar('p_rgh',merged_p_rgh,str(p_rgh),m_new,n_new)
    writeVector('U',merged_U,'(0 0 0)',m_new,n_new)

def downScale(old,m_old,n_old,m_new,n_new,default_value): # old is one dimentional array m_old*n_old
    # this is a simple downscaling method, which creates new points by averaging neighbors
    new=[default_value]*m_new*n_new
    for xi in range(m_old-1):
        for yi in range(n_old-1):
            new[xi*2+yi*2*m_new]=old[xi+yi*m_old]
            new[xi*2+1+yi*2*m_new]=0.5*(old[xi+yi*m_old]+old[xi+1+yi*m_old])
            new[xi*2+(yi*2+1)*m_new]=0.5*(old[xi+yi*m_old]+old[xi+(yi+1)*m_old])
            new[xi*2+1+(yi*2+1)*m_new]=0.25*(old[xi+yi*m_old]+old[xi+1+yi*m_old]+old[xi+(yi+1)*m_old]+old[xi+1+(yi+1)*m_old])
    for xi in range(m_new-1):
        new[xi+(n_new-1)*m_new]=new[xi+(n_new-2)*m_new]
    for yi in range(n_new):
        new[m_new-1+yi*m_new]=new[m_new-2+yi*m_new]
    return new

def pieceInfo(aw,threShold,m): # m is horizontal grid number of one dimentional array aw
    # get alpha.water piece information (slice center and half width)
    aw_xi=[]
    aw_zi=[]
    for ind in range(len(aw)):
        if aw[ind]>threShold:
            aw_xi.append(ind%m)
            aw_zi.append(int(ind/m))
    centerxi=int(sum(aw_xi)/len(aw_xi))
    centerzi=int(sum(aw_zi)/len(aw_zi))
    slice_half_width=min(centerxi,abs(centerxi-m))
    return [centerxi,centerzi,slice_half_width,m]

def readFoamScalarOutput(filePath): # read OpenFOAM scalar output files only
    with open(filePath, 'r') as file:
        data_temp = file.readlines()
    data=[]
    for point in data_temp[22:int(data_temp[20])+22]:
        data.append(float(point))
    return data

def readFoamVectorOutput(filePath): # read OpenFOAM scalar output files only
    with open(filePath, 'r') as file:
        data_temp = file.readlines()
    [U0,U1,U2]=[[],[],[]]
    for point in data_temp[22:int(data_temp[20])+22]:
        U0.append([float(a) for a in point[1:-2].split(' ')][0])
        U1.append([float(a) for a in point[1:-2].split(' ')][1])
        U2.append([float(a) for a in point[1:-2].split(' ')][2])
    return [U0,U1,U2]

def mergeTwoPieces(data_d2,data_d1,m_new,n_new,transfer_2,transfer_1,default_value,hor_adjust):
    block=[default_value]*m_new*n_new
    [centerxi_d1,centerzi_d1,slice_half_width_d1,m_d1]=transfer_1
    [centerxi_d2,centerzi_d2,slice_half_width_d2,m_d2]=transfer_2

    for xi_d2 in range(slice_half_width_d2*2+1):
        for zi_d2 in range(int(n_new/2)):
            block_xi=int(xi_d2+m_new/2-slice_half_width_d2)
            block_zi=int(zi_d2+n_new/2)
            slice_xi=int(centerxi_d2-slice_half_width_d2+xi_d2)
            slice_zi=int(centerzi_d2-n_new/16+zi_d2)
            block[block_xi+block_zi*m_new]=data_d2[slice_zi*m_d2+slice_xi]
        
    for xi_d1 in range(slice_half_width_d1*2+1):
        for zi_d1 in range(-int(n_new/4),int(min(n_new/4,n_new/2-centerzi_d1+n_new/16))):
            block_xi=int(xi_d1+m_new/2-slice_half_width_d1+hor_adjust)
            block_zi=int(zi_d1+n_new/4)
            slice_xi=int(centerxi_d1-slice_half_width_d1+xi_d1)
            slice_zi=int(centerzi_d1-n_new/16+zi_d1)
            block[block_xi+block_zi*m_new]=data_d1[slice_zi*m_d1+slice_xi] 
    return block

def writeScalar(fileName,merged_data,default_value,m,n):
    import os
    with open('0/'+fileName+'.orig', 'r') as file:
        data_write = file.read()
    data_tem='internalField   nonuniform List<scalar>\n'+str(m*n)+'\n(\n'
    for i in range(m*n):
        data_tem=data_tem+str(merged_data[i])+'\n'
    data_tem=data_tem+')\n;\n'
    data_write=data_write.replace('internalField   uniform '+default_value+';\n',data_tem)
    os.system("touch 0/temp.txt")
    with open('0/temp.txt', 'w') as file:
        file.write(data_write)
    os.system("cp 0/temp.txt 0/"+str(fileName))
    os.system("rm 0/temp.txt")

def writeVector(fileName,merged_data,default_value,m,n):
    import os
    with open('0/'+fileName+'.orig', 'r') as file:
        data_write = file.read()
    data_tem='internalField   nonuniform List<vector>\n'+str(m*n)+'\n(\n'
    for i in range(m*n):
        data_tem=data_tem+'('+str(merged_data[0][i])+' '+str(merged_data[1][i])+' '+str(merged_data[2][i])+')\n'
    data_tem=data_tem+')\n;\n'
    data_write=data_write.replace('internalField   uniform '+default_value+';\n',data_tem)
    os.system("touch 0/temp.txt")
    filepath='0/temp.txt'
    with open(filepath, 'w') as file:
        file.write(data_write)
    os.system("cp 0/temp.txt 0/U")
    os.system("rm 0/temp.txt")