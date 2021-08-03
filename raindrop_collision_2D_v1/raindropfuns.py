## OpenFOAM setting functions
def writeblockMeshDict(domainwidth,m_grid,n_grid):
    import os
    os.system('cp system/blockMeshDict system/blockMeshDict.txt')    
    with open('system/blockMeshDict.txt', 'r') as file:
        data = file.readlines()
    data[16]='convertToMeters '+str(domainwidth)+';\n'
    data[35]='    hex (0 1 2 3 4 5 6 7) ('+str(m_grid)+' '+str(n_grid)+' 1) simpleGrading (1 1 1)\n'
    with open('system/blockMeshDict.txt', 'w') as file:
        file.writelines(data)
    os.system('cp system/blockMeshDict.txt system/blockMeshDict')
    os.system('rm system/blockMeshDict.txt')

def writesetFieldsDict(dropdiameter,dropcenterx,dropcentery):
    import os
    dropradius=float(dropdiameter)/2
    os.system('cp system/setFieldsDict system/setFieldsDict.txt')    
    with open('system/setFieldsDict.txt', 'r') as file:
        data = file.readlines()
    data[26]='        p1 ('+str(dropcenterx)+' '+str(dropcentery)+' -1); p2 ('+str(dropcenterx)+' '+str(dropcentery)+' 1); radius '+str(dropradius)+';\n'
    with open('system/setFieldsDict.txt', 'w') as file:
        file.writelines(data)
    os.system('cp system/setFieldsDict.txt system/setFieldsDict')
    os.system('rm system/setFieldsDict.txt')

def writecontrolDict(writeInterval,t0,t1):
    import os
    with open('system/controlDict.txt', 'r') as file:
        data = file.readlines()
    data[31]='writeInterval       '+str(writeInterval)+';\n' 
    data[21]='startTime       '+str(t0)+';\n' 
    data[25]='endTime         '+str(t1)+';\n'
    with open('system/controlDict.txt', 'w') as file:
        file.writelines(data)
    os.system('cp system/controlDict.txt system/controlDict')

## OpenFOAM input writing functions

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

## Sliding mesh functions

def prep_processor01(domainwidth,m_grid,n_grid):
    import os
    writeblockMeshDict(domainwidth,m_grid,n_grid)
    writecontrolDict(0.0001,str(0),str(0.0001))
    os.system('blockMesh > log.blockMesh')
    os.system('decomposePar > log.decomposePar')
    os.system('./runParallel')
    os.system('mkdir prepp0')
    os.system('mkdir prepp1')
    os.system('cp -r processor0/0.0001 prepp0')
    os.system('cp -r processor1/0.0001 prepp1')
    os.system('rm log*')
    os.system('rm -r processor*')


    os.system('cp system/blockMeshDict system/blockMeshDict.txt')    
    with open('system/blockMeshDict.txt', 'r') as file:
        data = file.readlines()
    data[16]='convertToMeters '+str(domainwidth)+';\n'
    data[35]='    hex (0 1 2 3 4 5 6 7) ('+str(m_grid)+' '+str(n_grid)+' 1) simpleGrading (1 1 1)\n'
    with open('system/blockMeshDict.txt', 'w') as file:
        file.writelines(data)
    os.system('cp system/blockMeshDict.txt system/blockMeshDict')
    os.system('rm system/blockMeshDict.txt')

def slidevertical(writeInterval,keyword,pressure):
    import os
    os.system('cp processor0/'+keyword+'/U processor2/'+keyword+'/U.txt')
    os.system('cp processor0/'+keyword+'/alpha.water processor2/'+keyword+'/aw.txt')
    #os.system('cp processor0/'+keyword+'/alphaPhi0.water processor2/'+keyword+'/ap0w.txt')
    os.system('cp processor0/'+keyword+'/p processor2/'+keyword+'/p.txt')
    os.system('cp processor0/'+keyword+'/p_rgh processor2/'+keyword+'/p_rgh.txt')
    #os.system('cp processor0/'+keyword+'/phi processor2/'+keyword+'/phi.txt')

    os.system('cp processor1/'+keyword+'/U processor3/'+keyword+'/U.txt')
    os.system('cp processor1/'+keyword+'/alpha.water processor3/'+keyword+'/aw.txt')
    #os.system('cp processor1/'+keyword+'/alphaPhi0.water processor3/'+keyword+'/ap0w.txt')
    os.system('cp processor1/'+keyword+'/p processor3/'+keyword+'/p.txt')
    os.system('cp processor1/'+keyword+'/p_rgh processor3/'+keyword+'/p_rgh.txt')
    #os.system('cp processor1/'+keyword+'/phi processor3/'+keyword+'/phi.txt')

    os.system('cp prepp0/0.0001'+'/alpha.water processor0/'+keyword)
    #os.system('cp prepp0/0.0001'+'/alphaPhi0.water processor0/'+keyword)
    os.system('cp prepp0/0.0001'+'/p processor0/'+keyword)
    os.system('cp prepp0/0.0001'+'/p_rgh processor0/'+keyword)
    #os.system('cp prepp0/0.0001'+'/phi processor0/'+keyword)
    os.system('cp prepp0/0.0001'+'/U processor0/'+keyword)

    os.system('cp prepp1/0.0001'+'/alpha.water processor1/'+keyword)
    #os.system('cp prepp1/0.0001'+'/alphaPhi0.water processor1/'+keyword)
    os.system('cp prepp1/0.0001'+'/p processor1/'+keyword)
    os.system('cp prepp1/0.0001'+'/p_rgh processor1/'+keyword)
    #os.system('cp prepp1/0.0001'+'/phi processor1/'+keyword)
    os.system('cp prepp1/0.0001'+'/U processor1/'+keyword)
    

    with open('processor3/'+keyword+'/U.txt', 'r') as file:
            data = file.read()
    data=data.replace('    top\n    {\n        type            pressureInletOutletVelocity;\n        value           nonuniform List<vector> 0();\n    }\n','')
    data=data.replace('procBoundary1to0','procBoundary3to2')
    data=data.replace('    procBoundary1to2\n    {\n        type            processor;\n','    top\n    {\n        type            pressureInletOutletVelocity;\n')    
    with open('processor3/'+keyword+'/U.txt', 'w') as file:
            file.write(data)

    with open('processor3/'+keyword+'/aw.txt', 'r') as file:
            data = file.read()
    data=data.replace('    top\n    {\n        type            inletOutlet;\n        inletValue      nonuniform List<scalar> 0();\n        value           nonuniform List<scalar> 0();\n    }\n','')
    data=data.replace('procBoundary1to0','procBoundary3to2')
    data=data.replace('    procBoundary1to2\n    {\n        type            processor;\n','    top\n    {\n        type            inletOutlet;\n        inletValue      uniform 0;\n')    
    with open('processor3/'+keyword+'/aw.txt', 'w') as file:
            file.write(data)

    #with open('processor3/'+keyword+'/ap0w.txt', 'r') as file:
    #        data = file.read()
    #data=data.replace('    top\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n','')
    #data=data.replace('procBoundary1to0','procBoundary3to2')
    #data=data.replace('    procBoundary1to2\n    {\n        type            processor;\n','    top\n    {\n        type            calculated;\n')    
    #with open('processor3/'+keyword+'/ap0w.txt', 'w') as file:
    #        file.write(data)

    with open('processor3/'+keyword+'/p.txt', 'r') as file:
            data = file.read()
    data=data.replace('    top\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n','')
    data=data.replace('procBoundary1to0','procBoundary3to2')
    data=data.replace('    procBoundary1to2\n    {\n        type            processor;\n','    top\n    {\n        type            calculated;\n')    
    with open('processor3/'+keyword+'/p.txt', 'w') as file:
            file.write(data)

    with open('processor3/'+keyword+'/p_rgh.txt', 'r') as file:
            data = file.read()
    data=data.replace('    top\n    {\n        type            totalPressure;\n        rho             rho;\n        psi             none;\n        gamma           1;\n        p0              nonuniform List<scalar> 0();\n        value           nonuniform List<scalar> 0();\n    }\n','')
    data=data.replace('procBoundary1to0','procBoundary3to2')
    data=data.replace('    procBoundary1to2\n    {\n        type            processor;\n','    top\n    {\n        type            totalPressure;\n        rho             rho;\n        psi             none;\n        gamma           1;\n        p0              uniform '+str(pressure)+';\n')    
    with open('processor3/'+keyword+'/p_rgh.txt', 'w') as file:
            file.write(data)

    #with open('processor3/'+keyword+'/phi.txt', 'r') as file:
    #        data = file.read()
    #data=data.replace('    top\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n','')
    #data=data.replace('procBoundary1to0','procBoundary3to2')
    #data=data.replace('    procBoundary1to2\n    {\n        type            processor;\n','    top\n    {\n        type            calculated;\n')    
    #with open('processor3/'+keyword+'/phi.txt', 'w') as file:
    #        file.write(data)

    with open('processor2/'+keyword+'/U.txt', 'r') as file:
            data = file.read()
    data=data.replace('procBoundary0to1','procBoundary2to3')
    data=data.replace('    bot\n    {\n        type            pressureInletOutletVelocity;\n','    bot\n    {\n        type            pressureInletOutletVelocity;\n        value           nonuniform List<vector> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n')    
    with open('processor2/'+keyword+'/U.txt', 'w') as file:
            file.write(data)

    with open('processor2/'+keyword+'/aw.txt', 'r') as file:
            data = file.read()
    data=data.replace('procBoundary0to1','procBoundary2to3')
    data=data.replace('    bot\n    {\n        type            inletOutlet;\n        inletValue      uniform 0;\n','    bot\n    {\n        type            inletOutlet;\n        inletValue      nonuniform List<scalar> 0();\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n')    
    with open('processor2/'+keyword+'/aw.txt', 'w') as file:
            file.write(data)

    #with open('processor2/'+keyword+'/ap0w.txt', 'r') as file:
    #        data = file.read()
    #data=data.replace('procBoundary0to1','procBoundary2to3')
    #data=data.replace('    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> ','    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n        value           nonuniform List<scalar> ')     
    #data=data.replace('    bot\n    {\n        type            calculated;\n        value           uniform 0;\n','    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n        value           uniform 0;\n')
    #with open('processor2/'+keyword+'/ap0w.txt', 'w') as file:
    #        file.write(data)

    with open('processor2/'+keyword+'/p.txt', 'r') as file:
            data = file.read()
    data=data.replace('procBoundary0to1','procBoundary2to3')
    data=data.replace('    bot\n    {\n        type            calculated;\n','    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n')    
    with open('processor2/'+keyword+'/p.txt', 'w') as file:
            file.write(data)

    with open('processor2/'+keyword+'/p_rgh.txt', 'r') as file:
            data = file.read()
    data=data.replace('procBoundary0to1','procBoundary2to3')
    data=data.replace('    bot\n    {\n        type            totalPressure;\n        rho             rho;\n        psi             none;\n        gamma           1;\n        p0              uniform '+str(pressure)+';\n        value           uniform '+str(pressure)+';\n    }\n','    bot\n    {\n        type            totalPressure;\n        rho             rho;\n        psi             none;\n        gamma           1;\n        p0              nonuniform List<scalar> 0();\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n        value           uniform '+str(pressure)+';\n    }\n')    
    data=data.replace('    bot\n    {\n        type            totalPressure;\n        rho             rho;\n        psi             none;\n        gamma           1;\n        p0              uniform '+str(pressure)+';\n        value           nonuniform List<scalar> \n','    bot\n    {\n        type            totalPressure;\n        rho             rho;\n        psi             none;\n        gamma           1;\n        p0              nonuniform List<scalar> 0();\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n        value           nonuniform List<scalar> \n')    
    with open('processor2/'+keyword+'/p_rgh.txt', 'w') as file:
            file.write(data)

    #with open('processor2/'+keyword+'/phi.txt', 'r') as file:
    #        data = file.read()
    #data=data.replace('procBoundary0to1','procBoundary2to3')
    #data=data.replace('    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> \n','    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n        value           nonuniform List<scalar> \n')    
    #with open('processor2/'+keyword+'/phi.txt', 'w') as file:
    #        file.write(data)

    os.system('cp processor3/'+keyword+'/U.txt processor3/'+keyword+'/U')
    os.system('cp processor3/'+keyword+'/aw.txt processor3/'+keyword+'/alpha.water')
    #os.system('cp processor3/'+keyword+'/ap0w.txt processor3/'+keyword+'/alphaPhi0.water')
    os.system('cp processor3/'+keyword+'/p.txt processor3/'+keyword+'/p')
    os.system('cp processor3/'+keyword+'/p_rgh.txt processor3/'+keyword+'/p_rgh')
    #os.system('cp processor3/'+keyword+'/phi.txt processor3/'+keyword+'/phi')

    os.system('cp processor2/'+keyword+'/U.txt processor2/'+keyword+'/U')
    os.system('cp processor2/'+keyword+'/aw.txt processor2/'+keyword+'/alpha.water')
    #os.system('cp processor2/'+keyword+'/ap0w.txt processor2/'+keyword+'/alphaPhi0.water')
    os.system('cp processor2/'+keyword+'/p.txt processor2/'+keyword+'/p')
    os.system('cp processor2/'+keyword+'/p_rgh.txt processor2/'+keyword+'/p_rgh')
    #os.system('cp processor2/'+keyword+'/phi.txt processor2/'+keyword+'/phi')

    os.system('rm processor3/'+keyword+'/*.txt')
    os.system('rm processor2/'+keyword+'/*.txt')

   
    # dist=[0, 0, 0, 0]
    # maxaw=[0, 0, 0, 0]
    # maxnx=[0, 0, 0, 0]
    # maxny=[0, 0, 0, 0]
    # for i in range(4):
    #     with open('processor'+str(i)+'/'+keyword+'/alpha.water', 'r') as file:
    #         data = file.readlines()
    #     if any(char.isdigit() for char in data[20]):
    #         floatdata=[]
    #         for ele in data[22:int(data[20])+22]:
    #             floatdata.append(float(ele))
    #         dist[i]=int(floatdata.index(max(floatdata))%m-(m/2))
    #         maxnx[i]=int(floatdata.index(max(floatdata))%m)
    #         maxny[i]=int(floatdata.index(max(floatdata))/m)+pn*i
    #         maxaw[i]=max(floatdata)
    # maxindex=maxaw.index(max(maxaw))
    # return [dist[maxindex],maxnx[maxindex],maxny[maxindex],maxaw[maxindex]]         

def slidehorizontal(keyword,dist,m,pn,pressure):
    import os
    for i in range(4):
        os.system('cp processor'+str(i)+'/'+keyword+'/U processor'+str(i)+'/'+keyword+'/U.txt')
        os.system('cp processor'+str(i)+'/'+keyword+'/alpha.water processor'+str(i)+'/'+keyword+'/aw.txt')
        os.system('cp processor'+str(i)+'/'+keyword+'/p processor'+str(i)+'/'+keyword+'/p.txt')
        os.system('cp processor'+str(i)+'/'+keyword+'/p_rgh processor'+str(i)+'/'+keyword+'/p_rgh.txt')

        with open('processor'+str(i)+'/'+keyword+'/aw.txt', 'r') as file:
            data = file.readlines()
        if any(char.isdigit() for char in data[20]):
            for j in range(pn):
                if dist>0:                
                    data[22+m*j:22+m*j+m-dist]=data[22+m*j+dist:22+m*j+m]
                    data[22+m*j+m-dist:22+m*j+m]=['0\n']*dist
                if dist<0:
                    data[22+m*j-dist:22+m*j+m]=data[22+m*j:22+m*j+m+dist]
                    data[22+m*j:22+m*j-dist]=['0\n']*(-dist)            
            for startindex in [startindexi+2 for startindexi, ele in enumerate(data) if ele == str(m)+'\n']:
                if dist>0:
                    data[startindex:startindex+m-dist]=data[startindex+dist:startindex+m]
                    data[startindex+m-dist:startindex+m]=['0\n']*dist
                if dist<0:
                    data[startindex-dist:startindex+m]=data[startindex:startindex+m+dist]
                    data[startindex:startindex-dist]=['0\n']*(-dist)                   
            with open('processor'+str(i)+'/'+keyword+'/aw.txt', 'w') as file:
                file.writelines(data)
            os.system('cp processor'+str(i)+'/'+keyword+'/aw.txt processor'+str(i)+'/'+keyword+'/alpha.water')
 
        with open('processor'+str(i)+'/'+keyword+'/U.txt', 'r') as file:
            data = file.readlines()         
        if any(char.isdigit() for char in data[20]):
            for j in range(pn):
                if dist>0:                
                    data[22+m*j:22+m*j+m-dist]=data[22+m*j+dist:22+m*j+m]
                    data[22+m*j+m-dist:22+m*j+m]=['(0 0 0)\n']*dist
                if dist<0:
                    data[22+m*j-dist:22+m*j+m]=data[22+m*j:22+m*j+m+dist]
                    data[22+m*j:22+m*j-dist]=['(0 0 0)\n']*(-dist)            
            for startindex in [startindexi+2 for startindexi, ele in enumerate(data) if ele == str(m)+'\n']:
                if dist>0:
                    data[startindex:startindex+m-dist]=data[startindex+dist:startindex+m]
                    data[startindex+m-dist:startindex+m]=['(0 0 0)\n']*dist
                if dist<0:
                    data[startindex-dist:startindex+m]=data[startindex:startindex+m+dist]
                    data[startindex:startindex-dist]=['(0 0 0)\n']*(-dist)        
            with open('processor'+str(i)+'/'+keyword+'/U.txt', 'w') as file:
                file.writelines(data)
            os.system('cp processor'+str(i)+'/'+keyword+'/U.txt processor'+str(i)+'/'+keyword+'/U')

        with open('processor'+str(i)+'/'+keyword+'/p.txt', 'r') as file:
            data = file.readlines()
        if any(char.isdigit() for char in data[20]):
            for j in range(pn):
                if dist>0:                
                    data[22+m*j:22+m*j+m-dist]=data[22+m*j+dist:22+m*j+m]
                    data[22+m*j+m-dist:22+m*j+m]=[str(pressure)+'\n']*dist
                if dist<0:
                    data[22+m*j-dist:22+m*j+m]=data[22+m*j:22+m*j+m+dist]
                    data[22+m*j:22+m*j-dist]=[str(pressure)+'\n']*(-dist)            
            for startindex in [startindexi+2 for startindexi, ele in enumerate(data) if ele == str(m)+'\n']:
                if dist>0:
                    data[startindex:startindex+m-dist]=data[startindex+dist:startindex+m]
                    data[startindex+m-dist:startindex+m]=[str(pressure)+'\n']*dist
                if dist<0:
                    data[startindex-dist:startindex+m]=data[startindex:startindex+m+dist]
                    data[startindex:startindex-dist]=[str(pressure)+'\n']*(-dist)        
            with open('processor'+str(i)+'/'+keyword+'/p.txt', 'w') as file:
                file.writelines(data)
            os.system('cp processor'+str(i)+'/'+keyword+'/p.txt processor'+str(i)+'/'+keyword+'/p')

        with open('processor'+str(i)+'/'+keyword+'/p_rgh.txt', 'r') as file:
            data = file.readlines()
        if any(char.isdigit() for char in data[20]):
            for j in range(pn):
                if dist>0:                
                    data[22+m*j:22+m*j+m-dist]=data[22+m*j+dist:22+m*j+m]
                    data[22+m*j+m-dist:22+m*j+m]=[str(pressure)+'\n']*dist
                if dist<0:
                    data[22+m*j-dist:22+m*j+m]=data[22+m*j:22+m*j+m+dist]
                    data[22+m*j:22+m*j-dist]=[str(pressure)+'\n']*(-dist)            
            for startindex in [startindexi+2 for startindexi, ele in enumerate(data) if ele == str(m)+'\n']:
                if dist>0:
                    data[startindex:startindex+m-dist]=data[startindex+dist:startindex+m]
                    data[startindex+m-dist:startindex+m]=[str(pressure)+'\n']*dist
                if dist<0:
                    data[startindex-dist:startindex+m]=data[startindex:startindex+m+dist]
                    data[startindex:startindex-dist]=[str(pressure)+'\n']*(-dist)        
        with open('processor'+str(i)+'/'+keyword+'/p_rgh.txt', 'w') as file:
            file.writelines(data)
        os.system('cp processor'+str(i)+'/'+keyword+'/p_rgh.txt processor'+str(i)+'/'+keyword+'/p_rgh')

    #os.system('rm processor*/'+keyword+'/*.txt')      

## OpenFOAM output reading functions

def dropinfo(keyword,m,pn,threshold):
    import math
    import os
    tempaw=[]
    tempnx=[]
    tempny=[]
    tempUx=[]
    tempUz=[]
    temp=[]
    for pi in range(4):
        with open('processor'+str(pi)+'/'+keyword+'/alpha.water', 'r') as file:
            awdata = file.readlines()
        with open('processor'+str(pi)+'/'+keyword+'/U', 'r') as file:
            Udata = file.readlines()
        if any(char.isdigit() for char in awdata[20]):
            for awi in range(int(awdata[20])):
                if float(awdata[awi+22])>=threshold:
                    tempaw.append(float(awdata[awi+22]))
                    tempnx.append(awi%m)
                    tempny.append(int(awi/m)+pn*pi)
                    temp.append([awi%m,int(awi/m)])
                    tempUx.append(float((Udata[awi+22])[1:-2].split(' ')[0]))
                    tempUz.append(float((Udata[awi+22])[1:-2].split(' ')[1]))
    if len(tempnx)==0:
        os.system("reconstructPar")
        os.system("foamToVTK")
        exit()

    centernx=sum(tempnx)/len(tempnx)
    centerny=sum(tempny)/len(tempny)
    averageaw=sum(tempaw)/len(tempaw)
    averageUx=sum(tempUx)/len(tempUx)
    averageUz=sum(tempUz)/len(tempUz)

    group=[]
    while len(temp)!=0:
        subgroup=[]
        subgroup.append(temp[0])
        temp.remove(temp[0])

        subni=0
        
        while subni!=len(subgroup):
            curr=subgroup[subni]
            neighbors=[[curr[0]-1,curr[1]], [curr[0]+1,curr[1]], [curr[0],curr[1]-1], [curr[0],curr[1]+1]]
            for neighbor in neighbors:
                if neighbor in temp:
                    subgroup.append(neighbor)
                    temp.remove(neighbor)
            subni=subni+1 
        group.append(subgroup)
    number = len(group)
    if number == 2:
        nx0, ny0=[loc[0] for loc in group[0]],[loc[1] for loc in group[0]]
        nx1, ny1=[loc[0] for loc in group[1]],[loc[1] for loc in group[1]]
        centernx0, centerny0 = float(sum(nx0))/len(nx0), sum(ny0)/len(ny0)
        centernx1, centerny1 = float(sum(nx1))/len(nx1), sum(ny1)/len(ny1)
        angle = math.atan(abs((centernx0-centernx1)/(centerny0-centerny1)))/math.pi*180
    else:
        angle =-1

    datatowrite='nx ny alpha.water Ux(m/s) Uz(m/s)\n'
    for li in range(len(tempnx)):
        datatowrite+=str(tempnx[li])+' '+str(tempny[li])+' '+str(tempaw[li])+' '+str(tempUx[li])+' '+str(tempUz[li])+'\n'
    os.system('touch droploc/drop'+keyword)

    with open('droploc/drop'+keyword,'w') as file:
        file.write(datatowrite)

    return [centernx,centerny,averageaw,averageUx,averageUz,number,angle] 
  
def readFoamScalarOutput(filePath): # read OpenFOAM scalar output files only
    import os
    with open(filePath, 'r') as file:
        data_temp = file.readlines()
    data=[]
    for point in data_temp[22:int(data_temp[20])+22]:
        data.append(float(point))
    return data

def readFoamVectorOutput(filePath): # read OpenFOAM scalar output files only
    import os
    with open(filePath, 'r') as file:
        data_temp = file.readlines()
    [U0,U1,U2]=[[],[],[]]
    for point in data_temp[22:int(data_temp[20])+22]:
        U0.append([float(a) for a in point[1:-2].split(' ')][0])
        U1.append([float(a) for a in point[1:-2].split(' ')][1])
        U2.append([float(a) for a in point[1:-2].split(' ')][2])
    return [U0,U1,U2] 

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

def dropab(filePath,m,n,domainwidth,threShold):
    aw=readFoamScalarOutput(filePath+'/alpha.water')
    drop_mi=[]
    drop_ni=[]
    for m_i in range(m):
        for n_i in range(n):
            if aw[m_i+n_i*m]>threShold:
                drop_mi.append(m_i)
                drop_ni.append(n_i)
    a=float(max(drop_mi)-min(drop_mi))/float(m)*domainwidth
    b=float(max(drop_ni)-min(drop_ni))/float(m)*domainwidth
    return [a,b]

## Data processing functions

def downScale(old,m_old,n_old,downScaleFactor,default_value): # old is one dimentional array m_old*n_old
    # this is a simple downscaling method, which creates new points by averaging neighbors
    import math
    m_new=int(m_old*downScaleFactor)
    n_new=int(n_old*downScaleFactor)
    new=[default_value]*m_new*n_new

    for xi_new in range(m_new-1):
        for yi_new in range(n_new-1):
            x = xi_new / downScaleFactor
            y = yi_new / downScaleFactor
            [x1, x2, y1, y2] = [math.floor(x), math.ceil(x), math.floor(y), math.ceil(y)]
            [Q11, Q12, Q21, Q22] = [old[int(x1+y1*m_old)], old[int(x1+y2*m_old)], old[int(x2+y1*m_old)], old[int(x2+y2*m_old)]]
            if x1 == x2:
                [Q_y1, Q_y2] = [Q11, Q22]
            else:
                [Q_y1, Q_y2] = [(x2-x)*Q11/(x2-x1) + (x-x1)*Q21/(x2-x1) , (x2-x)*Q12/(x2-x1) + (x-x1)*Q22/(x2-x1)]
            if y1 == y2:
                Q = Q_y1
            else:
                Q = (y2-y)*Q_y1/(y2-y1) + (y-y1)*Q_y2/(y2-y1)
            new[xi_new + m_new*yi_new] = Q
     
    #for xi in range(m_old-1):
    #    for yi in range(n_old-1):
    #        new[xi*2+yi*2*m_new]=old[xi+yi*m_old]
    #        new[xi*2+1+yi*2*m_new]=0.5*(old[xi+yi*m_old]+old[xi+1+yi*m_old])
    #        new[xi*2+(yi*2+1)*m_new]=0.5*(old[xi+yi*m_old]+old[xi+(yi+1)*m_old])
    #        new[xi*2+1+(yi*2+1)*m_new]=0.25*(old[xi+yi*m_old]+old[xi+1+yi*m_old]+old[xi+(yi+1)*m_old]+old[xi+1+(yi+1)*m_old])
    #for xi in range(m_new-1):
    #    new[xi+(n_new-1)*m_new]=new[xi+(n_new-2)*m_new]
    #for yi in range(n_new):
    #    new[m_new-1+yi*m_new]=new[m_new-2+yi*m_new]
    return new

## Merging functions, which combine two output profiles to one for collision analysis

def mergeInit(d2_drop_path,d1_drop_path,m_2,m_1,m_merged,downScaleFactor,threShold,p_rgh,hor_adjust,warmperiod,reverse):
    import os.path
    n_1=m_1*2
    n_2=m_2*2
    m_new=m_merged
    n_new=m_new*2
    n_merged=m_merged*2    
    # Calculate transfer arrays
    aw_d1=readFoamScalarOutput(d1_drop_path+'/alpha.water')    
    aw_d2_raw=readFoamScalarOutput(d2_drop_path+'/alpha.water')
    aw_d2=downScale(aw_d2_raw,m_2,n_2,downScaleFactor,0)
    transfer_1=pieceInfo(aw_d1,threShold,m_1)
    transfer_2=pieceInfo(aw_d2,threShold,int(m_2*downScaleFactor))
    if reverse==True:
        [transfer_1,transfer_2]=[transfer_2,transfer_1]
    # alpha.water file
    aw_d1=readFoamScalarOutput(d1_drop_path+'/alpha.water')
    aw_d2_raw=readFoamScalarOutput(d2_drop_path+'/alpha.water')
    aw_d2=downScale(aw_d2_raw,m_2,n_2,downScaleFactor,0)
    if reverse==True:
        [aw_d1, aw_d2] = [aw_d2, aw_d1]
    # p_rgh file
    p_rgh_d1=readFoamScalarOutput(d1_drop_path+'/p_rgh')
    p_rgh_d2_raw=readFoamScalarOutput(d2_drop_path+'/p_rgh')
    p_rgh_d2=downScale(p_rgh_d2_raw,m_2,n_2,downScaleFactor,p_rgh)
    if reverse==True:
        [p_rgh_d1, p_rgh_d2] = [p_rgh_d2, p_rgh_d1]    

    # U
    [U0_d1,U1_d1,U2_d1]=readFoamVectorOutput(d1_drop_path+'/U')
    [U0_d2_raw,U1_d2_raw,U2_d2_raw]=readFoamVectorOutput(d2_drop_path+'/U')
    U0_d2=downScale(U0_d2_raw,m_2,n_2,downScaleFactor,0)
    U1_d2=downScale(U1_d2_raw,m_2,n_2,downScaleFactor,0)
    U2_d2=downScale(U2_d2_raw,m_2,n_2,downScaleFactor,0)
    if reverse==True:
        [[U0_d1,U1_d1,U2_d1], [U0_d2,U1_d2,U2_d2]] = [[U0_d2,U1_d2,U2_d2], [U0_d1,U1_d1,U2_d1]]
    # Merge
    merged_aw=mergeTwoPieces(aw_d2,aw_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)
    merged_p_rgh=mergeTwoPieces(p_rgh_d2,p_rgh_d1,m_new,n_new,transfer_2,transfer_1,p_rgh,hor_adjust)
    merged_U0=mergeTwoPieces(U0_d2,U0_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)
    merged_U1=mergeTwoPieces(U1_d2,U1_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)
    merged_U2=mergeTwoPieces(U2_d2,U1_d1,m_new,n_new,transfer_2,transfer_1,0,hor_adjust)

    merged_U=[merged_U0,merged_U1,merged_U2]
    [merged_aw, merged_p_rgh,merged_U]=warmup(merged_aw, merged_p_rgh,merged_U,m_merged,n_merged,p_rgh,warmperiod)

    # Write files
    writeScalar('alpha.water',merged_aw,'0',m_new,n_new)
    writeScalar('p_rgh',merged_p_rgh,str(p_rgh),m_new,n_new)
    writeVector('U',merged_U,'(0 0 0)',m_new,n_new)

def mergeTwoPieces(data_d2,data_d1,m_new,n_new,transfer_2,transfer_1,default_value,hor_adjust):
    block=[default_value]*m_new*n_new
    [centerxi_d1,centerzi_d1,slice_half_width_d1,m_d1]=transfer_1
    [centerxi_d2,centerzi_d2,slice_half_width_d2,m_d2]=transfer_2
    #print(str(centerxi_d1)+' '+str(centerzi_d1)+' '+str(slice_half_width_d1)+' '+str(m_d1))
    #print(str(centerxi_d2)+' '+str(centerzi_d2)+' '+str(slice_half_width_d2)+' '+str(m_d2))
    for xi_d2 in range(slice_half_width_d2*2+1):
        for zi_d2 in range(int(n_new/2)):
            block_xi=int(xi_d2+m_new/2-slice_half_width_d2-hor_adjust/2)
            block_zi=int(zi_d2+n_new/2)
            slice_xi=int(centerxi_d2-slice_half_width_d2+xi_d2)
            slice_zi=int(centerzi_d2-n_new/16+zi_d2)
            if block_xi<m_new and block_xi>0:
                if slice_zi*m_d2+slice_xi<len(data_d2):
                    block[block_xi+block_zi*m_new]=data_d2[slice_zi*m_d2+slice_xi]
      
    for xi_d1 in range(slice_half_width_d1*2+1):
        for zi_d1 in range(-int(n_new/2),0):
            block_xi=int(xi_d1+m_new/2-slice_half_width_d1+hor_adjust/2)
            block_zi=int(zi_d1+n_new/2)
            slice_xi=int(centerxi_d1-slice_half_width_d1+xi_d1)
            slice_zi=int(centerzi_d1-n_new/4+n_new/2-n_new/8+zi_d1)
            if block_xi<m_new and block_xi>0:
                if slice_zi*m_d1+slice_xi<len(data_d1):
                    block[block_xi+block_zi*m_new]=data_d1[slice_zi*m_d1+slice_xi] 
    return block

def shiftOnePiece(data,m,n,transfer,default_value,hor_adjust,ver_adjust):

    block=[default_value]*m*n
    [centerxi,centerzi,slice_half_width,m]=transfer

    for xi in range(slice_half_width*2+1):
        for zi in range(int(n/2)):
            block_xi=int(xi+m/2-slice_half_width-hor_adjust)
            block_zi=int(zi+n/2+ver_adjust)
            slice_xi=int(centerxi-slice_half_width+xi)
            slice_zi=int(centerzi-n/16+zi)
            if block_xi<m and block_xi>0:
                block[block_zi*m+block_xi]=data[slice_zi*m+slice_xi]
    return block

def shiftInit(d2_drop_path,m,n,m_merged,downScaleFactor,threShold,p_rgh,hor_adjust,warmperiod):
    ver_adjust = 0
    n_merged=m_merged*2
    # Downscale alpha.water file
    aw_d2_raw=readFoamScalarOutput(d2_drop_path+'/alpha.water')
    aw_d2=downScale(aw_d2_raw,m,n,downScaleFactor,0)
    transfer_2=pieceInfo(aw_d2,threShold,m*downScaleFactor)
    
    # Downscale p_rgh file
    p_rgh_d2_raw=readFoamScalarOutput(d2_drop_path+'/p_rgh')
    p_rgh_d2=downScale(p_rgh_d2_raw,m,n,downScaleFactor,p_rgh)

    # Downscale U
    [U0_d2_raw,U1_d2_raw,U2_d2_raw]=readFoamVectorOutput(d2_drop_path+'/U')
    U0_d2=downScale(U0_d2_raw,m,n,downScaleFactor,0)
    U1_d2=downScale(U1_d2_raw,m,n,downScaleFactor,0)
    U2_d2=downScale(U2_d2_raw,m,n,downScaleFactor,0)
    
    #shift 
    shifted_aw=shiftOnePiece(aw_d2,m_merged,n_merged,transfer_2,0,hor_adjust,ver_adjust)
    shifted_p_rgh=shiftOnePiece(p_rgh_d2,m_merged,n_merged,transfer_2,p_rgh,hor_adjust,ver_adjust)
    shifted_U0=shiftOnePiece(U0_d2,m_merged,n_merged,transfer_2,0,hor_adjust,ver_adjust)
    shifted_U1=shiftOnePiece(U1_d2,m_merged,n_merged,transfer_2,0,hor_adjust,ver_adjust)
    shifted_U2=shiftOnePiece(U2_d2,m_merged,n_merged,transfer_2,0,hor_adjust,ver_adjust)
    shifted_U=[shifted_U0,shifted_U1,shifted_U2]
    [shifted_aw, shifted_p_rgh,shifted_U]=warmup(shifted_aw, shifted_p_rgh,shifted_U,m_merged,n_merged,p_rgh,warmperiod)
    # Write files
    writeScalar('alpha.water',shifted_aw,'0',m_merged,n_merged)
    writeScalar('p_rgh',shifted_p_rgh,str(p_rgh),m_merged,n_merged)
    writeVector('U',shifted_U,'(0 0 0)',m_merged,n_merged)

def warmup(aw_orig, p_rgh_orig, U_orig, m_new, n_new, p_rgh, warmperiod):
    import os
    writeScalar('alpha.water',aw_orig,'0',m_new,n_new)
    writeScalar('p_rgh',p_rgh_orig,str(p_rgh),m_new,n_new)
    writeVector('U',U_orig,'(0 0 0)',m_new,n_new)    
    writecontrolDict(warmperiod,str(0),str(warmperiod))
    os.system ('interFoam  > log.temp')
    os.system ('rm log.temp')
    aw=readFoamScalarOutput(str(warmperiod)+'/alpha.water')
    p_rgh=readFoamScalarOutput(str(warmperiod)+'/p_rgh')
    U=readFoamVectorOutput(str(warmperiod)+'/U')
    os.system ('rm -rf '+str(warmperiod))
    return [aw, p_rgh, U]


