import os

def writeblockMeshDict(domainwidth,m_grid,n_grid):
    os.system('cp system/blockMeshDict system/blockMeshDict.txt')    
    with open('system/blockMeshDict.txt', 'r') as file:
        data = file.readlines()
    data[16]='convertToMeters '+str(domainwidth)+';\n'
    data[35]='    hex (0 1 2 3 4 5 6 7) ('+str(m_grid)+' '+str(n_grid)+' 1) simpleGrading (1 1 1)\n'
    with open('system/blockMeshDict.txt', 'w') as file:
        file.writelines(data)
    os.system('cp system/blockMeshDict.txt system/blockMeshDict')
    os.system('rm system/blockMeshDict.txt')

def writesetFieldsDict(domainwidth,dropdiameter):
    dropradius=dropdiameter/2
    dropcenterx=domainwidth/2
    dropcentery=-1.5*dropdiameter
    os.system('cp system/setFieldsDict system/setFieldsDict.txt')    
    with open('system/setFieldsDict.txt', 'r') as file:
        data = file.readlines()
    data[26]='        p1 ('+str(dropcenterx)+' '+str(dropcentery)+' -1); p2 ('+str(dropcenterx)+' '+str(dropcentery)+' 1); radius '+str(dropradius)+';\n'
    with open('system/setFieldsDict.txt', 'w') as file:
        file.writelines(data)
    os.system('cp system/setFieldsDict.txt system/setFieldsDict')
    os.system('rm system/setFieldsDict.txt')

def writecontrolDict(writeInterval,t0,t1):
    with open('system/controlDict.txt', 'r') as file:
        data = file.readlines()
    data[31]='writeInterval       '+str(writeInterval)+';\n' 
    data[21]='startTime       '+str(t0)+';\n' 
    data[25]='endTime         '+str(t1)+';\n'
    with open('system/controlDict.txt', 'w') as file:
        file.writelines(data)
    os.system('cp system/controlDict.txt system/controlDict')

def slidevertical(writeInterval,keyword,pressure):

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

    os.system('cp processor0/'+str(writeInterval)+'/alpha.water processor0/'+keyword)
    #os.system('cp processor0/'+str(writeInterval)+'/alphaPhi0.water processor0/'+keyword)
    os.system('cp processor0/'+str(writeInterval)+'/p processor0/'+keyword)
    os.system('cp processor0/'+str(writeInterval)+'/p_rgh processor0/'+keyword)
    #os.system('cp processor0/'+str(writeInterval)+'/phi processor0/'+keyword)
    os.system('cp processor0/'+str(writeInterval)+'/U processor0/'+keyword)

    os.system('cp processor1/'+str(writeInterval)+'/alpha.water processor1/'+keyword)
    #os.system('cp processor1/'+str(writeInterval)+'/alphaPhi0.water processor1/'+keyword)
    os.system('cp processor1/'+str(writeInterval)+'/p processor1/'+keyword)
    os.system('cp processor1/'+str(writeInterval)+'/p_rgh processor1/'+keyword)
    #os.system('cp processor1/'+str(writeInterval)+'/phi processor1/'+keyword)
    os.system('cp processor1/'+str(writeInterval)+'/U processor1/'+keyword)
    

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
    data=data.replace('    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> \n','    bot\n    {\n        type            calculated;\n        value           nonuniform List<scalar> 0();\n    }\n    procBoundary2to1\n    {\n        type            processor;\n        value           nonuniform List<scalar> \n')    
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

def dropinfo(keyword,m,pn,threshold):
    tempaw=[]
    tempnx=[]
    tempny=[]
    tempUx=[]
    tempUz=[]
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
                    tempUx.append(float((Udata[awi+22])[1:-2].split(' ')[0]))
                    tempUz.append(float((Udata[awi+22])[1:-2].split(' ')[1]))

    centernx=sum(tempnx)/len(tempnx)
    centerny=sum(tempny)/len(tempny)
    averageaw=sum(tempaw)/len(tempaw)
    averageUx=sum(tempUx)/len(tempUx)
    averageUz=sum(tempUz)/len(tempUz)
    datatowrite='nx ny alpha.water Ux(m/s) Uz(m/s)\n'
    for li in range(len(tempnx)):
        datatowrite+=str(tempnx[li])+' '+str(tempny[li])+' '+str(tempaw[li])+' '+str(tempUx[li])+' '+str(tempUz[li])+'\n'
    os.system('touch droploc/drop'+keyword)

    with open('droploc/drop'+keyword,'w') as file:
        file.write(datatowrite)

    return [centernx,centerny,averageaw,averageUx,averageUz]    
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



  

 