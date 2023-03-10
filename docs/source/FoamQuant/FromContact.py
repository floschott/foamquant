def ReadContactTable(nameread, dirread, imrange, verbose=False, endread='.tsv', n0=3, maximumCoordinationNumber=20):
    """
    Read Contact tables and return a dictionary
    
    :param nameread: read image name 
    :type nameread: str
    :param dirread: read image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param verbose: if True, print progression
    :type verbose: Bool
    :param endread: read RegionProperties file extension, default is '.tsv' 
    :type endread: str
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param maximumCoordinationNumber: number of digit for the saving index, default is 3
    :type maximumCoordinationNumber: int
    :return: contact table dictionary {'lab','lab_noedge','Z', 'z','y','x', 'labs', 'conts'}
    """ 
    
    import numpy as np 
    import pandas as pd
    from FoamQuant.Helper import strindex
    import csv
    
    LProperties = []
    
    #Batch loop
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)

        Regprops = pd.read_csv(dirread+nameread+imifordir+endread, sep='\t',engine="python",  quoting=csv.QUOTE_NONE)
        lab = np.asarray(Regprops['lab'])
        lab_noedge = np.asarray(Regprops['lab_noedge'])
        Z = np.asarray(Regprops['Z'])
        z = np.asarray(Regprops['z'])
        y = np.asarray(Regprops['y'])
        x = np.asarray(Regprops['x'])
        
        labi=[];conti=[]
        for i in range(1,maximumCoordinationNumber+1):
            labi.append(np.asarray(Regprops['lab'+str(i)]))
            conti.append(np.asarray(Regprops['cont'+str(i)]))
        labs = np.transpose(labi)
        conts = np.transpose(conti)
        
        if verbose:
                print(nameread+imifordir+': done')
                
        Properties={'lab':lab,'lab_noedge':lab_noedge,'Z':Z, 'z':z,'y':y,'x':x, 'labs':labs, 'conts':conts}
        LProperties.append(Properties)
    
    return LProperties



def Texture(Table, verbose=False):
    """
    From a contact table dictionary, compute the texture
    
    :param Table: contact table dictionary
    :type Table: dict
    :param verbose: if True, print progression
    :type verbose: Bool
    :return: lab withedge, lab withoutedge, centroid, radius, S1,S2,S3, S1z,S1y,S1x, S2z,S2y,S2x, S3z,S3y,S3x, U1,U2,U3, U, type
    """
    
    
    import numpy as np
    
    # lab
    lab=Table['lab']
    lab_noedge=Table['lab_noedge']
    # centroid
    centroid=[]; centroid_noedge = []
    for i in range(len(Table['z'])):
        centroid.append([Table['z'][i],Table['y'][i],Table['x'][i]])
        if Table['lab_noedge'][i]>0:
            centroid_noedge.append(centroid[i])
    # labs
    labs = Table['labs']
    
    # Init properties list
    rad = []
    La=[];Lb=[];Lc=[] 
    Laz=[];Lay=[];Lax=[]
    Lbz=[];Lby=[];Lbx=[]
    Lcz=[];Lcy=[];Lcx=[]
    
    LUa=[];LUb=[];LUc=[]; LU=[]; Ltype=[]
    
    # Texture tensor M
    for h in range(len(lab)):  # line index
        z1,y1,x1 = centroid[h]
        M = np.zeros((3,3))
        count = 0
        
        for k in range(len(labs[h])): # column index
            if labs[h][k] != 0:
                z2,y2,x2 = centroid[labs[h][k]-1]
                [Z,Y,X] = np.asarray([z2,y2,x2]) - np.asarray([z1,y1,x1])
                m = [[Z**2,Z*Y,Z*X],
                     [Y*Z,Y**2,Y*X],
                     [X*Z,X*Y,X**2]]
                M = M + m
                count = count + 1
        
        if count>2:
            M = np.asarray(M)/count   # M = <m> = sum(m)/count
            
            # eig values and vectors
            Val, Vect = np.linalg.eig(M)
            a,b,c = np.sort([Val[0],Val[1],Val[2]])
            oVect = [];eig=[a,b,c]
            for eigi in range(3):
                for val in Val:
                    if val == eig[eigi]:
                        oVect.append(Vect[eigi])
            La.append(a)
            Lb.append(b)
            Lc.append(c)
            Laz.append(oVect[0][0]); Lay.append(oVect[0][1]); Lax.append(oVect[0][2])
            Lbz.append(oVect[1][0]); Lby.append(oVect[1][1]); Lbx.append(oVect[1][2])
            Lcz.append(oVect[2][0]); Lcy.append(oVect[2][1]); Lcx.append(oVect[2][2])
            
            # req for texture (carreful in m^2)
            req = np.power(a*b*c, 1/3)
            rad.append(np.sqrt(req))
            # Strain
            Ua = np.log(a/req)*0.5 # 0.5 for m^2 to m
            Ub = np.log(b/req)*0.5 # 0.5 for m^2 to m
            Uc = np.log(c/req)*0.5 # 0.5 for m^2 to m
            LUa.append(Ua)
            LUb.append(Ub)
            LUc.append(Uc)
            # Strain vM invariant
            LU.append(np.sqrt(0.5*(np.power(Ua-Ub, 2)+np.power(Ua-Uc,2)+np.power(Ub-Uc,2))))
            # Type: oblate or prolate
            if np.abs(Ua) > np.abs(Uc):
                Ltype.append(1)
            else:
                Ltype.append(-1)
            
        else:
            if verbose:
                print('negative texture, not enough contacts:', count)
            La.append(-1)
            Lb.append(-1)
            Lc.append(-1)
            
            rad.append(-1)
                
            Laz.append(-1)
            Lay.append(-1)
            Lax.append(-1)
                
            Lbz.append(-1)
            Lby.append(-1)
            Lbx.append(-1)
                
            Lcz.append(-1)
            Lcy.append(-1)
            Lcx.append(-1)
            
            LUa.append(-1)
            LUb.append(-1)
            LUc.append(-1)
            LU.append(-1)
            Ltype.append(0)
            
    # Keep texture not at the edge
    La_noedge = np.zeros((len(centroid_noedge),1))
    Lb_noedge = np.zeros((len(centroid_noedge),1))
    Lc_noedge = np.zeros((len(centroid_noedge),1))
    
    Lrad_noedge = np.zeros((len(centroid_noedge),1))
                
    Laz_noedge = np.zeros((len(centroid_noedge),1))
    Lay_noedge = np.zeros((len(centroid_noedge),1))
    Lax_noedge = np.zeros((len(centroid_noedge),1))
                
    Lbz_noedge = np.zeros((len(centroid_noedge),1))
    Lby_noedge = np.zeros((len(centroid_noedge),1))
    Lbx_noedge = np.zeros((len(centroid_noedge),1))
                
    Lcz_noedge = np.zeros((len(centroid_noedge),1))
    Lcy_noedge = np.zeros((len(centroid_noedge),1))
    Lcx_noedge = np.zeros((len(centroid_noedge),1))
            
    LUa_noedge = np.zeros((len(centroid_noedge),1))
    LUb_noedge = np.zeros((len(centroid_noedge),1))
    LUc_noedge = np.zeros((len(centroid_noedge),1))
    LU_noedge = np.zeros((len(centroid_noedge),1))
    Ltype_noedge = np.zeros((len(centroid_noedge),1))
    
    labnoedgefromwithedge = np.zeros((len(centroid_noedge),1))
    labnoedgefromwithoutedge = np.zeros((len(centroid_noedge),1))
    
    for h in range(len(La)):
        count=0
        for s in range(len(lab)):
            if lab[s] == lab[h]:
                count+=1
                index=s
        if lab_noedge[index] > 0:  
            La_noedge[lab_noedge[index]-1] = La[h]
            Lb_noedge[lab_noedge[index]-1] = Lb[h]
            Lc_noedge[lab_noedge[index]-1] = Lc[h]
            
            Lrad_noedge[lab_noedge[index]-1] = rad[h]
                
            Laz_noedge[lab_noedge[index]-1] = Laz[h]
            Lay_noedge[lab_noedge[index]-1] = Lay[h]
            Lax_noedge[lab_noedge[index]-1] = Lax[h]
                
            Lbz_noedge[lab_noedge[index]-1] = Lbz[h]
            Lby_noedge[lab_noedge[index]-1] = Lby[h]
            Lbx_noedge[lab_noedge[index]-1] = Lbx[h]
                        
            Lcz_noedge[lab_noedge[index]-1] = Lcz[h]
            Lcy_noedge[lab_noedge[index]-1] = Lcy[h]
            Lcx_noedge[lab_noedge[index]-1] = Lcx[h]
            
            LUa_noedge[lab_noedge[index]-1] = LUa[h]
            LUb_noedge[lab_noedge[index]-1] = LUb[h]
            LUc_noedge[lab_noedge[index]-1] = LUc[h]
            LU_noedge[lab_noedge[index]-1] = LU[h]
            Ltype_noedge[lab_noedge[index]-1] = Ltype[h]
            labnoedgefromwithedge[lab_noedge[index]-1] = lab[h]
            labnoedgefromwithoutedge[lab_noedge[index]-1] = lab_noedge[h]
    
    return labnoedgefromwithedge, labnoedgefromwithoutedge, centroid_noedge, Lrad_noedge, La_noedge,Lb_noedge,Lc_noedge, Laz_noedge,Lay_noedge,Lax_noedge, Lbz_noedge,Lby_noedge,Lbx_noedge, Lcz_noedge,Lcy_noedge,Lcx_noedge, LUa_noedge,LUb_noedge,LUc_noedge, LU_noedge, Ltype_noedge
    
    
def Texture_Batch(nameread, namesave, dirread, dirsave, imrange, verbose=False, endsave='.tsv', n0=3,field=False):
    """
    Run Texture function on a batch of images and save the outputs as .tsv
    
    :param readdir: Labeled images folder
    :type readdir: str
    :param readdir: folder to save the .tsv doc
    :type readdir: str
    :param readend: tiff image saving end, default is '.tiff'
    :type readend: str
    :param imrange: list of image indexes
    :type imrange: int array
    :param IncludeInertia: if True, also return inertial components
    :type IncludeInertia: Bool
    :param verbose: if True, print verbose including the number of labels
    :type verbose: Bool
    """   
    
    import numpy as np 
    from tifffile import imread
    import csv
    from FoamQuant.FromContact import Texture
    from FoamQuant.Helper import strindex
    import os
    
    #Check directory
    isExist = os.path.exists(dirsave)
    print('Path exist:', isExist)
    if not isExist:
        print('Saving path does not exist:\n', isExist)
        return
    
    #Batch loop
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)
        # read contact table
        table = ReadContactTable(nameread, dirread, [imi], verbose=False)
        # texture
        lab, labnoedge, centroid, rad, La,Lb,Lc, Laz,Lay,Lax, Lbz,Lby,Lbx, Lcz,Lcy,Lcx, LUa,LUb,LUc, LU, Ltype = Texture(table[0])
        # Save as TSV
        with open(dirsave + namesave + imifordir + endsave, 'w', newline='') as csvfile:        
            writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['lab','labnoedge','z','y','x', 'rad', 'M1','M2','M3','e1z','e1y','e1x','e2z','e2y','e2x','e3z','e3y','e3x','U1','U2','U3','U','type'])
            for i in range(len(lab)):
                writer.writerow([int(lab[i][0]), 
                                 int(labnoedge[i][0]), 
                                 centroid[i][0],centroid[i][1],centroid[i][2], 
                                 rad[i][0],
                                 La[i][0],Lb[i][0],Lc[i][0],
                                 Laz[i][0],Lay[i][0],Lax[i][0],
                                 Lbz[i][0],Lby[i][0],Lbx[i][0],
                                 Lcz[i][0],Lcy[i][0],Lcx[i][0],
                                 LUa[i][0],LUb[i][0],LUc[i][0],
                                 LU[i][0],Ltype[i][0]])
        if verbose:
            print(namesave+imifordir+': done')
