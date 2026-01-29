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
            
            Laz.append(M[0][0]); Lay.append(M[0][1]); Lax.append(M[0][2])
            Lbz.append(M[1][0]); Lby.append(M[1][1]); Lbx.append(M[1][2])
            Lcz.append(M[2][0]); Lcy.append(M[2][1]); Lcx.append(M[2][2])
            
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
            
def Read_Texture(nameread, dirread, imrange, verbose=False, endread='.tsv', n0=3):            
    import numpy as np 
    import pandas as pd
    from FoamQuant.Helper import strindex
    import csv
    
    lab=[]; labnoedge=[]; z=[];y=[];x=[]; rad=[]
    La=[];Lb=[];Lc=[] 
    Laz=[];Lay=[];Lax=[]
    Lbz=[];Lby=[];Lbx=[]
    Lcz=[];Lcy=[];Lcx=[]
    LUa=[];LUb=[];LUc=[]; LU=[]; Ltype=[]
    
    #Batch loop
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)

        Texture = pd.read_csv(dirread+nameread+imifordir+endread, sep='\t',engine="python",  quoting=csv.QUOTE_NONE)
        lab.append(np.asarray(Texture['lab']))
        labnoedge.append(np.asarray(Texture['labnoedge']))
        z.append(np.asarray(Texture['z']))
        y.append(np.asarray(Texture['y']))
        x.append(np.asarray(Texture['x']))
        rad.append(np.asarray(Texture['rad']))
        La.append(np.asarray(Texture['M1']))
        Lb.append(np.asarray(Texture['M2']))
        Lc.append(np.asarray(Texture['M3'])) 
            
        Laz.append(np.asarray(Texture['e1z']))
        Lay.append(np.asarray(Texture['e1y']))
        Lax.append(np.asarray(Texture['e1x']))
        Lbz.append(np.asarray(Texture['e2z']))
        Lby.append(np.asarray(Texture['e2y']))
        Lbx.append(np.asarray(Texture['e2x']))
        Lcz.append(np.asarray(Texture['e3z']))
        Lcy.append(np.asarray(Texture['e3y']))
        Lcx.append(np.asarray(Texture['e3x']))
        
        LUa.append(np.asarray(Texture['U1']))
        LUb.append(np.asarray(Texture['U2']))
        LUc.append(np.asarray(Texture['U3']))
        LU.append(np.asarray(Texture['U']))
        Ltype.append(np.asarray(Texture['type']))
        
        if verbose:
            print(nameread+imifordir+': done')
        
    lab=np.concatenate(lab)
    labnoedge=np.concatenate(labnoedge)
    z=np.concatenate(z)
    y=np.concatenate(y)
    x=np.concatenate(x)
    rad=np.concatenate(rad)
    La=np.concatenate(La)
    Lb=np.concatenate(Lb)
    Lc=np.concatenate(Lc)
    Laz=np.concatenate(Laz)
    Lay=np.concatenate(Lay)
    Lax=np.concatenate(Lax)
    Lbz=np.concatenate(Lbz)
    Lby=np.concatenate(Lby)
    Lbx=np.concatenate(Lbx)
    Lcz=np.concatenate(Lcz)
    Lcy=np.concatenate(Lcy)
    Lcx=np.concatenate(Lcx)
    LUa=np.concatenate(LUa)
    LUb=np.concatenate(LUb)
    LUc=np.concatenate(LUc)
    LU=np.concatenate(LU)
    Ltype=np.concatenate(Ltype)
    
    Properties={'lab':lab,'labnoedge':labnoedge,'z':z,'y':y,'x':x,'rad':rad,'S1':La,'S2':Lb,'S3':Lc,
                'e1z':Laz,'e1y':Lay,'e1x':Lax,'e2z':Lbz,'e2y':Lby,'e2x':Lbx,'e3z':Lcz,'e3y':Lcy,'e3x':Lcx,
                'U1':LUa,'U2':LUb,'U3':LUc,'U':LU,'type':Ltype}
    return Properties

def ContactProp(image, field=False, verbose=False):
    """
    Return basic region properties from a labeled contact image: labels [0], centroids [1], volumes [2], (inertia components [3])
    
    :param image: 3D image 
    :type image: int numpy array
    :param IncludeInertia: if True, also return inertial components
    :type IncludeInertia: Bool
    :return: array of labels [0], centroids [1], volumes [2], (inertia components [3])
    """
    import numpy as np
    from skimage.measure import regionprops
    
    if not field:
        field = [0,np.inf,0,np.inf,0,np.inf]
    Reg = regionprops(image)
    
    lab = []; centroid = []; vol = []; areafit = []; radfit = []; ecc = []

    La=[];Lb=[];Lc=[] 
    Laz=[];Lay=[];Lax=[]
    Lbz=[];Lby=[];Lbx=[]
    Lcz=[];Lcy=[];Lcx=[]
    
    for reg in Reg:
        z,y,x = reg.centroid
        if z>=field[0] and z<=field[1] and y>=field[2] and y<=field[3] and x>=field[4] and x<=field[5]:
            lab.append(reg.label)
            centroid.append([z,y,x])
            V = reg.area
            vol.append(V)

            eig1,eig2,eig3 = reg.inertia_tensor_eigvals
            I = reg.inertia_tensor
            Val, Vect = np.linalg.eig(I)
            
            oVect = [];eig=[eig1,eig2,eig3]
            for eigi in range(3):
                for val in Val:
                    if round(val) == round(eig[eigi]):
                        oVect.append(Vect[eigi])
            
            if eig2+eig3-eig1>0 and eig1+eig3-eig2>0 and eig1+eig2-eig3>0 and len(oVect)==3:
                # Semi-axes
                a = np.sqrt(5/2*(eig2+eig3-eig1))
                b = np.sqrt(5/2*(eig1+eig3-eig2))
                c = np.sqrt(5/2*(eig1+eig2-eig3))
                # Ordered semi-axes
                sa,sb,sc = np.sort([a,b,c])
                La.append(sa)
                Lb.append(sb)
                Lc.append(sc)
                # Ordered corresponding eig-vectors
                ooVect=[]; eig=[a,b,c]
                for si in [sa,sb,sc]:
                    for eigi in range(3):
                        if si == eig[eigi]:
                            ooVect.append(oVect[eigi])
                oVect=ooVect
                a=sa; b=sb; c=sc
                Laz.append(oVect[0][0]); Lay.append(oVect[0][1]); Lax.append(oVect[0][2])
                Lbz.append(oVect[1][0]); Lby.append(oVect[1][1]); Lbx.append(oVect[1][2])
                Lcz.append(oVect[2][0]); Lcy.append(oVect[2][1]); Lcx.append(oVect[2][2])
                
                # Area
                areafit.append(2*np.pi*c*b)
                #radius
                radfit.append(np.mean([b,c]))
                #eccentricity
                ecc.append(np.sqrt(np.power(c,2)-np.power(b,2)))
                
            else:
                La.append(-1)
                Lb.append(-1)
                Lc.append(-1)
                
                Laz.append(-1)
                Lay.append(-1)
                Lax.append(-1)
                
                Lbz.append(-1)
                Lby.append(-1)
                Lbx.append(-1)
                
                Lcz.append(-1)
                Lcy.append(-1)
                Lcx.append(-1)
                
                areafit.append(-1)     
                radfit.append(-1)
                ecc.append(-1)
                
                if verbose:
                    print('Error: toos small contact, cannot define properties of 0 thickness surface: ', eig2+eig3-eig1, eig1+eig3-eig2, eig1+eig2-eig3)
                #print('Error: negative shape value, coord',centroid[-1], 'vol', vol)
        
    return lab,centroid,vol, areafit,radfit,ecc, La,Lb,Lc, Laz,Lay,Lax, Lbz,Lby,Lbx, Lcz,Lcy,Lcx



def ContactProp_Batch(nameread, namesave, dirread, dirsave, imrange, closing=False, verbose=False, endread='.tiff', endsave='.tsv', n0=3,field=False):
    """
    Run RegionProp function on a batch of images and save the outputs as .tsv
    
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
    from skimage.morphology import closing
    from skimage.morphology import ball 
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
        # read image
        image = imread(dirread + nameread + imifordir + endread)
        if closing:
            image = closing(image, ball(1))
        
        lab,centroid,vol, areafit,radfit,ecc, La,Lb,Lc, Laz,Lay,Lax, Lbz,Lby,Lbx, Lcz,Lcy,Lcx = ContactProp(image, field=field, verbose=verbose)
        # Save as TSV
        with open(dirsave + namesave + imifordir + endsave, 'w', newline='') as csvfile:        
            writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['lab','z','y','x', 'vol','areafit','radfit',
                             'ecc','a','B','C','e1z','e1y','e1x','e2z','e2y','e2x','e3z','e3y','e3x'])
            for i in range(len(lab)):
                writer.writerow([lab[i], centroid[i][0],centroid[i][1],centroid[i][2], 
                                 vol[i],areafit[i], radfit[i],ecc[i],
                                 La[i],Lb[i],Lc[i],
                                 Laz[i],Lay[i],Lax[i],
                                 Lbz[i],Lby[i],Lbx[i],
                                 Lcz[i],Lcy[i],Lcx[i]])
        if verbose:
            print(namesave+imifordir+': done')

            
def Read_ContactProp(nameread, dirread, imrange, verbose=False, endread='.tsv', n0=3):            
    import numpy as np 
    import pandas as pd
    from FoamQuant.Helper import strindex
    import csv
    
    lab=[]; z=[];y=[];x=[]; vol=[]; areafit=[]
    radfit=[]; ecc=[]
    
    La=[];Lb=[];Lc=[] 
    Laz=[];Lay=[];Lax=[]
    Lbz=[];Lby=[];Lbx=[]
    Lcz=[];Lcy=[];Lcx=[]
    
    #Batch loop
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)

        Contactprops = pd.read_csv(dirread+nameread+imifordir+endread, sep='\t',engine="python",  quoting=csv.QUOTE_NONE)
        lab.append(np.asarray(Contactprops['lab']))
        z.append(np.asarray(Contactprops['z']))
        y.append(np.asarray(Contactprops['y']))
        x.append(np.asarray(Contactprops['x']))
        vol.append(np.asarray(Contactprops['vol']))
        areafit.append(np.asarray(Contactprops['areafit']))
        radfit.append(np.asarray(Contactprops['radfit']))
        ecc.append(np.asarray(Contactprops['ecc']))
        La.append(np.asarray(Contactprops['a']))
        Lb.append(np.asarray(Contactprops['B']))
        Lc.append(np.asarray(Contactprops['C']))
            
        Laz.append(np.asarray(Contactprops['e1z']))
        Lay.append(np.asarray(Contactprops['e1y']))
        Lax.append(np.asarray(Contactprops['e1x']))
        Lbz.append(np.asarray(Contactprops['e2z']))
        Lby.append(np.asarray(Contactprops['e2y']))
        Lbx.append(np.asarray(Contactprops['e2x']))
        Lcz.append(np.asarray(Contactprops['e3z']))
        Lcy.append(np.asarray(Contactprops['e3y']))
        Lcx.append(np.asarray(Contactprops['e3x']))
        
        if verbose:
                print(nameread+imifordir+': done')
        
    lab=np.concatenate(lab)
    z=np.concatenate(z)
    y=np.concatenate(y)
    x=np.concatenate(x)
    vol=np.concatenate(vol)
    areafit=np.concatenate(areafit)
    radfit=np.concatenate(radfit)
    ecc=np.concatenate(ecc)
    La=np.concatenate(La)
    Lb=np.concatenate(Lb)
    Lc=np.concatenate(Lc)
    Laz=np.concatenate(Laz)
    Lay=np.concatenate(Lay)
    Lax=np.concatenate(Lax)
    Lbz=np.concatenate(Lbz)
    Lby=np.concatenate(Lby)
    Lbx=np.concatenate(Lbx)
    Lcz=np.concatenate(Lcz)
    Lcy=np.concatenate(Lcy)
    Lcx=np.concatenate(Lcx)
    
    Properties={'lab':lab,'z':z,'y':y,'x':x,'vol':vol,'areafit':areafit,'radfit':radfit,'ecc':ecc,'a':La,'B':Lb,'C':Lc,
                'e1z':Laz,'e1y':Lay,'e1x':Lax,'e2z':Lbz,'e2y':Lby,'e2x':Lbx,'e3z':Lcz,'e3y':Lcy,'e3x':Lcx}
    return Properties

def Read_ContactPair(nameread, dirread, imrange, verbose=False, endread='.tsv', n0=3, extended=False):            
    import numpy as np 
    import pandas as pd
    from FoamQuant.Helper import strindex
    import csv
    
    LPairs=[]
    
    #Batch loop
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)

        Contactprops = pd.read_csv(dirread+nameread+imifordir+endread, sep='\t',engine="python",  quoting=csv.QUOTE_NONE)
        cont = np.asarray(Contactprops['cont'])
        lab1 = np.asarray(Contactprops['lab1'])
        lab2 = np.asarray(Contactprops['lab2'])
        
        if extended:
            z1 = np.asarray(Contactprops['z1'])
            z2 = np.asarray(Contactprops['z2'])
            y1 = np.asarray(Contactprops['y1'])
            y2 = np.asarray(Contactprops['y2'])
            x1 = np.asarray(Contactprops['x1'])
            x2 = np.asarray(Contactprops['x2'])
            
            Pairs={'cont':cont,'lab1':lab1,'lab2':lab2,'z1':z1,'z2':z2,'y1':y1,'y2':y2,'x1':x1,'x2':x2}
            LPairs.append(Pairs)
        else:
            Pairs={'cont':cont,'lab1':lab1,'lab2':lab2}
            LPairs.append(Pairs)
        if verbose:
                print(nameread+imifordir+': done')
    return LPairs


def Translate_Pairs(Pairs2, tracking12):
    import numpy as np
    from tqdm import tqdm
    
    # List of pairs at time 2
    cont_t2 = Pairs2['cont']
    lab1_t2 = Pairs2['lab1'] # label n1 at time t2
    lab2_t2 = Pairs2['lab2'] # label n2 at time t2

    # Label tracking from time 1 to time 2
    labs_t1 = tracking12['lab1'] #labs at time t1
    labs_t2 = tracking12['lab2'] #labs at time t2
    
    x1_t1 = tracking12['x1']
    x2_t1 = tracking12['x2']
    y1_t1 = tracking12['y1']
    y2_t1 = tracking12['y2']
    z1_t1 = tracking12['z1']
    z2_t1 = tracking12['z2']

    cont_tsl = cont_t2
    lab1_tsl = []
    lab2_tsl = []
    x1_tsl = []
    x2_tsl = []
    y1_tsl = []
    y2_tsl = []
    z1_tsl = []
    z2_tsl = []

    for pi in tqdm(range(len(lab1_t2))):
        found1 = False
        found2 = False
        for tbli in range(len(labs_t2)):
            if not found1 and lab1_t2[pi] == labs_t2[tbli]: # if contact lab at time t2 == tracking lab at time t2
                lab1_tsl.append(labs_t1[tbli])
                
                x1_tsl.append(x1_t1[tbli])
                y1_tsl.append(y1_t1[tbli])
                z1_tsl.append(z1_t1[tbli])
                found1 = True # label n1 was traked
                
            if not found2 and lab2_t2[pi] == labs_t2[tbli]: # if contact lab at time t2 == tracking lab at time t2
                lab2_tsl.append(labs_t1[tbli])
                
                x2_tsl.append(x1_t1[tbli])
                y2_tsl.append(y1_t1[tbli])
                z2_tsl.append(z1_t1[tbli])
                found2 = True # label n2 was traked
                
        if not found1:
            lab1_tsl.append(-1)#np.nan)
            x1_tsl.append(-1)#np.nan)
            y1_tsl.append(-1)#np.nan)
            z1_tsl.append(-1)#np.nan)
        if not found2:
            lab2_tsl.append(-1)#np.nan)
            x2_tsl.append(-1)#np.nan)
            y2_tsl.append(-1)#np.nan)
            z2_tsl.append(-1)#np.nan)

    cont_tsl=np.asarray(cont_tsl)
    lab1_tsl=np.asarray(lab1_tsl)
    lab2_tsl=np.asarray(lab2_tsl)

    pairs_tsl = {'cont':cont_tsl,'lab1':lab1_tsl,'lab2':lab2_tsl, 'x1':x1_tsl, 'x2':x2_tsl,'y1':y1_tsl, 'y2':y2_tsl, 'z1':z1_tsl, 'z2':z2_tsl}
    #Pairs_tsl.append(pairs_tsl)
        
    return pairs_tsl

def Translate_Pairs_Batch(nameread, namesave, nameread_track, dirread_track, dirread, dirsave, imrange, endsave='.tsv', n0=3):
    import csv
    import numpy as np
    from FoamQuant.Helper import strindex
    from FoamQuant.Tracking import Read_LabelTracking
    from FoamQuant.FromContact import Read_ContactPair
    
    for i in range(len(imrange)-1):
        ind1 = imrange[i]
        ind2 = imrange[i+1]
        imifordir2 = strindex(ind2, n0=n0)

        tracking12 = Read_LabelTracking(nameread_track, dirread_track, [ind1, ind2], verbose=True)
        Pairs2 = Read_ContactPair(nameread, dirread, [ind2], verbose=True)

        Pairs2_tsl = Translate_Pairs(Pairs2[0], tracking12[0])

        with open(dirsave + namesave + imifordir2 + endsave, 'w', newline='') as csvfile:        
                writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

                #first line
                writer.writerow(['cont','lab1','lab2','z1','z2','y1','y2','x1','x2'])
                #data line by line
                for i in range(len(Pairs2_tsl['lab1'])):
                    writer.writerow([Pairs2_tsl['cont'][i], 
                                     Pairs2_tsl['lab1'][i],Pairs2_tsl['lab2'][i], 
                                     Pairs2_tsl['z1'][i],Pairs2_tsl['z2'][i],
                                     Pairs2_tsl['y1'][i],Pairs2_tsl['y2'][i],
                                     Pairs2_tsl['x1'][i],Pairs2_tsl['x2'][i]])