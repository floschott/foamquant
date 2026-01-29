def CutMovie(series, imrange, readdir, savedir, zcut=None,ycut=None,xcut=None, binvalue=1, verbose=False, EndName='_CutMovie_'):
    import numpy as np
    from tifffile import imsave, imread
    from spam.DIC.deform import binning
    
    from Package.Basic.RangeList import RangeList
    import spam
    
    for imi in imrange:
        # image string index
        imistr = str(imi)
        imistrlen = len(imistr)
        imifordir = (3-imistrlen)*'0'+imistr

        # read image
        if binvalue>1:
            image = spam.DIC.deform.binning(imread(readdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_'+imifordir+'.tif'),binvalue)
        else:
            image = imread(readdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_'+imifordir+'.tif')
            
        # Check case of mid cut: =-1
        if zcut == -1:
            zcut=Z//2
            cutim=image[zcut,:,:]
        elif ycut == -1:
            ycut=Y//2
            cutim=image[:,ycut,:]
        elif xcut == -1:
            xcut=X//2
            cutim=image[:,:,xcut]

        # Check which cut direction
        elif zcut != None:
            cutim=image[zcut,:,:]
        elif ycut != None:
            cutim=image[:,ycut,:]
        elif xcut != None:
            cutim=image[:,:,xcut]
        
        # save movie, image per image
        imsave(series+EndName+ str(imi) + '.tif', cutim, bigtiff=True)
        
        #Verbose
        if verbose:
            print(imi, ': done')
            
def CylMovie(imrange, dirread, dirsave, nameread, namesave, CylRadius, binvalue=1, verbose=False,endread='.tiff',endsave='.tiff'):
    import numpy as np
    from tifffile import imwrite, imread
    from spam.DIC.deform import binning
    
    #from Package.Basic.RangeList import RangeList
    #from Package.Basic.InterpolateCylinder import InterpolateCylinder
    import spam

    # image string index
    imistr = str(imrange[0])
    imistrlen = len(imistr)
    imifordir = (3-imistrlen)*'0'+imistr
    
    # read first image
    if binvalue>1:
        image = spam.DIC.deform.binning(imread(dirread + nameread + imifordir+endread),binvalue)
    else:
        image = imread(dirread + nameread + imifordir+endread)

    Z,Y,X = np.shape(image)
    Npoint = np.int(2*np.pi*CylRadius)

    for imi in imrange:
        if imi > imrange[0]:
            # image string index
            imistr = str(imi)
            imistrlen = len(imistr)
            imifordir = (3-imistrlen)*'0'+imistr
            
            # read image
            if binvalue>1:
                image = spam.DIC.deform.binning(imread(dirread + nameread + imifordir+endread),binvalue)
            else:
                image = imread(dirread + nameread + imifordir+endread)
        
        #Cylindrical interpolation
        interpolated = InterpolateCylinder(image, CylRadius, verbose=True)

        # save movie image per image
        imwrite(dirsave + namesave + imifordir+endsave, interpolated, bigtiff=True)
        
        #Verbose
        if verbose:
            print(imi, ': done')
            
            
def InterpolateCylinder(image, CylRadius, verbose=False, plotfigure=False, nearest=False):
    import numpy as np
    from scipy.interpolate import RectBivariateSpline, NearestNDInterpolator
    import matplotlib.pyplot as plt

    Z,Y,X = image.shape
    # Position of the simulation data
    Npoint = np.int32(2*np.pi*CylRadius)
    angles = np.linspace(0, 2*np.pi, Npoint)
    x = X//2 + np.cos(angles)*CylRadius
    y = Y//2 + np.sin(angles)*CylRadius

    interpolated=np.zeros((Z, Npoint))
    verbi=0
    for zi in range(Z):
        if nearest: # if nearest intepolation
            MeshY,MeshX = np.meshgrid(range(Y),range(X))
            points = np.asarray(list(zip(np.reshape(MeshY, (X*Y)), np.reshape(MeshX, (X*Y)))))
            value = np.reshape(image[zi], (X*Y))
            spline = NearestNDInterpolator(points, value)
            interpolated_list = spline(np.asarray(list(zip(y, x))))
            interpolated[zi] = interpolated_list
        else: # if spline interpolation
            spline = RectBivariateSpline(np.arange(Y), np.arange(X), image[zi])
            interpolated[zi] = spline(y, x, grid=False)
        verbi+=1
        
    if verbose and verbi == 10:
        print(np.round(zi/Z*100), '%')
        verbi=0
    
    if plotfigure:
        fig, ax = plt.subplots(1,2, figsize = (20, 20))
        ax[0].imshow(interpolated, 'bone')
        ax[1].imshow(image[Z//2,:,:], 'bone')
        ax[1].plot(y,x, 'r')
        ax[1].plot(Y//2,X//2, 'xr', markersize=20)
        
        return interpolated, fig
    
    return interpolated

def AssembleMovie(imrange, dirread, dirsave, nameread, namesave,endread='.tiff',endsave='.tiff'):
    import numpy as np
    from tifffile import imwrite, imread
    
    # image string index
    imistr = str(imrange[0])
    imistrlen = len(imistr)
    imifordir = (3-imistrlen)*'0'+imistr
    image = imread(dirread+nameread+ imifordir + endread)
    Y,X = np.shape(image)
    
    Movie=np.zeros((len(imrange),Y,X))

    for i in range(len(imrange)):
        
        # image string index
        imistr = str(imrange[i])
        imistrlen = len(imistr)
        imifordir = (3-imistrlen)*'0'+imistr
        
        Movie[i] = imread(dirread+nameread+ imifordir + endread)
        print('Image',i,': done')

    imwrite(dirsave+namesave+ imifordir + endsave, Movie, bigtiff=True)