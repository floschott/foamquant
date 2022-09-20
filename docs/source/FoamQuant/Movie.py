def CutMovie(series, imrange, readdir, savedir, zcut=None,ycut=None,xcut=None, binvalue=1, verbose=False, suffix='_CutMovie_'):
    """
    Save 2D cuts of a series of 3D images for creating a cross-section movie.
    
    :param series: Name of the series
    :type series: str
    :param imrange: image index range
    :type imrange: int numpy array
    :param readdir: read image directory
    :type readdir: str
    :param savedir: save image directory
    :type savedir: str
    :param zcut: Optional, z cut position
    :type zcut: int
    :param ycut: Optional, y cut position
    :type zcut: int
    :param xcut: Optional, x cut position
    :type xcut: int
    :param binvalue: Optional, movie binning
    :type binvalue: int
    :param verbose: Optional, if True, print for loop progress
    :type verbose: Bool
    :param suffix: Saving suffix name
    :type suffix: str
    :return: None
    """
    
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
            
def InterpolateCylinder(image, CylRadius, verbose=False, plotfigure=False):
    """
    Return a 2D cylindrical interpolation of a 3D image
    
    :param image: 3D image
    :type image: numpy array
    :param CylRadius: Cylindrical radius for interpolation
    :type CylRadius: int
    :param verbose: Optional, if True, print the interpolation progress
    :type verbose: Bool
    :param verbose: Optional, if True, print the interpolation progress
    :type verbose: Bool
    :param plotfigure: Optional, if True, plot the interpolation
    :type plotfigure: Bool
    :return: 2D image numoy array
    """
    
    import numpy as np
    from scipy.interpolate import RectBivariateSpline

    Z,Y,X = image.shape
    # Position of the simulation data
    Npoint = np.int(2*np.pi*CylRadius)
    angles = np.linspace(0, 2*np.pi, Npoint)
    x = X//2 + np.cos(angles)*CylRadius
    y = Y//2 + np.sin(angles)*CylRadius

    interpolated=np.zeros((Z, Npoint))
    verbi=0
    for zi in range(Z):
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

def CylMovie(series, imrange, readdir, savedir, CylRadius, binvalue=1, verbose=False, EndName='_MovieCylinder_'):
    """
    Save series of 2D cylindrical interpolated images of a series of 3D images for creating a cross-section movie.
    
    :param series: Name of the series
    :type series: str
    :param imrange: image index range
    :type imrange: int numpy array
    :param readdir: read image directory
    :type readdir: str
    :param savedir: save image directory
    :type savedir: str
    :param CylRadius: Cylindrical radius for interpolation
    :type CylRadius: int
    :param binvalue: Optional, movie binning
    :type binvalue: int
    :param verbose: Optional, if True, print for loop progress
    :type verbose: Bool
    :param suffix: Saving suffix name
    :type suffix: str
    :return: None
    """
    
    import numpy as np
    from tifffile import imsave, imread
    from spam.DIC.deform import binning
    
    from Package.Basic.RangeList import RangeList
    from Package.Basic.InterpolateCylinder import InterpolateCylinder
    import spam

    # image string index
    imistr = str(imrange[0])
    imistrlen = len(imistr)
    imifordir = (3-imistrlen)*'0'+imistr
    
    # read first image
    if binvalue>1:
        image = spam.DIC.deform.binning(imread(readdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_'+imifordir+'.tif'),binvalue)
    else:
        image = imread(readdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_'+imifordir+'.tif')

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
                image = spam.DIC.deform.binning(imread(readdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_'+imifordir+'.tif'),binvalue)
            else:
                image = imread(readdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_'+imifordir+'.tif')
        
        #Cylindrical interpolation
        interpolated = InterpolateCylinder(image, CylRadius, verbose=True)

        # save movie image per image
        imsave(series+EndName+ str(imi) + '.tif', interpolated, bigtiff=True)
        
        #Verbose
        if verbose:
            print(imi, ': done')
            
            
def AssembleMovie(series, imrange, readdir, savedir, suffixread='_CutMovie_',suffixsave='_FullCutMovie'):
    """
    Assemble series of 2D images ans save one 3D numpy image (movie).
    
    :param series: Name of the series
    :type series: str
    :param imrange: image index range
    :type imrange: int numpy array
    :param readdir: read image directory
    :type readdir: str
    :param savedir: save image directory
    :type savedir: str
    :param suffixread: Read suffix name
    :type suffixread: str
    :param suffixsave: Saving suffix name
    :type suffixsave: str
    :return: None
    """
    
    import numpy as np
    from tifffile import imsave, imread
    
    image = imread(series+EndNameread+ str(imrange[i]) + '.tif')
    Y,X = np.shape(image)
    
    Movie=np.zeros((len(imrange),Y,X)

    for i in range(len(imrange)):
        Movie[i] = imread(series+EndNameread+ str(imrange[i]) + '.tif')
        print('Image',i,': done')

    imsave(series+EndNamesave+'.tif', Movie, bigtiff=True)            
         
