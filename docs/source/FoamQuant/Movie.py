def CutMovie(series, imrange, readdir, savedir,
             zcut=None, ycut=None, xcut=None,
             binvalue=1, verbose=False, EndName='_CutMovie_'):
    """
    CutMovie
    --------
    Extracts 2D slices from a 3D image sequence and saves them as individual images
    (movie frames).

    Parameters
    ----------
    series : str
        Name of the image series.
    imrange : iterable of int
        Image indices to process.
    readdir : str
        Directory containing the input images.
    savedir : str
        Directory where the cut images will be saved.
    zcut : int or None, optional
        Z index for slicing. If -1, the mid-plane is selected.
    ycut : int or None, optional
        Y index for slicing. If -1, the mid-plane is selected.
    xcut : int or None, optional
        X index for slicing. If -1, the mid-plane is selected.
    binvalue : int, optional
        Binning factor applied before slicing.
    verbose : bool, optional
        If True, prints progress information.
    EndName : str, optional
        Suffix added to the output filenames.

    Outputs
    -------
    None
        Saves one TIFF image per input image.
    """

    import numpy as np
    from tifffile import imsave, imread
    from spam.DIC.deform import binning
    import spam

    for imi in imrange:
        imistr = str(imi)
        imifordir = (3 - len(imistr)) * '0' + imistr

        if binvalue > 1:
            image = spam.DIC.deform.binning(
                imread(readdir + '/2_RemoveBackground/' + series + '/' +
                       series + '_RemoveBackground_' + imifordir + '.tif'),
                binvalue
            )
        else:
            image = imread(readdir + '/2_RemoveBackground/' + series + '/' +
                           series + '_RemoveBackground_' + imifordir + '.tif')

        Z, Y, X = image.shape

        if zcut == -1:
            cutim = image[Z // 2, :, :]
        elif ycut == -1:
            cutim = image[:, Y // 2, :]
        elif xcut == -1:
            cutim = image[:, :, X // 2]
        elif zcut is not None:
            cutim = image[zcut, :, :]
        elif ycut is not None:
            cutim = image[:, ycut, :]
        elif xcut is not None:
            cutim = image[:, :, xcut]

        imsave(savedir + series + EndName + str(imi) + '.tif',
               cutim, bigtiff=True)

        if verbose:
            print(imi, ': done')

def CylMovie(imrange, dirread, dirsave,
             nameread, namesave,
             CylRadius, binvalue=1,
             verbose=False, endread='.tiff', endsave='.tiff'):
    """
    CylMovie
    --------
    Generates a cylindrical projection movie from a sequence of 3D images.

    Parameters
    ----------
    imrange : iterable of int
        Image indices to process.
    dirread : str
        Directory containing the input images.
    dirsave : str
        Directory where the output images will be saved.
    nameread : str
        Base name of the input images.
    namesave : str
        Base name of the output images.
    CylRadius : float
        Radius of the cylinder used for interpolation.
    binvalue : int, optional
        Binning factor applied before interpolation.
    verbose : bool, optional
        If True, prints progress information.
    endread : str, optional
        File extension of input images.
    endsave : str, optional
        File extension of output images.

    Outputs
    -------
    None
        Saves one interpolated cylindrical image per input image.
    """

    import numpy as np
    from tifffile import imwrite, imread
    from spam.DIC.deform import binning
    import spam

    for imi in imrange:
        imistr = str(imi)
        imifordir = (3 - len(imistr)) * '0' + imistr

        if binvalue > 1:
            image = spam.DIC.deform.binning(
                imread(dirread + nameread + imifordir + endread),
                binvalue
            )
        else:
            image = imread(dirread + nameread + imifordir + endread)

        interpolated = InterpolateCylinder(image, CylRadius, verbose=True)

        imwrite(dirsave + namesave + imifordir + endsave,
                interpolated, bigtiff=True)

        if verbose:
            print(imi, ': done')

def InterpolateCylinder(image, CylRadius,
                        verbose=False, plotfigure=False, nearest=False):
    """
    InterpolateCylinder
    -------------------
    Interpolates a 3D image onto a cylindrical surface.

    Parameters
    ----------
    image : ndarray (Z, Y, X)
        Input 3D image.
    CylRadius : float
        Radius of the cylinder.
    verbose : bool, optional
        If True, prints progress information.
    plotfigure : bool, optional
        If True, returns a matplotlib figure showing the interpolation.
    nearest : bool, optional
        If True, uses nearest-neighbor interpolation instead of spline.

    Outputs
    -------
    interpolated : ndarray (Z, N)
        Cylindrical interpolated image.
    fig : matplotlib.figure.Figure, optional
        Returned only if plotfigure=True.
    """

    import numpy as np
    from scipy.interpolate import RectBivariateSpline, NearestNDInterpolator
    import matplotlib.pyplot as plt

    Z, Y, X = image.shape
    Npoint = np.int32(2 * np.pi * CylRadius)

    angles = np.linspace(0, 2 * np.pi, Npoint)
    x = X // 2 + np.cos(angles) * CylRadius
    y = Y // 2 + np.sin(angles) * CylRadius

    interpolated = np.zeros((Z, Npoint))

    for zi in range(Z):
        if nearest:
            MeshY, MeshX = np.meshgrid(range(Y), range(X))
            points = np.column_stack((MeshY.ravel(), MeshX.ravel()))
            values = image[zi].ravel()
            interp = NearestNDInterpolator(points, values)
            interpolated[zi] = interp(np.column_stack((y, x)))
        else:
            spline = RectBivariateSpline(np.arange(Y),
                                         np.arange(X),
                                         image[zi])
            interpolated[zi] = spline(y, x, grid=False)

        if verbose and zi % 10 == 0:
            print(np.round(zi / Z * 100), '%')

    if plotfigure:
        fig, ax = plt.subplots(1, 2, figsize=(20, 20))
        ax[0].imshow(interpolated, cmap='bone')
        ax[1].imshow(image[Z // 2], cmap='bone')
        ax[1].plot(y, x, 'r')
        ax[1].plot(Y // 2, X // 2, 'xr', markersize=20)
        return interpolated, fig

    return interpolated

def AssembleMovie(imrange, dirread, dirsave,
                  nameread, namesave,
                  endread='.tiff', endsave='.tiff'):
    """
    AssembleMovie
    -------------
    Assembles a sequence of 2D images into a 3D movie stack.

    Parameters
    ----------
    imrange : iterable of int
        Image indices to assemble.
    dirread : str
        Directory containing the input images.
    dirsave : str
        Directory where the movie will be saved.
    nameread : str
        Base name of the input images.
    namesave : str
        Base name of the output movie.
    endread : str, optional
        File extension of input images.
    endsave : str, optional
        File extension of output movie.

    Outputs
    -------
    None
        Saves a 3D TIFF stack (movie).
    """

    import numpy as np
    from tifffile import imwrite, imread

    imistr = str(imrange[0])
    imifordir = (3 - len(imistr)) * '0' + imistr
    image = imread(dirread + nameread + imifordir + endread)

    Y, X = image.shape
    Movie = np.zeros((len(imrange), Y, X))

    for i, imi in enumerate(imrange):
        imistr = str(imi)
        imifordir = (3 - len(imistr)) * '0' + imistr
        Movie[i] = imread(dirread + nameread + imifordir + endread)
        print('Image', i, ': done')

    imwrite(dirsave + namesave + imifordir + endsave,
            Movie, bigtiff=True)
