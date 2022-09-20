def ReadRaw(series, imi, rawdir, zN=800, top=0, bottom=None):
    """
    Read raw images saved as 2D tiffs.
    
    :param series: Name of the series
    :type series: str
    :param imi: image index
    :type imi: int
    :param rawdir: raw image directory
    :type rawdir: str
    :param zN: number of 2D images
    :type zN: int
    :param top: starting z index
    :type top: int
    :param bottom: ending z index
    :type bottom: int
    :return: 3D numpy array
    """
    
    import numpy as np
    from tifffile import imread
    
    # image string index
    imistr = str(imi)
    imistrlen = len(imistr)
    imifordir = (5-imistrlen)*'0'+imistr
    
    
    if bottom != None:
        zN=bottom
    
    # Init image
    image = np.zeros((zN-top,2016,2016))
    
    for zi in range(top, zN):
        # horyzontal slice string index
        zistr = str(zi+1)
        zistrlen = len(zistr)
        zifordir = (3-zistrlen)*'0'+zistr
        
        # horyzontal slice directory
        imdir = rawdir + '/' + series + '/' + 'rec_8bit_phase_'+imifordir + '/' + series+'_'+imifordir+'_'+zifordir+'.rec.8bit.tif'
        
        # read slice and put into 3D image
        image[zi-top] = imread(imdir)
        
    # return 3D image
    return image




def RangeList(i1, i2, verbose=False):
    """
    Create a int numpy array between i1 and i2.
    
    :param i1: begin index
    :type i1: str
    :param i2: end index
    :type i2: int
    :param verbose: if True, print the array
    :type verbose: Bool
    :return: numpy int array
    """
    
    import numpy as np
    List = np.uint8(np.linspace(i1,i2,i2-i1+1))
    if verbose:
        print(List)
    return List
