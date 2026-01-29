def strindex(i, n0):
    """
    Return str index written on n0 digit
    
    :param i: index 
    :type i: int
    :param n0: number of 0 digit
    :type n0: int
    :return: str index
    """    
    
    istr = str(i)
    istrlen = len(istr)
    fullistr = (n0-istrlen)*'0'+istr
    return fullistr

def RangeList(i1, i2, verbose=False):
    """
    Return an index range
    
    :param i1: first index 
    :type i1: int
    :param i2: last index 
    :type i2: int
    :param verbose: If True, print the range
    :type verbose: Bool
    :return: int numpy array
    """
    
    import numpy as np
    List = np.uint8(np.linspace(i1,i2,i2-i1+1))
    if verbose:
        print(List)
    return List

def ReadRaw(series, imi, rawdir, crop=None, n0f=3,n0s=3):
    import numpy as np
    from tifffile import imread
    
    # image string index
    imifordir = strindex(imi, n0f)
    
    if len(np.shape(crop))> 0:
        zm=crop[0]
        zM=crop[1]
        ym=crop[2]
        yM=crop[3]
        xm=crop[4]
        xM=crop[5]
    else:
        zm=0
        zM=800
        ym=0
        yM=2016
        xm=0
        xM=2016
        
    # Init image
    image = np.zeros((zM-zm,yM-ym,xM-xm))
    for zi in range(zm, zM):
        # horyzontal slice string index
        zifordir = strindex(zi+1, n0s)
        # horyzontal slice directory
        imdir = rawdir + '/' + series + '/' + 'rec_8bit_phase_'+imifordir + '/' + series+'_'+imifordir+'_'+zifordir+'.rec.8bit.tif'
        # read slice and put into 3D image
        image[zi-zm] = imread(imdir)[ym:yM,xm:xM]
    # return 3D image
    return image

