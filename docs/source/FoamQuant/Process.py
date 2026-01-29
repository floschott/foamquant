def RemoveBackground(image, method='white_tophat', radius=5):
    """
    Remove grey-scale image low frequency background
    
    :param image: 3D image 
    :type image: float numpy array
    :param method: method for removing the background, either 'white_tophat':white tophat filter or 'remove_gaussian': remove the Gaussian filtered image
    :type method: str
    :param radius: white_tophat kernel radius or sigma gaussian filter radius
    :type radius: int
    :return: float numpy array
    """
    
    if method == 'white_tophat':
        from skimage.morphology import white_tophat
        from skimage.morphology import ball
        filtered = white_tophat(image, ball(radius)) # radius
        return filtered
        
    if method == 'remove_gaussian':
        from skimage.filters import gaussian
        import numpy as np
        filtered = gaussian(image, sigma=radius, preserve_range=True) # radius
        filtered = image-filtered
        filtered = (filtered-np.min(filtered))/(np.max(filtered)-np.min(filtered))
        return filtered
    
    
def RemoveSpeckle(image, method='median', radius=1, weight=0.1):
    """
    Remove speckle from the image
    
    :param image: 3D image 
    :type image: float numpy array
    :param method: method for removing the speckle, either 'median', 'gaussian' or 'tv_chambolle'
    :type method: str
    :param radius: kernel radius or sigma gaussian filter radius
    :type radius: int
    :param weight: weight for tv_chambolle
    :type weight: int
    :return: float numpy array
    """
    
    import numpy as np
    if method == 'median':
        from scipy import ndimage
        filtered = ndimage.median_filter(image, size=(radius, radius, radius)) # radius
        return filtered
        
    if method == 'gaussian':
        from skimage.filters import gaussian
        filtered = gaussian(image, sigma=radius) # radius
        return filtered
    
    if method == 'tv_chambolle':
        from skimage.restoration import denoise_tv_chambolle
        filtered = denoise_tv_chambolle(image,weight=weight) # weight
        return filtered


def PhaseSegmentation(image, method='ostu_global', th=0.5, th0=0.3, th1=0.7, returnOtsu=False):
    """
    Perform the phase segmentation
    
    :param image: 3D image 
    :type image: float numpy array
    :param method: Optional, method for segmenting the phase, either 'simple' for simple threshold, 'ostu_global' for a global Otsu threshold, 'niblack', 'sauvola', or 'sobel', Default is 'ostu_global'
    :type method: str
    :param th: Optional, given threshold for 'simple' method
    :type th: float
    :param th0: Optional, given low threshold for 'sobel' method
    :type th0: float
    :param th1: Optional, given high threshold for 'sobel' method
    :type th1: float
    :param returnotsu: Optional, if True, returns Otsu threshold for 'ostu_global' method
    :type returnotsu: Bool
    :return: int numpy array and float
    """
    
    import numpy as np
    
    if method == 'simple':
        segmented = image <= th
        return np.asarray(segmented,dtype='uint8')
        
    if method == 'ostu_global':
        from skimage.filters import threshold_otsu
        t_glob = threshold_otsu(image)
        segmented = image <= t_glob
        if returnOtsu:
            return np.asarray(segmented,dtype='uint8'), t_glob
        return np.asarray(segmented,dtype='uint8')
    
    if method == 'niblack':
        from skimage.filters import threshold_niblack
        t_loc = threshold_niblack(image)
        segmented = image <= t_loc
        return np.asarray(segmented,dtype='uint8')
    
    if method == 'sauvola':
        from skimage.filters import threshold_sauvola
        t_loc = threshold_sauvola(image)
        segmented = image <= t_loc
        return np.asarray(segmented,dtype='uint8')
        
    if method == 'sobel':
        from skimage.filters import sobel
        from skimage.segmentation import watershed
        edges = sobel(image)
        markers = np.zeros_like(image) 
        foreground, background = 1, 2
        markers[image < th0] = background
        markers[image > th1] = foreground
        segmented = watershed(edges, markers)
        segmented = segmented > 1
        return np.asarray(segmented,dtype='uint8')
    
    
def MaskCyl(image, rpercent=None, out=True):
    """ 
    Create a cylindrical mask of the size of the image along the Z axis
    
    :param image: 3D image 
    :type image: float numpy array
    :return: int numpy array
    """    
    
    import numpy as np
    from spam.mesh.structured import createCylindricalMask
    
    # if given radius
    if rpercent != None:
        cyl = createCylindricalMask(np.shape(image), np.int32((np.shape(image)[1]-2)//2*rpercent), voxSize=1.0, centre=None)
    else:
        cyl = createCylindricalMask(np.shape(image), (np.shape(image)[1]-2)//2, voxSize=1.0, centre=None)
    # if middle or outside mask
    if out:
        return cyl
    else:
        return cyl<1


def RemoveSpeckleBin(image, RemoveObjects=True, RemoveHoles=True, BinClosing=False, ClosingRadius=None, GiveVolumes=False, verbose=False, Vminobj=None, Vminhole=None):
    """
    Remove small objects and holes in binary image
    
    :param image: 3D image 
    :type image: int numpy array
    :param RemoveObjects: if True, removes the small objects
    :type RemoveObjects: Bool
    :param RemoveHoles: if True, removes the small holes
    :type RemoveHoles: Bool
    :param BinClosing: if True, perform a binnary closing of radius ClosingRadius
    :type BinClosing: Bool
    :param ClosingRadius: radius of the binnary closing
    :type ClosingRadius: int
    :param GiveVolumes: if True, returns in addition the used min volume thresholds for objects and holes
    :type GiveVolumes: Bool
    :param verbose: if True, print progression steps of the cleaning
    :type verbose: Bool
    :param Vminobj: if given the min volume threshold for the objects is not computed, Vminobj is used as the threshold 
    :type Vminobj: int
    :param Vminhole: if given the min volume threshold for the holes is not computed, Vminobj is used as the threshold 
    :type Vminhole: int
    :return: int numpy array, int, int
    """
    
    import numpy as np
    from skimage.measure import label
    from skimage.measure import regionprops
    from skimage.morphology import remove_small_objects
    
    image = (image > 0)*1
    
    if RemoveObjects:
        if Vminobj == None:
            regions_obj=regionprops(label(image))
            v_obj_beg=[]
            for i in range(len(regions_obj)):
                v_obj_beg.append(regions_obj[i].area)
            if verbose:
                NumberOfObjects_beg = len(regions_obj)
                MaxVolObjects_beg = np.max(v_obj_beg)
                print('Before: Nobj',NumberOfObjects_beg)
            del(regions_obj)
            
            if len(v_obj_beg)>1:
                image = remove_small_objects(label(image), min_size=np.int32(np.max(v_obj_beg)-2))
                if verbose:
                    regions_obj=regionprops(label(image))
                    v_obj_beg=[]
                    for i in range(len(regions_obj)):
                        v_obj_beg.append(regions_obj[i].area)
                    NumberOfObjects_beg = len(regions_obj)
                    print('After: Nobj',NumberOfObjects_beg)
        else:
            if verbose:
                regions_obj=regionprops(label(image))
                v_obj_beg=[]
                for i in range(len(regions_obj)):
                    v_obj_beg.append(regions_obj[i].area)
                NumberOfObjects_beg = len(regions_obj)
                MaxVolObjects_beg = np.max(v_obj_beg)
                print('Before: Nobj',NumberOfObjects_beg)
            image = remove_small_objects(label(image), min_size=Vminobj)
            if verbose:
                regions_obj=regionprops(label(image))
                v_obj_beg=[]
                for i in range(len(regions_obj)):
                    v_obj_beg.append(regions_obj[i].area)
                NumberOfObjects_beg = len(regions_obj)
                print('After: Nobj',NumberOfObjects_beg)
    
    if RemoveHoles:
        if Vminhole == None:
            image = (image < 1)*1
            regions_hol=regionprops(label(image))
            v_hol_beg=[]
            for i in range(len(regions_hol)):
                v_hol_beg.append(regions_hol[i].area)
            if verbose:
                NumberOfHoles_beg = len(regions_hol)
                MaxVolHoles_beg = np.max(v_hol_beg)
                print('Before: Nhol',NumberOfHoles_beg)
            del(regions_hol)
            
            if len(v_hol_beg)>1:
                image = remove_small_objects(label(image), min_size=np.int32(np.max(v_hol_beg)-2))
                if verbose:
                    regions_hol=regionprops(label(image))
                    v_hol_beg=[]
                    for i in range(len(regions_hol)):
                        v_hol_beg.append(regions_hol[i].area)
                    NumberOfHoles_beg = len(regions_hol)
                    print('After: Nhol',NumberOfHoles_beg)
                image = (image < 1)*1
        else:
            image = (image < 1)*1
            if verbose:
                regions_hol=regionprops(label(image))
                v_hol_beg=[]
                for i in range(len(regions_hol)):
                    v_hol_beg.append(regions_hol[i].area)
                NumberOfHoles_beg = len(regions_hol)
                MaxVolHoles_beg = np.max(v_hol_beg)
                print('Before: Nhol',NumberOfHoles_beg)
            image = remove_small_objects(label(image), min_size=Vminobj)
            if verbose:
                regions_hol=regionprops(label(image))
                v_hol_beg=[]
                for i in range(len(regions_hol)):
                    v_hol_beg.append(regions_hol[i].area)
                NumberOfHoles_beg = len(regions_hol)
                print('After: Nhol',NumberOfHoles_beg)
            image = (image < 1)*1
    
    #Bin closing
    if BinClosing:
        from skimage.morphology import closing
        from skimage.morphology import ball
        if closingradius == None:
            image = closing(image)
        else:
            image = closing(image, ball(ClosingRadius))
        if verbose:
            print('Closing done')
    
    image = (image > 0)*1
    
    # If return the threshold volumes for objects and holes
    if GiveVolumes:
        return np.asarray(image, dtype='uint8'), np.int32(np.max(v_obj_beg)-2), np.int32(np.max(v_hol_beg)-2)
    
    return np.asarray(image, dtype='uint8')


def BubbleSegmentation(image, SigSeeds=1, SigWatershed=1, watershed_line=False, radius_opening=None, verbose=False, esti_min_dist=None, compactness=None):
    """
    Perform the bubble watershed segmentation
    
    :param image: 3D image 
    :type image: int numpy array
    :param SigSeeds: Optional, Gaussian filter Sigma for the seeds
    :type SigSeeds: int
    :param SigWatershed: Optional, Gaussian filter Sigma for the watershed
    :type SigWatershed: int
    :param watershed_line: Optional, If True keep the 0-label surfaces between the segmented bubble regions
    :type watershed_line: Bool
    :param radius_opening: Optional, if not None, perform a radius opening operation on the labelled image with the given radius
    :type radius_opening: None or int
    :param verbose: Optional, if True, print progression steps of the segmentation
    :type verbose: Bool
    :param esti_min_dist: Optional, min distance between the seeds
    :type esti_min_dist: None or float
    :param compactness: Optional, compactness of the diffusion
    :type compactness: None or float
    :return: int numpy array
    """
    
    import numpy as np
    from scipy import ndimage
    from skimage.segmentation import watershed
    from skimage.feature import peak_local_max
    from skimage.filters import gaussian
    from skimage.morphology import opening, ball
    
    # Distance map
    image = np.float16(image)
    Distmap = ndimage.distance_transform_edt(image)
    if verbose:
        print('Distance map: done')
    SmoothDistmap = gaussian(Distmap, sigma=SigSeeds);
    if verbose:
        print('Seeds distance map: done')
    
    # Extracting the seeds
    if esti_min_dist != None:
        local_max_coord = peak_local_max(SmoothDistmap, min_distance = int(esti_min_dist), exclude_border=False)
    else:
        local_max_coord = peak_local_max(SmoothDistmap, exclude_border=False)
    local_max_im = np.zeros(np.shape(image))
    for locmax in local_max_coord:
        local_max_im[locmax[0],locmax[1],locmax[2]] = 1
    del(local_max_coord)
    if verbose:
        print('Seeds: done')
    
    # Smooth distance map for watershed
    SmoothDistmap = gaussian(Distmap, sigma=SigWatershed)
    del(Distmap)
    if verbose:
        print('Watershed distance map: done')
    
    # Watershed
    if compactness != None:
        labelled_im = watershed(-SmoothDistmap, ndimage.label(local_max_im)[0], mask = image, compactness = compactness, watershed_line=watershed_line)
    else:
        labelled_im = watershed(-SmoothDistmap, ndimage.label(local_max_im)[0], mask = image, watershed_line=watershed_line)
    del(SmoothDistmap); del(local_max_im)
    if verbose:
        print('Watershed: done')
    
    #Opening
    if radius_opening!=None:
        labelled_im = opening(labelled_im, ball(radius_opening))
        if verbose:
            print('Opening: done')
            
    return labelled_im


def RemoveEdgeBubble(image, mask=None, verbose=False):
    """
    Remove the bubbles on the image edges and in intersection with the masks (if given)
    
    :param image: 3D image 
    :type image: int numpy array
    :param mask: 3D image, if given, removes also the labels at the intersection with the mask
    :type mask: None, Bool or int numpy array
    :param rpercent: If mask is True, will create a cylindrical mask with rpercent of the half-size of the image 
    :type rpercent: float
    :param masktopbottom: list as [zmin, zmax] for top and bottom edges masking
    :type masktopbottom: array
    :param returnmask: if True, aditionally returns the mask
    :type returnmask: Bool
    :param verbose: Optional, if True, print progression steps of the segmentation
    :type verbose: Bool
    :return: int numpy array image or [image, mask] arrays if returnmask is True
    """
    
    from spam.label.label import labelsOnEdges, removeLabels, makeLabelsSequential
    from skimage.measure import regionprops
    #from Package.Process.MaskCyl import MaskCyl
    import numpy as np
    
    # if mask
    if len(np.shape(mask))>0:
        # Masked labels
        imlabtouchmask = image*(1-mask)
        Reg = regionprops(imlabtouchmask)
        labtouchmask=[]
        for reg in Reg:
            labtouchmask.append(reg.label)
        # Remove masked labels
        image = removeLabels(image, labtouchmask)
        # Make sequential
        image = makeLabelsSequential(image)
        if verbose:
            print('Bubbles removed at the mask edges')
    
    # Remove top-bottom edge labels
    labedge = labelsOnEdges(image)
    image = removeLabels(image, labedge)
    image = makeLabelsSequential(image)
    if verbose:
        print('Bubbles removed at the top and bottom edges')
    
    return image


#------------------------------------------------------------

def RemoveBackground_BatchHome(series, rawdir, prossdir, imrange, method='white_tophat', radius=5,  weight=0.1, crop=None, bottom=None, verbose=False, Binning=None, n0rf=3,n0rs=3, n0w=3):
       
    from tifffile import imsave
    from FoamQuant.Helper import ReadRaw
    from FoamQuant.Helper import strindex
    
    
    for imi in imrange:
        image = ReadRaw(series, imi, rawdir, crop=crop, n0f=n0rf, n0s=n0rs)
        
        # image string index
        imifordir = strindex(imi, n0w)
        
        if Binning!=None:
            import spam.DIC
            image = spam.DIC.deform.binning(image, Binning)
            imsave(prossdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_Bin'+str(Binning)+'_'+imifordir, image, bigtiff=True)
        else:
            imsave(prossdir + '/2_RemoveBackground/' + series + '/' + series+'_RemoveBackground_'+imifordir, image, bigtiff=True)
        
        if verbose:
            print(series+' '+str(imi)+': done\n')
            
def RemoveBackground_Batch(nameread, namesave, dirread, dirsave, imrange, method='white_tophat', radius=5, verbose=False, Binning=None, n0=3, endread='.tiff',endsave='.tiff'):
    """
    Uses RemoveBackground function batchwise
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param method: method for removing the background, either 'white_tophat':white tophat filter or 'remove_gaussian': remove the Gaussian filtered image
    :type method: str
    :param radius: white_tophat kernel radius or sigma gaussian filter radius
    :type radius: int
    :param verbose: if True, print progression
    :type verbose: Bool
    :param Binning: the binning number
    :type Binning: None or int
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param endread: read image file extension, default is '.tiff' 
    :type endread: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    """
    
    from tifffile import imread, imsave
    from FoamQuant.Helper import ReadRaw
    from FoamQuant.Helper import strindex
    from FoamQuant.Process import RemoveBackground
    
    for imi in imrange:
        imifordir = strindex(imi, n0)
        # read image
        image = imread(dirread + nameread + imifordir + endread)
        
        if Binning!=None:
            import spam.DIC
            image = spam.DIC.deform.binning(image, Binning)
        
        # remove background
        image = RemoveBackground(image, method=method, radius=radius)
        
        imsave(dirsave + namesave + imifordir + endsave, image, bigtiff=True)
        if verbose:
            print(namesave+' '+str(imi)+': done\n')

def PhaseSegmentation_Batch(nameread, namesave, dirread, dirsave, imrange, method='ostu_global', th=None, th0=None, th1=None, returnOtsu=False, verbose=False, ROIotsu=False, n0=3, endread='.tiff',endsave='.tiff'):
    """
    Uses PhaseSegmentation function batchwise
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param method: Optional, method for segmenting the phase, either 'simple' for simple threshold, 'ostu_global' for a global Otsu threshold, 'niblack', 'sauvola', or 'sobel', Default is 'ostu_global'
    :type method: str
    :param th: Optional, given threshold for 'simple' method
    :type th: float
    :param th0: Optional, given low threshold for 'sobel' method
    :type th0: float
    :param th1: Optional, given high threshold for 'sobel' method
    :type th1: float
    :param returnOtsu: Optional, if True, returns Otsu threshold for 'ostu_global' method
    :type returnOtsu: Bool
    :param verbose: if True, print progression steps of the cleaning
    :type verbose: Bool
    :param ROIotsu: list of length 6 defining the region of interest for determining the single Otsu threshold such as [zmin,zmax,ymin,ymax,xmin,xmax]
    :type ROIotsu: list or numpy array
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param endread: read image file extension, default is '.tiff' 
    :type endread: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    """
    
    from tifffile import imread, imsave
    import numpy as np
    from FoamQuant.Helper import strindex
    
    Lth=[]
    for imi in imrange:
        imifordir = strindex(imi, n0)
        # read image
        image = imread(dirread + nameread + imifordir + endread)
        
        # if keep otsu threshold
        if returnOtsu:
            # region of interest for Otsu
            if len(np.shape(ROIotsu))>0:
                print('ROIotsu', ROIotsu)
                image, th = PhaseSegmentation(image[ROIotsu[0]:ROIotsu[1],ROIotsu[2]:ROIotsu[3],ROIotsu[4]:ROIotsu[5]], 
                                              method='ostu_global',
                                              returnOtsu=returnOtsu)
            else:
                image, th = PhaseSegmentation(image, 
                                              method='ostu_global',
                                              returnOtsu=returnOtsu)
            Lth.append(th)
        else:
            image = PhaseSegmentation(image, method=method, th=th, th0=th0, th1=th1, returnOtsu=False)
            
        # save image
        imsave(dirsave + namesave + imifordir + endsave, np.asarray(image, dtype='uint8'), bigtiff=True)
        # if verbose
        if verbose:
            print(namesave+' '+str(imi)+': done\n')
    
    # return otsu threshold if true
    if returnOtsu:
        return Lth
    
def Masking_Batch(nameread, namesave, dirread, dirsave, imrange, verbose=False, n0=3, endread='.tiff',endsave='.tiff', Masktype=[False,True]):
    """
    Uses MaskCyl function batchwise
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param verbose: if True, print progression steps of the cleaning
    :type verbose: Bool
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param endread: read image file extension, default is '.tiff' 
    :type endread: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    """
    
    from tifffile import imread, imsave
    from FoamQuant.Helper import strindex
    import numpy as np
    
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)
        # read image
        image = imread(dirread + nameread + imifordir + endread)
        
        # if first image, create the mask
        if imi == imrange[0]:
            Mask = np.ones(np.shape(image), dtype='uint8')
            # outside mask
            if  Masktype[1]>0: #want an outside mask
                if Masktype[1]<1: #gave the rpencent
                    Mask = Mask*MaskCyl(image, out=True, rpercent=Masktype[1])
                else: #gave True
                    Mask = Mask*MaskCyl(image, out=True)
            # central mask
            if  Masktype[0]>0: #want an inside mask
                if Masktype[0]<1: #gave the rpencent
                    Mask = Mask*MaskCyl(image, out=False, rpercent=Masktype[0])
                else: #gave True
                    Mask = Mask*MaskCyl(image, out=False)
            
        image = Mask*image
        # save image
        imsave(dirsave + namesave + imifordir + endsave, image, bigtiff=True)
        # if verbose
        if verbose:
            print(namesave+imifordir+': done\n')
            
    
def RemoveSpeckleBin_Batch(nameread, namesave, dirread, dirsave, imrange, verbose=False, RemoveObjects=True, RemoveHoles=True, Cobj=0.5, Chole=0.5,  BinClosing=False, ClosingRadius=None, n0=3,endread='.tiff',endsave='.tiff'):
    """
    Uses RemoveSpeckleBin function batchwise
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param verbose: if True, print the progression
    :type verbose: Bool
    :param RemoveObjects: If True, remove the small objects with a volume below Cobj*MaxObjVol
    :type RemoveObjects: Bool
    :param RemoveHoles: If True, remove the small holes with a volume below Chole*MaxholeVol
    :type RemoveHoles: Bool
    :param Cobj: volume thresholding coefficient for removing the small objects
    :type Cobj: float
    :param Chole: volume thresholding coefficient for removing the small holes
    :type Chole: float
    :param BinClosing: volume thresholding coefficient for removing the small holes
    :type BinClosing: float
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param endread: read image file extension, default is '.tiff' 
    :type endread: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    """
    
    from tifffile import imread, imsave
    from FoamQuant.Helper import strindex
    
    # read first image
    imifordir = strindex(imrange[0], n0)
    image = imread(dirread + nameread + imifordir + endread)  
    image, Vobj, Vhole = RemoveSpeckleBin(image, 
                                          GiveVolumes=True, 
                                          verbose=verbose)
    print('First image (vox): maxObj',Vobj, 'maxHole',Vhole)
    # Volume thesholds
    Vminobj = round(Vobj*Cobj)
    Vminhole = round(Vhole*Chole)
    print('Thresholds (vox): thrObj',Vminobj, 'thrHole',Vminhole,'\n')
    
    for imi in imrange:
        imifordir = strindex(imi, n0)
        # read image
        image = imread(dirread + nameread + imifordir + endread)        
        # remove the small holes and objects
        image = RemoveSpeckleBin(image, 
                                 Vminobj=Vminobj, 
                                 Vminhole=Vminhole, 
                                 verbose=verbose, 
                                 RemoveObjects=RemoveObjects, 
                                 RemoveHoles=RemoveHoles, 
                                 BinClosing=BinClosing, 
                                 ClosingRadius=ClosingRadius)
        # save image
        imsave(dirsave + namesave + imifordir + endsave, image, bigtiff=True)
        # if verbose
        if verbose:
            print(namesave+imifordir+': done\n')
            
            
def BubbleSegmentation_Batch(nameread, namesave, dirread, dirsave, imrange, verbose=False, n0=3, endread='.tiff', endsave='.tiff', writeparameters=False, Binning=None, SigSeeds=1, SigWatershed=1, watershed_line=False, esti_min_dist=None, compactness=None, radius_opening=None, ITK=False, ITKLevel=1):
    """
    Uses BubbleSegmentation function batchwise
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param verbose: if True, print the progression
    :type verbose: Bool
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param endread: read image file extension, default is '.tiff' 
    :type endread: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    :param writeparameters: saved in a text file the segmentation parameters
    :type writeparameters: str
    :param Binning: saved image file extension, default is '.tiff' 
    :type Binning: str
    :param SigSeeds: Optional, Gaussian filter Sigma for the seeds
    :type SigSeeds: int
    :param SigWatershed: Optional, Gaussian filter Sigma for the watershed
    :type SigWatershed: int
    :param watershed_line: Optional, If True keep the 0-label surfaces between the segmented bubble regions
    :type watershed_line: Bool
    :param esti_min_dist: Optional, min distance between the seeds
    :type esti_min_dist: None or float
    :param compactness: Optional, compactness of the diffusion
    :type compactness: None or float
    :param radius_opening: Optional, if not None, perform a radius opening operation on the labelled image with the given radius
    :type radius_opening: None or int
    :param ITK: If True, the ITKwatershed from SPAM is used
    :type ITK: Bool
    :param ITKLevel: Optional, default is 1
    :type ITKLevel: float
    """
    
    from tifffile import imread, imsave
    import os
    from FoamQuant.Helper import strindex
    
    #Check saving directory
    isExist = os.path.exists(dirsave)
    print('Path exist:', isExist)
    if not isExist:
        print('Error: saving path does not exist', dirsave)
        return
    if ITK:
        import spam.label
    else:
        # write parameters
        if writeparameters:
            file1 = open(dirsave + namesave + '_WatershedParameters.txt',"w")
            L = ["nameread \n",str(nameread),
                 "\n imrange \n",str(imrange),
                 "\n SigSeeds \n",str(SigSeeds),
                 "\n SigWatershed \n",str(SigWatershed),
                 "\n watershed_line \n",str(watershed_line),
                 "\n radius_opening \n",str(radius_opening),
                 "\n Twatershed_line \n",str(watershed_line),
                 "\n radius_opening \n",str(radius_opening),
                 "\n esti_min_dist \n",str(esti_min_dist),
                 "\n compactness \n",str(compactness)] 
            file1.writelines(L)
            file1.close() 

    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)
        # read image
        image = imread(dirread + nameread + imifordir + endread)
        # if binning
        if Binning != None:
            import spam.DIC
            image = spam.DIC.deform.binning(image, Binning)
        # Bubble segmentation
        if ITK:
            image = spam.label.ITKwatershed.watershed(image, markers=None, watershedLevel=ITKLevel)
        else:
            SigSeeds = SigSeeds//Binning
            SigWatershed = SigWatershed//Binning
            esti_min_dist = esti_min_dist//Binning
            image = BubbleSegmentation(image, 
                                       SigSeeds, 
                                       SigWatershed, 
                                       watershed_line, 
                                       radius_opening, 
                                       verbose, 
                                       esti_min_dist, 
                                       compactness)
        # save image
        imsave(dirsave + namesave + imifordir + endsave, image, bigtiff=True)
        # if verbose
        if verbose:
            print(namesave+imifordir+': done\n')
            
            
def RemoveEdgeBubble_Batch(nameread, namesave, dirread, dirsave, imrange, verbose=False, endread='.tiff', endsave='.tiff', n0=3, masktopbottom=None, Masktype=[False,True]):
    """
    Uses RemoveEdgeBubble function batchwise
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param verbose: if True, print the progression
    :type verbose: Bool
    :param endread: read image file extension, default is '.tiff' 
    :type endread: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param maskcyl: if True, create a cylindrical mask
    :type maskcyl: Bool
    :param rpercent: If mask is True, will create a cylindrical mask with rpercent of the half-size of the image 
    :type rpercent: float
    :param masktopbottom: list as [zmin, zmax] for top and bottom edges masking
    :type masktopbottom: array
    """
    
    from tifffile import imread, imsave
    #from Package.Process.RemoveEdgeBubble import RemoveEdgeBubble
    import os
    from FoamQuant.Helper import strindex
    import numpy as np
    
    #Check directory
    isExist = os.path.exists(dirsave)
    print('Path exist:', isExist)
    if not isExist:
        print('Error: Saving path does not exist', dirsave)
        return
        
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)
        # read image
        image = imread(dirread + nameread + imifordir + endread)
        
        
        # if the first image, create the mask
        if imi == imrange[0]:
            if Masktype[1]>0 or Masktype[0]>0 or masktopbottom != None: #if should have a mask?
                Mask=np.ones(np.shape(image), dtype='uint8')
                # outside mask
                if  Masktype[1]>0: #want an outside mask
                    if Masktype[1]<1: #gave the rpencent
                        Mask = Mask*MaskCyl(image, out=True, rpercent=Masktype[1])
                    else: #gave True
                        Mask = Mask*MaskCyl(image, out=True)
                    if verbose:
                        print('Outside mask was created')
                # central mask
                if  Masktype[0]>0: #want an inside mask
                    if Masktype[0]<1: #gave the rpencent
                        Mask = Mask*MaskCyl(image, out=False, rpercent=Masktype[0])
                    else: #gave True
                        Mask = Mask*MaskCyl(image, out=False)
                    if verbose:
                        print('Central mask was created')
                # top-bottom mask
                if masktopbottom != None:
                    Mask[:masktopbottom[0]] = 0
                    Mask[masktopbottom[1]:] = 0
                    if verbose:
                        print('top-bottom edges mask was created')
        
        if Masktype[1]>0 or Masktype[0]>0 or masktopbottom != None: #if should have a mask?
            image = RemoveEdgeBubble(image, mask=Mask) 
        else:
            image = RemoveEdgeBubble(image)
            
        # save image
        imsave(dirsave + namesave + imifordir + endsave, image, bigtiff=True)

        # if verbose
        if verbose:
            print(namesave+imifordir+': done\n')
            
def GetContacts_Batch(nameread, nameread_noedge, namesave, dirread, dirread_noedge, dirsave, imrange, verbose=False, endread='.tiff',endread_noedge='.tiff',endsave='.tiff', n0=3, save='all', maximumCoordinationNumber=20):
    """
    Uses GetContacts function batchwise
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param verbose: if True, print the progression
    :type verbose: Bool
    :param endread: read labeled image file extension, default is '.tiff' 
    :type endread: str
    :param endread_noedge: read labeled no-edge image file extension, default is '.tiff' 
    :type endread_noedge: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param save: 'all': save all, 'coord':save the coordination image, 'image': save contact image, table': save the contact table or 'pair': save the contact pairs
    :type save: str or list of str
    :param maximumCoordinationNumber: the maximum coordination number, default 20
    :type maximumCoordinationNumber: int
    """
    
    from tifffile import imread, imsave
    import os
    from FoamQuant.Helper import strindex
    import csv
    
    #Check directory
    isExist = os.path.exists(dirread)
    print('Path exist:', isExist)
    if not isExist:
        print('Error: Saving path does not exist', dirread)
        return
        
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)
        # read images
        image = imread(dirread + nameread + imifordir + endread)
        image_noedge = imread(dirread_noedge + nameread_noedge + imifordir + endread_noedge)
        
        # GetContacts 
        if save == 'all' or 'coord' in save:
            lab, labnoedge_ext, centroid, image,Z,contactTable,contactingLabels, CoordinationImage = GetContacts(image, 
                                                                                                  image_noedge, 
                                                                                                  maximumCoordinationNumber=maximumCoordinationNumber, 
                                                                                                  returnCoordImage=True)
            # Coordination image
            imsave(dirsave + 'Coordination_' + imifordir + endsave, CoordinationImage, bigtiff=True)
        else:
            
            lab, labnoedge_ext, centroid, image,Z,contactTable,contactingLabels = GetContacts(image, 
                                                                               image_noedge, 
                                                                               maximumCoordinationNumber=maximumCoordinationNumber, 
                                                                               returnCoordImage=False)            
        # save contact image
        if save=='all' or 'image' in save:
            # save image
            imsave(dirsave + namesave + imifordir + endsave, image, bigtiff=True)
        # save contact table
        if save=='all' or 'table' in save:
            with open(dirsave + namesave+ 'table_' + imifordir + '.tsv', 'w', newline='') as csvfile:        
                writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

                #first line
                line1 = ['lab','lab_noedge','Z','z','y','x']
                for labcontacti in range(1, maximumCoordinationNumber+1):
                    line1.append('lab'+str(labcontacti))
                    line1.append('cont'+str(labcontacti))
                writer.writerow(line1)

                #data line by line
                for i in range(len(Z)):
                    linei = [i+1, labnoedge_ext[i], Z[i],centroid[i][0],centroid[i][1],centroid[i][2]]
                    for labcontacti in range(maximumCoordinationNumber):
                        linei.append(contactTable[i][2*labcontacti])
                        linei.append(contactTable[i][2*labcontacti+1])
                    writer.writerow(linei)
                    
        # save pair table
        if save=='all' or save=='pair':
            with open(dirsave + namesave + 'pair_' + imifordir + '.tsv', 'w', newline='') as csvfile:        
                writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

                #first line
                writer.writerow(['cont','lab1','lab2'])#,'z1','z2','y1','y2','x1','x2'])
                #data line by line
                for i in range(1,len(contactingLabels)):
                    writer.writerow([i, contactingLabels[i-1][0],contactingLabels[i-1][1]])#, 
                                     #contactingLabels[i-1][2],contactingLabels[i-1][3],
                                     #contactingLabels[i-1][4],contactingLabels[i-1][5],
                                     #contactingLabels[i-1][6],contactingLabels[i-1][7]])
        # if verbose
        if verbose:
            print(namesave+imifordir+': done')
            
            
            
def GetContacts(image, image_noedge, maximumCoordinationNumber=20, returnCoordImage=False):
    """
    Return labels [0], centroids [1], coordinations [2] of the no-edge image, (and coordination image [3])
    
    :param image: full 3D image 
    :type image: int numpy array
    :param image_noedge: 3D image with removed label at the edges
    :type image_noedge: int numpy array
    :param maximumCoordinationNumber: the maximum coordination number, default 20
    :type maximumCoordinationNumber: int
    :param returnCoordinationImage: if True, additionally returns image_noedge coordination image
    :type returnCoordinationImage: Bool
    :return: labels [0], centroids [1], coordinations [2] of the no-edge image, (and coordination image [3])
    """    
    
    import numpy as np
    from skimage.measure import regionprops
    import spam.label
    
    # Regions Im edge
    reg = regionprops(image)
    centroid = []; lab = []
    for reg in reg:
        lab.append(reg.label)
        centroid.append(reg.centroid)
    # Regions Im noedge
    reg_noedge = regionprops(image_noedge)
    centroid_noedge = []; lab_noedge=[]
    for reg in reg_noedge:
        lab_noedge.append(reg.label)
        centroid_noedge.append(reg.centroid)

    # Contacts Im edge
    contactVolume, Z, contactTable, contactingLabels = spam.label.contacts.labelledContacts(image, maximumCoordinationNumber = maximumCoordinationNumber)
    del(contactingLabels)
    
    # Contacts Im no-edge
    contactVolume_ned, Z_ned, contactTable_ned, contactingLabels_ned = spam.label.contacts.labelledContacts(image_noedge, maximumCoordinationNumber = maximumCoordinationNumber)
    del(contactVolume_ned); del(Z_ned); del(contactTable_ned)
    
    Z = Z[1:] #remove background label (0 label)
    contactTable = contactTable[1:] #remove background label (0 label)
    
    #lengthcontactingLabels = len(contactingLabels)
    #CoordcontactingLabels = np.zeros((lengthcontactingLabels,8))
    #CoordcontactingLabels[:,0:2] = contactingLabels
    
    #for tabi in range(lengthcontactingLabels):
    #    for labi in range(len(lab)):
    #        if lab[labi] == contactingLabels[tabi][0]:
    #            CoordcontactingLabels[tabi][2] = centroid[labi][0]
    #            CoordcontactingLabels[tabi][3] = centroid[labi][1]
    #            CoordcontactingLabels[tabi][4] = centroid[labi][2]
    #        if lab[labi] == contactingLabels[tabi][1]:
    #            CoordcontactingLabels[tabi][5] = centroid[labi][0]
    #            CoordcontactingLabels[tabi][6] = centroid[labi][1]
    #            CoordcontactingLabels[tabi][7] = centroid[labi][2]
        
    # Keep Coordination at the edge
    Z_noegde=[0]
    labnoedge_ext=[]
    for h in range(len(Z)):
        z,y,x = centroid[h]                     
        z = round(z); y = round(y); x = round(x)
        if image_noedge[z,y,x] != 0: 
            Z_noegde.append(Z[h])
            labnoedge_ext.append(image_noedge[z,y,x])
        else:
            labnoedge_ext.append(-1)

    if returnCoordImage:
        CoordinationImage = spam.label.label.convertLabelToFloat(image_noedge, np.asarray(Z_noegde))
        return lab, labnoedge_ext, centroid, contactVolume,Z,contactTable,contactingLabels_ned, CoordinationImage
    
    return lab, labnoedge_ext, centroid, contactVolume,Z,contactTable,contactingLabels_ned




# -------------------------------

def FastLocalThickness_Batch(nameread, namesave, dirread, dirsave, imrange, verbose=False, endread='.tiff', endsave='.tiff', n0=3, WalThickness=True, Separation=False, scale=1):
    """
    Uses localthickness function batchwise. IMPORTANT, please refer to Dahl, V. A. and Dahl A. B. work: Git-link February 2023: https://github.com/vedranaa/local-thickness.git
    
    :param nameread: read image name 
    :type nameread: str
    :param namesave: saved image name 
    :type namesave: str
    :param dirread: read image directory 
    :type dirread: str
    :param dirread: saved image directory 
    :type dirread: str
    :param imrange: image indexes array
    :type imrange: list or numpy array
    :param verbose: if True, print the progression
    :type verbose: Bool
    :param endread: read image file extension, default is '.tiff' 
    :type endread: str
    :param endsave: saved image file extension, default is '.tiff' 
    :type endsave: str
    :param n0: number of digit for the saving index, default is 3
    :type n0: int
    :param WalThickness: If True, save the wall thickness (zeros in the binary image)
    :type WalThickness: Bool
    :param Separation: If True, save the separation thickess (ones in the binary image)
    :type Separation: Bool
    :param scale: Optional downscaling factor, default is 1
    :type scale: float
    """
    
    from tifffile import imread, imsave
    import localthickness as lt
    import os
    from FoamQuant.Helper import strindex
    
    #Check directory
    isExist = os.path.exists(dirsave)
    print('Path exist:', isExist)
    if not isExist:
        print('Error: Saving path does not exist', dirsave)
        return
        
    for imi in imrange:
        # image string index
        imifordir = strindex(imi, n0)
        # read image
        image = imread(dirread + nameread + imifordir + endread)
        # Wall thickness
        if WalThickness:
            thickness = lt.local_thickness(image<1, scale=scale)
            imsave(dirsave+namesave+imifordir+'_WT'+endsave, thickness, bigtiff=True)
        # Separation
        if Separation:
            separation = lt.local_thickness(image>0, scale=scale)
            imsave(dirsave+namesave+imifordir+'_SEP'+endsave, separation, bigtiff=True)
        # if verbose
        if verbose:
            print(namesave+imifordir+': done\n')
            
            
