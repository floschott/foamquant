def RemoveBackground(image, method='white_tophat', radius=5):
    
    ''' Remove grey-scale image low frequency background
    
    :param image: 3D image 
    :type image: float array
    
    :param method: method for removing the background, either 'white_tophat':white tophat filter or 'remove_gaussian': remove the Gaussian filtered image
    :type method: str
    
    :param radius: white_tophat kernel radius or sigma gaussian filter radius
    :type radius: int

    '''    
    
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


def RemoveSpeckle(image, method='median', radius=5, weight=0.1):
    import numpy as np
    if method == 'median':
        from scipy import ndimage
        filtered = ndimage.median_filter(crop, size=(radius, radius, radius)) # radius
        return filtered
        
    if method == 'gaussian':
        from skimage.filters import gaussian
        filtered = gaussian(image, sigma=radius) # radius
        return filtered
    
    if method == 'tv_chambolle':
        from skimage.restoration import denoise_tv_chambolle
        filtered = denoise_tv_chambolle(image,weight=weight) # weight
        return filtered

def PhaseSegmentation(image, method='ostu_global', th=0.5, radius=5, th0=0.3, th1=0.7, returnotsu=False):
    import numpy as np
    
    if method == 'simple':
        segmented = image <= th
        return np.asarray(segmented,dtype='uint8')
        
    if method == 'ostu_global':
        from skimage.filters import threshold_otsu
        t_glob = threshold_otsu(image)
        segmented = image <= t_glob
        if returnotsu:
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

def MaskCyl(image):
    import numpy as np
    from spam.mesh.structured import createCylindricalMask
    cyl = createCylindricalMask(np.shape(image), (np.shape(image)[1]-2)//2, voxSize=1.0, centre=None)
    return cyl
  
def RemoveSpeckleBin(image, RemoveObjects=True, RemoveHoles=True, BinClosing=False, ClosingRadius=None, GiveVolumes=False, Verbose=True, Vminobj=None, Vminhole=None):
    import numpy as np
    from skimage.measure import label
    from skimage.measure import regionprops
    
    image = np.asarray(image,dtype='uint8')
    
    if RemoveObjects:
        if Vminobj == None:
            #Remove small objects
            regions_obj=regionprops(label(image))
            v_obj_beg=[]
            for i in range(len(regions_obj)):
                v_obj_beg.append(regions_obj[i].area)
            NumberOfObjects_beg = len(regions_obj)
            MaxVolObjects_beg = np.max(v_obj_beg)
            del(regions_obj)
            if len(v_obj_beg)>1:
                from skimage.morphology import remove_small_objects
                image = remove_small_objects(label(image), min_size=np.int(np.max(v_obj_beg)-2))
                if Verbose:
                    print('Small object removed')
        else:
            #Remove small objects
            from skimage.morphology import remove_small_objects
            image = remove_small_objects(label(image), min_size=Vminobj)
            if Verbose:
                    print('Small object removed')
    
    if RemoveHoles:
        if Vminhole == None:
        #Remove small holes
            regions_hol=regionprops(label(1-image))
            v_hol_beg=[]
            for i in range(len(regions_hol)):
                v_hol_beg.append(regions_hol[i].area)
            NumberOfHoles_beg = len(regions_hol)
            MaxVolHoles_beg = np.max(v_hol_beg)
            del(regions_hol)
            if len(v_hol_beg)>1:
                from skimage.morphology import remove_small_objects
                image = 1-remove_small_objects(label(1-image), min_size=np.int(np.max(v_hol_beg)-2))
                if Verbose:    
                    print('Small holes removed')
        else:
            from skimage.morphology import remove_small_objects
            image = 1-remove_small_objects(label(1-image), min_size=Vminhole)
            if Verbose:    
                print('Small holes removed')
    
    #Bin closing
    if BinClosing:
        from skimage.morphology import closing
        from skimage.morphology import ball
        if closingradius == None:
            image = closing(image)
        else:
            image = closing(image, ball(ClosingRadius))
        if Verbose:
            print('Closing done')
    
    image = (image > 0)*1
    
    # If return the threshold volumes for objects and holes
    if GiveVolumes:
        return np.asarray(image, dtype='uint8'), np.int(np.max(v_obj_beg)-2), np.int(np.max(v_hol_beg)-2)
    
    return np.asarray(image, dtype='uint8')

def BubbleSegmentation(image, SigSeeds=1, SigWatershed=1, watershed_line=False, radius_opening=None, verbose=False):
    import numpy as np
    from scipy import ndimage as ndi
    from skimage.segmentation import watershed
    from skimage.feature import peak_local_max
    from skimage.filters import gaussian
    from skimage.morphology import opening, ball
    
    #Marker seeds
    distance = gaussian(ndi.distance_transform_edt(image),sigma=SigSeeds)
    if verbose:
        print('seeds distance map: done')
    coords = peak_local_max(distance, labels=image)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)
    del(mask)
    if verbose:
        print('seeds: done')
        
    #Distance map
    distance = gaussian(ndi.distance_transform_edt(image),sigma=SigWatershed)
    if verbose:
        print('watershed distance map: done')
    labels = watershed(-distance, markers, mask=image,watershed_line=watershed_line)
    del(distance); del(markers)
    if verbose:
        print('watershed: done')
    
    #Opening
    if radius_opening!=None:
        labels = opening(labels, ball(radius_opening))
        if verbose:
            print('opening: done')
            
    return labels

def RemoveEdgeBubble(image, mask=None):
    from spam.label.label import labelsOnEdges, removeLabels, makeLabelsSequential
    from skimage.measure import regionprops
    
    # If: a mask is given
    if mask != None:
        # Masked labels
        imlabtouchmask = image*(1-mask)
        Reg = regionprops(imlabtouchmask)
        labtouchmask=[]
        for reg in Reg:
            labtouchmask.append(reg.label)
        # Remove masked labels
        nomasked = removeLabels(image, Ledge)
        # Make sequential
        nomasked = makeLabelsSequential(nomasked)
        
        return nomasked
    
    # Else: Remove edge labels
    labedge = labelsOnEdges(image)
    noedge = removeLabels(image, labedge)
    noedge = makeLabelsSequential(noedge)
    
    return noedge
    
    
      
      
      
      
      
