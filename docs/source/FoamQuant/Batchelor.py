# COPIED and MODIFIED FUNCTIONS FROM PoreSpy and SCIKIT-IMAGE

# Function 1: Region surface areas - COPIED
# Source: https://porespy.org/_modules/porespy/metrics/_meshtools.html#region_surface_areas
def region_surface_areas(regions, voxel_size=1, strel=None): 
    """
    Compute the surface area of each labeled region in a 2D or 3D image.

    Optionally computes interfacial areas by meshing each region
    using the marching cubes algorithm.

    :param regions: Labeled image where each integer label defines a region
                    (label 0 is ignored)
    :type regions: numpy.ndarray
    :param voxel_size: Physical size of one voxel edge
    :type voxel_size: float
    :param strel: Structuring element used to blur the region prior to meshing
    :type strel: numpy.ndarray or None

    :return: Surface area of each region (index = label - 1)
    :rtype: numpy.ndarray
    """
    
    #logger.info('Finding surface area of each region')
    im = regions
    if strel is None:
        strel = ps_round(1, im.ndim, smooth=False)
    # Get 'slices' into im for each pore region
    slices = spim.find_objects(im)
    # Initialize arrays
    Ps = np.arange(1, np.amax(im) + 1)
    sa = np.zeros_like(Ps, dtype=float)
    # Start extracting marching cube area from im
    msg = "Computing region surface area".ljust(60)
    for i in Ps:
        reg = i - 1
        if slices[reg] is not None:
            s = extend_slice(slices[reg], im.shape)
            sub_im = im[s]
            mask_im = sub_im == i
            mesh = mesh_region(region=mask_im, strel=strel)
            sa[reg] = mesh_surface_area1(mesh)
    result = sa * voxel_size**2
    return result

# Function 2: ps_round - COPIED
# Source: https://porespy.org/_modules/porespy/tools/_funcs.html#ps_round
def ps_round(r, ndim, smooth=True):
    """
    Create a spherical (or circular) structuring element.

    :param r: Radius of the structuring element
    :type r: float
    :param ndim: Dimensionality of the element (2 or 3)
    :type ndim: int
    :param smooth: If True, use smooth sphere faces
    :type smooth: bool

    :return: Structuring element
    :rtype: numpy.ndarray
    """
    
    rad = int(np.ceil(r))
    other = np.ones([2*rad + 1 for i in range(ndim)], dtype=bool)
    other[tuple(rad for i in range(ndim))] = False
    if smooth:
        ball = edt(other) < r
    else:
        ball = edt(other) <= r
    return ball

# Function 3: extend_slice - COPIED
# Source: https://porespy.org/_modules/porespy/tools/_funcs.html#extend_slice
def extend_slice(slices, shape, pad=1):
    """
    Extend slice indices by a fixed padding while respecting image bounds.

    :param slices: Slice objects defining the region of interest
    :type slices: tuple of slice
    :param shape: Shape of the full image
    :type shape: tuple or list
    :param pad: Number of voxels to extend in each direction
    :type pad: int or array-like

    :return: Extended slice objects
    :rtype: tuple of slice
    """
    
    shape = np.array(shape)
    pad = np.array(pad).astype(int)*(shape > 0)
    a = []
    for i, s in enumerate(slices):
        start = 0
        stop = shape[i]
        start = max(s.start - pad[i], 0)
        stop = min(s.stop + pad[i], shape[i])
        a.append(slice(start, stop, None))
    return tuple(a)


# Function 3: mesh_region - COPIED
# Source: https://porespy.org/_modules/porespy/tools/_funcs.html#mesh_region
def mesh_region(region: bool, strel=None):
    """
    Generate a triangular surface mesh of a region using marching cubes.

    :param region: Boolean image defining the region of interest
    :type region: numpy.ndarray
    :param strel: Structuring element used to blur the region prior to meshing
    :type strel: numpy.ndarray or None

    :return: Mesh object containing vertices, faces, normals and values
    :rtype: Results
    """
    
    im = region
    _check_for_singleton_axes(im)
    if strel is None:
        if region.ndim == 3:
            strel = ball(1)
        if region.ndim == 2:
            strel = disk(1)
    pad_width = np.amax(strel.shape)
    if im.ndim == 3:
        padded_mask = np.pad(im, pad_width=pad_width, mode='constant')
        padded_mask = spim.convolve(padded_mask * 1.0,
                                    weights=strel) / np.sum(strel)
    else:
        padded_mask = np.reshape(im, (1,) + im.shape)
        padded_mask = np.pad(padded_mask, pad_width=pad_width, mode='constant')
    verts, faces, norm, val = marching_cubes(padded_mask)
    result = Results()
    result.verts = verts - pad_width
    result.faces = faces
    result.norm = norm
    result.val = val
    return result

# Function 4: _check_for_singleton_axes - COPIED
# Source: https://porespy.org/_modules/porespy/tools/_funcs.html#mesh_region
def _check_for_singleton_axes(im):  # pragma: no cover
    """
    Check for singleton dimensions in an image and issue a warning.

    :param im: Input image
    :type im: numpy.ndarray

    :return: None
    """
    
    if im.ndim != im.squeeze().ndim:
        logger.warning("Input image conains a singleton axis. Reduce"
                       " dimensionality with np.squeeze(im) to avoid"
                       " unexpected behavior.")

# Class: COPIED
# Source: https://porespy.org/_modules/porespy/tools/_utils.html#Results
class Results:
    r"""
    A minimal class for use when returning multiple values from a function

    This class supports dict-like assignment and retrieval
    (``obj['im'] = im``), namedtuple-like attribute look-ups (``obj.im``),
    and generic class-like object assignment (``obj.im = im``)

    """

    def __init__(self, **kwargs):
        self._func = inspect.getouterframes(inspect.currentframe())[1].function
        self._time = time.asctime()

    def __iter__(self):
        for k, v in self.__dict__.items():
            if not k.startswith('_'):
                yield v

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __str__(self):
        header = "―" * 78
        lines = [
            header,
            f"Results of {self._func} generated at {self._time}",
            header,
        ]
        for item in list(self.__dict__.keys()):
            if item.startswith('_'):
                continue
            if (isinstance(self[item], np.ndarray)):
                s = np.shape(self[item])
                lines.append("{0:<25s} Array of size {1}".format(item, s))
            elif hasattr(self[item], 'keys'):
                N = len(self[item].keys())
                lines.append("{0:<25s} Dictionary with {1} items".format(item, N))
            else:
                lines.append("{0:<25s} {1}".format(item, self[item]))
        lines.append(header)
        return "\n".join(lines)

    
# Function 5: mesh_surface_area - COPIED
# Source: https://porespy.org/_modules/porespy/metrics/_meshtools.html#mesh_surface_area
def mesh_surface_area1(mesh=None, verts=None, faces=None):
    """
    Compute the surface area of a triangular mesh.

    :param mesh: Mesh object returned by ``mesh_region``
    :type mesh: Results or None
    :param verts: Mesh vertices
    :type verts: numpy.ndarray or None
    :param faces: Mesh faces
    :type faces: numpy.ndarray or None

    :return: Surface area of the mesh
    :rtype: float
    """
    if mesh:
        verts = mesh.verts
        faces = mesh.faces
    else:
        if (verts is None) or (faces is None):
            raise Exception('Either mesh or verts and faces must be given')
    surface_area = measure.mesh_surface_area(verts, faces)
    return surface_area

# Function 6: show_mesh - COPIED
# Source: https://porespy.org/_modules/porespy/visualization/_plots.html#show_mesh
def show_mesh(mesh):  # pragma: no cover
    """
    Visualize a triangular surface mesh in 3D.

    :param mesh: Mesh object returned by ``marching_cubes`` or ``mesh_region``
    :type mesh: Results

    :return: Matplotlib figure handle
    :rtype: matplotlib.figure.Figure
    """
    try:
        verts = mesh.vertices
    except AttributeError:
        verts = mesh.verts

    lim_max = np.amax(verts, axis=0)
    lim_min = np.amin(verts, axis=0)

    # Display resulting triangular mesh using Matplotlib.
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[mesh.faces])
    mesh.set_edgecolor('k')

    ax.add_collection3d(mesh)
    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")
    ax.set_xlim(lim_min[0], lim_max[0])
    ax.set_ylim(lim_min[1], lim_max[1])
    ax.set_zlim(lim_min[2], lim_max[2])

    return fig





#///////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////

# Function 7: Region_surface_areas MODIFIED as region_Batchelor_tensors
# Source: https://porespy.org/_modules/porespy/metrics/_meshtools.html#region_surface_areas
def region_Batchelor_tensors(regions, voxel_size=1, strel=None):
    """
    Compute the Batchelor tensor for each labeled region in an image.

    :param regions: Labeled image where each integer label defines a region
    :type regions: numpy.ndarray
    :param voxel_size: Physical voxel size
    :type voxel_size: float
    :param strel: Structuring element used to blur regions before meshing
    :type strel: numpy.ndarray or None

    :return:
        - Batchelor tensor for each region (N, 3, 3)
        - Surface area of each region
    :rtype: tuple (numpy.ndarray, numpy.ndarray)
    """
    from tqdm import tqdm
    im = regions
    if strel is None:
        strel = ps_round(1, im.ndim, smooth=False)
    # Get 'slices' into im for each pore region
    slices = spim.find_objects(im)
    # Initialize arrays
    Ps = np.arange(1, np.amax(im) + 1)
    B = np.zeros((len(Ps),3,3), dtype=float)          # Additionally the list of tensors
    Surf = np.zeros((len(Ps)), dtype=float)           # Modified as list of surfaces
    # Start extracting marching cube area from im
    msg = "Computing region surface area".ljust(60)
    for i in tqdm(Ps):
        reg = i-1
        if slices[reg] is not None:
            s = extend_slice(slices[reg], im.shape)
            sub_im = im[s]
            mask_im = sub_im == i
            mesh = mesh_region(region=mask_im, strel=strel) 
            B[reg], Surf[reg] = mesh_Batchelor_tensor(mesh)     # Modified wth Batchelor function
    return B, Surf

# Function 8: mesh_surface_area1 MODIFIED as mesh_Batchelor_tensor
# Source: https://porespy.org/_modules/porespy/metrics/_meshtools.html#region_surface_area
def mesh_Batchelor_tensor(mesh=None, verts=None, faces=None):
    """
    Compute the Batchelor tensor from a surface mesh.

    :param mesh: Mesh object returned by ``mesh_region``
    :type mesh: Results or None
    :param verts: Mesh vertices
    :type verts: numpy.ndarray or None
    :param faces: Mesh faces
    :type faces: numpy.ndarray or None

    :return:
        - Batchelor tensor (3x3)
        - Surface area of the mesh
    :rtype: tuple (numpy.ndarray, float)
    """
    if mesh:
        verts = mesh.verts
        faces = mesh.faces
    else:
        if (verts is None) or (faces is None):
            raise Exception('Either mesh or verts and faces must be given')
    mesh_for_batchelor = mesh_batchelor(verts, faces)
    return mesh_for_batchelor

# Function 9: skimage.measure.mesh_surface_area MODIFIED as mesh_batchelor
# Source: https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.mesh_surface_area
def mesh_batchelor(verts, faces):
    """
    Compute the Batchelor tensor from triangular mesh geometry.

    :param verts: Vertex coordinates (V, 3)
    :type verts: numpy.ndarray
    :param faces: Triangle indices (F, 3)
    :type faces: numpy.ndarray

    :return:
        - Batchelor tensor (3x3)
        - Surface area
    :rtype: tuple (numpy.ndarray, float)
    """
    # Fancy indexing to define two vector arrays from triangle vertices
    actual_verts = verts[faces]
    a = actual_verts[:, 0, :] - actual_verts[:, 1, :]
    b = actual_verts[:, 0, :] - actual_verts[:, 2, :]
    del actual_verts

    # Area of triangle in 3D = 1/2 * Euclidean norm of cross product
    #return ((np.cross(a, b) ** 2).sum(axis=1) ** 0.5).sum() / 2.

    #Area of triangle in 3D = 1/2 * Euclidean norm of cross product
    ds = (np.cross(a, b) ** 2).sum(axis=1) ** 0.5 / 2
    
    #Normal vector to the surface in 3D = cross product
    w = np.cross(a, b)    
    #local Batchelor tensor (deviatoric part)
    L = len(w)
    listB=np.zeros((L,3,3), dtype=float)
    for i in range(L):
        norm = np.linalg.norm(w[i])
        if np.linalg.norm(w[i])>0:
            ni=w[i]/np.linalg.norm(w[i])
        else:
            ni=w[i]
        N2=[[ni[0]**2,ni[0]*ni[1],ni[0]*ni[2]],
            [ni[0]*ni[1],ni[1]**2,ni[1]*ni[2]],
            [ni[0]*ni[2],ni[1]*ni[2],ni[2]**2]]
        listB[i] = (np.identity(3) - N2) * ds[i]
    B = listB.sum(axis=0)
    return B, np.sum(ds, axis=0)


# Batch Batchelor .tsv saving ------------------------------
def Batchelor_Batch(nameread, namesave, dirread, dirsave, imrange, verbose=False, endread='.tif', endsave='.tsv', n0=3):
    """
    Compute Batchelor tensors for a batch of labeled images and save to TSV.

    :param nameread: Base name of input images
    :type nameread: str
    :param namesave: Base name for output files
    :type namesave: str
    :param dirread: Directory containing labeled images
    :type dirread: str
    :param dirsave: Directory to save TSV files
    :type dirsave: str
    :param imrange: Image indices to process
    :type imrange: list or array
    :param verbose: Print progress information
    :type verbose: bool
    :param endread: Input image file extension
    :type endread: str
    :param endsave: Output file extension
    :type endsave: str
    :param n0: Zero-padding for image indices
    :type n0: int

    :return: None
    """ 
    
    import numpy as np 
    from tifffile import imread
    import csv
    from FoamQuant.FromLabelled import RegionProp
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
        # Region properties->lab, coord, volume
        Prop = RegionProp(image, field=False)
        # Batchelor tensor
        lB, lSurf = region_Batchelor_tensors(image, voxel_size=1, strel=None)
        
        
        
        # Save as TSV
        with open(dirsave + namesave + imifordir + endsave, 'w', newline='') as csvfile:        
            writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['lab','z','y','x', 'vol',
                             'mesharea', 
                             'B11','B12','B13','B22','B23','B33','b11','b12','b13','b22','b23','b33'])
            for i in range(len(lB)):
                writer.writerow([Prop[0][i], Prop[1][i][0],Prop[1][i][1],Prop[1][i][2], Prop[2][i],
                                 lSurf[i],
                                 lB[i][0][0], lB[i][0][1], lB[i][0][2], lB[i][1][1], lB[i][1][2], lB[i][2][2],
                                 lB[i][0][0]/Prop[2][i], lB[i][0][1]/Prop[2][i], lB[i][0][2]/Prop[2][i], 
                                 lB[i][1][1]/Prop[2][i], lB[i][1][2]/Prop[2][i], lB[i][2][2]/Prop[2][i]])
        if verbose:
            print(namesave+imifordir+': done')
            
        
            
#-------------------------------------------------------------
def Read_Batchelor(nameread, dirread, imrange, verbose=False, endread='.tsv', n0=3, normalised=True):        
    """
    Read Batchelor tensor TSV files generated by ``Batchelor_Batch``.

    :param nameread: Base file name
    :type nameread: str
    :param dirread: Directory containing TSV files
    :type dirread: str
    :param imrange: Image indices to read
    :type imrange: list or array
    :param verbose: Print progress information
    :type verbose: bool
    :param endread: File extension
    :type endread: str
    :param n0: Zero-padding for image indices
    :type n0: int
    :param normalised: If True, return volume-normalized tensors
    :type normalised: bool

    :return:
        - Labels
        - Coordinates
        - Volumes
        - Mesh surface areas
        - Batchelor tensors
    :rtype: tuple of lists
    """
    
    import numpy as np 
    import pandas as pd
    from FoamQuant.Helper import strindex
    import csv
    
    Llab=[]; LCoord=[]; Lvol=[]; Lmesharea=[]; LB=[]
    
    for imi in imrange: #image index loop
        lab=[]; Coord=[]; vol=[]; mesharea=[]; B=[]
        
        # image string index
        imifordir = strindex(imi, n0)

        props = pd.read_csv(dirread+nameread+imifordir+endread, sep='\t',engine="python",  quoting=csv.QUOTE_NONE)
        lab.append(np.asarray(props['lab']))
        z = np.asarray(props['z'])
        y = np.asarray(props['y'])
        x = np.asarray(props['x'])
        vol.append(np.asarray(props['vol']))
        mesharea.append(np.asarray(props['mesharea']))
        
        if normalised:
            b11 = np.asarray(props['b11'])
            b22 = np.asarray(props['b22'])
            b33 = np.asarray(props['b33'])
            b12 = np.asarray(props['b12'])
            b13 = np.asarray(props['b13'])
            b23 = np.asarray(props['b23'])
        else:
            b11 = np.asarray(props['B11'])
            b22 = np.asarray(props['B22'])
            b33 = np.asarray(props['B33'])
            b12 = np.asarray(props['B12'])
            b13 = np.asarray(props['B13'])
            b23 = np.asarray(props['B23'])

        for i in range(len(b11)): #bbl index loop
            B.append([[b11[i],b12[i],b13[i]],
                      [b12[i],b22[i],b23[i]],
                      [b13[i],b23[i],b33[i]]])
            Coord.append([z[i],y[i],x[i]])
        B = np.asarray(B)
        Coord = np.asarray(Coord)
        
        if verbose:
                print(nameread+imifordir+': done')
        
        Llab.append(lab)
        LCoord.append(Coord)
        Lvol.append(vol)
        Lmesharea.append(mesharea)
        LB.append(B)
    
    return Llab, LCoord, Lvol,Lmesharea, LB
