# Import dependencies
import logging
import numpy as np
import scipy.ndimage as spim
from scipy.special import erfc
from skimage.segmentation import relabel_sequential
from edt import edt
from skimage.morphology import ball, disk
try:
    from skimage.measure import marching_cubes
except ImportError:
    from skimage.measure import marching_cubes_lewiner as marching_cubes
import inspect
import time
from skimage import measure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# COPIED and MODIFIED FUNCTIONS FROM PoreSpy and SCIKIT-IMAGE

# Function 1: Region surface areas - COPIED
# Source: https://porespy.org/_modules/porespy/metrics/_meshtools.html#region_surface_areas
def region_surface_areas(regions, voxel_size=1, strel=None): 
    r"""
    Extract the surface area of each region in a labeled image.

    Optionally, it can also find the the interfacial area between all
    adjoining regions.

    Parameters
    ----------
    regions : ndarray
        An image of the pore space partitioned into individual pore regions.
        Note that zeros in the image will not be considered for area
        calculation.
    voxel_size : scalar
        The resolution of the image, expressed as the length of one side of a
        voxel, so the volume of a voxel would be **voxel_size**-cubed.  The
        default is 1.
    strel : array_like
        The structuring element used to blur the region.  If not provided,
        then a spherical element (or disk) with radius 1 is used.  See the
        docstring for ``mesh_region`` for more details, as this argument is
        passed to there.

    Returns
    -------
    areas : list
        A list containing the surface area of each region, offset by 1, such
        that the surface area of region 1 is stored in element 0 of the list.

    Examples
    --------
    `Click here
    <https://porespy.org/examples/metrics/reference/region_surface_areas.html>`_
    to view online example.

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
    r"""
    Creates round structuring element with the given radius and dimensionality

    Parameters
    ----------
    r : scalar
        The desired radius of the structuring element
    ndim : int
        The dimensionality of the element, either 2 or 3.
    smooth : boolean
        Indicates whether the faces of the sphere should have the little
        nibs (``True``) or not (``False``, default)

    Returns
    -------
    strel : ndarray
        A 3D numpy array of the structuring element

    Examples
    --------
    `Click here
    <https://porespy.org/examples/tools/reference/ps_round.html>`_
    to view online example.

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
    r"""
    Adjust slice indices to include additional voxles around the slice.

    This function does bounds checking to ensure the indices don't extend
    outside the image.

    Parameters
    ----------
    slices : list of slice objects
         A list (or tuple) of N slice objects, where N is the number of
         dimensions in the image.
    shape : array_like
        The shape of the image into which the slice objects apply.  This is
        used to check the bounds to prevent indexing beyond the image.
    pad : int or list of ints
        The number of voxels to expand in each direction.

    Returns
    -------
    slices : list of slice objects
        A list slice of objects with the start and stop attributes respectively
        incremented and decremented by 1, without extending beyond the image
        boundaries.

    Examples
    --------
    >>> from scipy.ndimage import label, find_objects
    >>> from porespy.tools import extend_slice
    >>> im = np.array([[1, 0, 0], [1, 0, 0], [0, 0, 1]])
    >>> labels = label(im)[0]
    >>> s = find_objects(labels)

    Using the slices returned by ``find_objects``, set the first label to 3

    >>> labels[s[0]] = 3
    >>> print(labels)
    [[3 0 0]
     [3 0 0]
     [0 0 2]]

    Next extend the slice, and use it to set the values to 4

    >>> s_ext = extend_slice(s[0], shape=im.shape, pad=1)
    >>> labels[s_ext] = 4
    >>> print(labels)
    [[4 4 0]
     [4 4 0]
     [4 4 2]]

    As can be seen by the location of the 4s, the slice was extended by 1, and
    also handled the extension beyond the boundary correctly.

    Examples
    --------
    `Click here
    <https://porespy.org/examples/tools/reference/extend_slice.html>`_
    to view online example.

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
    r"""
    Creates a tri-mesh of the provided region using the marching cubes
    algorithm

    Parameters
    ----------
    im : ndarray
        A boolean image with ``True`` values indicating the region of interest
    strel : ndarray
        The structuring element to use when blurring the region.  The blur is
        perfomed using a simple convolution filter.  The point is to create a
        greyscale region to allow the marching cubes algorithm some freedom
        to conform the mesh to the surface.  As the size of ``strel`` increases
        the region will become increasingly blurred and inaccurate. The default
        is a spherical element with a radius of 1.

    Returns
    -------
    mesh : tuple
        A named-tuple containing ``faces``, ``verts``, ``norm``, and ``val``
        as returned by ``scikit-image.measure.marching_cubes`` function.

    Examples
    --------
    `Click here
    <https://porespy.org/examples/tools/reference/mesh_region.html>`_
    to view online example.

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
    r"""
    Checks for whether the input image contains singleton axes and logs
    a proper warning in case found.

    Parameters
    ----------
    im : ndarray
        Input image.

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
    r"""
    Calculate the surface area of a meshed region

    Parameters
    ----------
    mesh : tuple
        The tuple returned from the ``mesh_region`` function
    verts : array
        An N-by-ND array containing the coordinates of each mesh vertex
    faces : array
        An N-by-ND array indicating which elements in ``verts`` form a mesh
        element.

    Returns
    -------
    surface_area : float
        The surface area of the mesh, calculated by
        ``skimage.measure.mesh_surface_area``

    Notes
    -----
    This function simply calls ``scikit-image.measure.mesh_surface_area``, but
    it allows for the passing of the ``mesh`` tuple returned by the
    ``mesh_region`` function, entirely for convenience.

    Examples
    --------
    `Click here
    <https://porespy.org/examples/metrics/reference/mesh_surface_area.html>`_
    to view online example.

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
    r"""
    Visualizes the mesh of a region as obtained by ``get_mesh`` function in
    the ``metrics`` submodule.

    Parameters
    ----------
    mesh : tuple
        A mesh returned by ``skimage.measure.marching_cubes``

    Returns
    -------
    fig : Matplotlib figure
        A handle to a matplotlib 3D axis

    Examples
    --------
    `Click here
    <https://porespy.org/examples/visualization/reference/show_mesh.html>`_
    to view online example.
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
    r"""
    Extract the Batchelor tensor of each region in a labeled image.

    Parameters
    ----------
    regions : ndarray
        An image of the pore space partitioned into individual pore regions.
        Note that zeros in the image will not be considered for area
        calculation.
    voxel_size : scalar
        The resolution of the image, expressed as the length of one side of a
        voxel, so the volume of a voxel would be **voxel_size**-cubed.  The
        default is 1.
    strel : array_like (NOT TESTED YET FOR BATCHELOR TENSOR!)
        The structuring element used to blur the region.  If not provided,
        then a spherical element (or disk) with radius 1 is used.  See the
        docstring for ``mesh_region`` for more details, as this argument is
        passed to there.

    Returns
    -------
    areas : list
        A list containing the Batchelor tensor of each region, offset by 1, such
        that the Batchelor tensor of region 1 is stored in element 0 of the list.

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
    r"""
    Calculate the Batchelor tensor of a meshed region

    Parameters
    ----------
    mesh : tuple
        The tuple returned from the ``mesh_region`` function
    verts : array
        An N-by-ND array containing the coordinates of each mesh vertex
    faces : array
        An N-by-ND array indicating which elements in ``verts`` form a mesh
        element.

    Returns
    -------
    mesh_for_batchelor : float 3x3 tensor
        The Batchelor tensor of the mesh, calculated by the modified
        ``skimage.measure.mesh_surface_area`` function called mesh_batchelor

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
    """Compute Batchelor tensor, given vertices and triangular faces.

    Parameters
    ----------
    verts : (V, 3) array of floats
        Array containing (x, y, z) coordinates for V unique mesh vertices.
    faces : (F, 3) array of ints
        List of length-3 lists of integers, referencing vertex coordinates as
        provided in `verts`.

    Returns
    -------
    B : float 3x3 tensor
        Batchelor tensor of mesh. Units now [coordinate units] ** 2.
    area : float
        Surface area of mesh. Units now [coordinate units] ** 2.

    Notes
    -----
    The arguments expected by this function are the first two outputs from
    `skimage.measure.marching_cubes`. For unit correct output, ensure correct
    `spacing` was passed to `skimage.measure.marching_cubes`.

    This algorithm works properly only if the ``faces`` provided are all
    triangles.

    See Also
    --------
    skimage.measure.marching_cubes

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
    Run region_Batchelor_tensors function on a batch of images and save the outputs as .tsv
    
    :param readdir: Labeled images folder
    :type readdir: str
    :param readdir: folder to save the .tsv doc
    :type readdir: str
    :param readend: tiff image saving end, default is '.tiff'
    :type readend: str
    :param imrange: list of image indexes
    :type imrange: int array
    :param verbose: if True, print verbose including the number of labels
    :type verbose: Bool
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
    Read the saved files generated by the .tsv Batchelor_Batch function on a batch of images
    
    :param nameread: file name (without the ending number)
    :type nameread: str
    :param dirread: folder directory
    :type dirread: str
    :param imrange: list of image indexes
    :type imrange: int array
    :param verbose: if True, print the image number
    :type verbose: Bool
    :param readend: ending
    :type readend: str
    :param n0: number of 0 in the indexing
    :type n0: int
    :return: label, Coordinates, volume, mesh_area, Batchelor_tensor
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