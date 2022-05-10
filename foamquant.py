"""
Foamquant - Python library for processing time series of large 3D images of evolving foams.
"""

__version__ = "0.1.0"

def RandomCmap(ncolors = 10000, namecmap = 'random_cmap', first_color_black=True):
    
    ''' 
    Return a random colormap
    Greatly inspired by: https://github.com/delestro/rand_cmap
    
    :param ncolors: The number of colors.
    :type ncolors: int
    :param namecmap: The random colormap name.
    :type namecmap: str
    :return: The random colormap random_cmap.
    '''
    
    #Libraries
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    
    #Random HSV list
    randHSVcolors = [(np.random.uniform(low=0.0, high=1), 
                      np.random.uniform(low=0.2, high=1), 
                      np.random.uniform(low=0.9, high=1)) for i in range(ncolors)]
    
    #HSV to RGB
    randRGBcolors = []
    for HSVcolor in randHSVcolors:
        randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))
    #First color black
    if first_color_black:
        randRGBcolors[0] = [0, 0, 0]
    #LinearSegmentedColormap
    random_cmap = LinearSegmentedColormap.from_list(name, randRGBcolors, N=ncolors)
    
    return(random_cmap)
