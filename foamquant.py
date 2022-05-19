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



def json_rand_dictionary(Ncolors, namecmap, first_color_black=True):
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np
    import json
    
    randHSVcolors = [(np.random.uniform(low=0.0, high=1), np.random.uniform(low=0.2, high=1), np.random.uniform(low=0.9, high=1)) for i in range(Ncolors)]
    
    randRGBcolors = []
    for i in range(len(randHSVcolors)):
        RGBcolors = colorsys.hsv_to_rgb(randHSVcolors[i][0], randHSVcolors[i][1], randHSVcolors[i][2])
        
        x0 = i/Ncolors
        if i >0:
            x0 = i/Ncolors#+0.001
        x1 = (i+1)/Ncolors
        
        randRGBcolors.append(x0)
        randRGBcolors.append(RGBcolors[0]); randRGBcolors.append(RGBcolors[1]); randRGBcolors.append(RGBcolors[2])
        randRGBcolors.append(x1)
        randRGBcolors.append(RGBcolors[0]); randRGBcolors.append(RGBcolors[1]); randRGBcolors.append(RGBcolors[2])
    
    #randRGBcolors = np.asarray(randRGBcolors)
    
    
    #print(randRGBcolors)
    if first_color_black == True:
        randRGBcolors[1:4] = [0,0,0]
        randRGBcolors[5:8] = [0,0,0]
    
    
    json_cmap = [
	{
		"ColorSpace" : "HSV",
		"Name" : namecmap,
		"RGBPoints" : randRGBcolors
    }
    ]
    
    with open(namecmap+".json", "w") as outfile:
        json.dump(json_cmap, outfile)
