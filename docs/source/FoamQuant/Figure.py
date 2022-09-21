def Cut3D(image, zcut=False, ycut=False, xcut=False, showcuts=False, showaxes=False, showhistogram=False, histtitle=False, cmap='gray', interpolation=None, figblocksize=5, returnfig=False):
    """
    Plot a 3x1 figure showing three orthogonal cross-sections of the 3D image.
    
    :param image: 3D numpy array
    :type image: int
    :param zcut: Optional z cut value
    :type zcut: int or False
    :param ycut: Optional y cut value
    :type ycut: int or False
    :param xcut: Optional x cut value
    :type xcut: int or False 
    :param showcuts: Optional plot the orthogonal cuts
    :type showcuts: Bool
    :param showaxes: Optional plot the axes
    :type showaxes: Bool
    :param showhistogram: Optional add a grey-scale histogram, the figure become 3x1
    :type showhistogram: Bool
    :param histtitle: Optional add a tittle to the histogram
    :type histtitle: str or False
    :param cmap: Optional the color map used for the cuts, Default cmap = 'gray' 
    :type cmap: str or cmap type
    :param interpolation: Optional the type of interpolation, Default interpolation = None 
    :type interpolation: str or None
    :param figblocksize: Optional size of the subfigure, Default figblocksize = 5 
    :type figblocksize: float
    :param returnfig: Optional if should return the figure, if not returns None
    :type returnfig: Bool
    :return: None or fig
    """
  
    import numpy as np
    import matplotlib.pyplot as plt
    from skimage.exposure import histogram
 
    
    shapezyx = np.shape(image)
    if zcut == False:
        zcut = shapezyx[0]//2
    if ycut == False:
        ycut = shapezyx[1]//2
    if xcut == False:
        xcut = shapezyx[2]//2
    
    if showhistogram:
        from skimage.exposure import histogram
        fig, ax = plt.subplots(ncols=4, figsize=(4*figblocksize, figblocksize))
        hist, hist_centers = skimage.exposure.histogram(image)
        ax[3].plot(hist_centers, hist, lw=2)
        ax[3].set_yscale('log')
        if histtitle != False:
            ax[3].set_title(histtitle)
        plt.tight_layout()
    else:
        fig, ax = plt.subplots(ncols=3, figsize=(3*figblocksize, figblocksize))
    
    ax[0].imshow(image[zcut,:,:], cmap=cmap, interpolation=interpolation)
    ax[1].imshow(image[:,ycut,:], cmap=cmap, interpolation=interpolation)
    ax[2].imshow(image[:,:,xcut], cmap=cmap, interpolation=interpolation)
    
    if showcuts:
        ax[0].plot([1,shapezyx[2]-1,shapezyx[2]-1,1,1],[shapezyx[1]-1,shapezyx[1]-1,1,1,shapezyx[1]-1],'r',linewidth=3) #zcut
        ax[1].plot([1,shapezyx[2]-1,shapezyx[2]-1,1,1],[shapezyx[0]-1,shapezyx[0]-1,1,1,shapezyx[0]-1],'b',linewidth=3) #ycut
        ax[2].plot([1,shapezyx[1]-1,shapezyx[1]-1,1,1],[shapezyx[0]-1,shapezyx[0]-1,1,1,shapezyx[0]-1],'g',linewidth=3) #xcut
        
        ax[0].plot([1,shapezyx[2]-1],[ycut,ycut],'b',linewidth=3) #ycut
        ax[0].plot([xcut,xcut],[1,shapezyx[1]-1],'g',linewidth=3) #xcut
        
        ax[1].plot([1,shapezyx[2]-1],[zcut,zcut],'r',linewidth=3) #zcut
        ax[1].plot([xcut,xcut],[1,shapezyx[0]-1],'g',linewidth=3) #xcut
        
        ax[2].plot([1,shapezyx[1]-1],[zcut,zcut],'r',linewidth=3) #zcut
        ax[2].plot([ycut,ycut],[1,shapezyx[0]-1],'b',linewidth=3) #ycut
    
    if showaxes:
        ax[0].set_xlabel('x'); ax[0].set_ylabel('y')
        ax[1].set_xlabel('x'); ax[1].set_ylabel('z')
        ax[2].set_xlabel('y'); ax[2].set_ylabel('z')
        
    plt.tight_layout()
    
    if returnfig:
        return fig

      
      
      
def Histogram(image, histtitle=False):
    """
    Plot a 1x1 grey value histogram.
    
    :param image: 3D image.
    :type image: numpy array
    :return: None
    """

    from skimage.exposure import histogram
    fig, ax = plt.subplots(ncols=1, figsize=(5, 5))
    hist, hist_centers = skimage.exposure.histogram(image)
    ax.plot(hist_centers, hist, lw=2)
    ax.set_yscale('log')
    if histtitle != False:
        ax.set_title(histtitle)


# Refenece: https://github.com/delestro/rand_cmap
def RandomCmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):
    """
    Creates a random colormap for matplotlib. Reference: copied from https://github.com/delestro/rand_cmap
    
    :param nlabels: Number of labels (size of colormap)
    :type nlabels: int
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :type type: str
    :param first_color_black: Option to use first color as black, True or False
    :type first_color_black: Bool
    :param last_color_black: Option to use last color as black, True or False
    :type last_color_black: Bool
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :type verbose: Bool
    :return: matplotlib colormap
    """

    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np


    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap



def RandomCmap_json(Ncolors, namecmap, first_color_black=True):
    """
    Save a json random colormap to be used with ParaView or Tomviz.
    
    :param Ncolors: Number of labels (size of colormap)
    :type Ncolors: int
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :type type: str
    :param first_color_black: Option to use first color as black, True or False
    :type first_color_black: Bool
    :return: None
    """      
    
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
            x0 = i/Ncolors
        x1 = (i+1)/Ncolors
        
        randRGBcolors.append(x0)
        randRGBcolors.append(RGBcolors[0]); randRGBcolors.append(RGBcolors[1]); randRGBcolors.append(RGBcolors[2])
        randRGBcolors.append(x1)
        randRGBcolors.append(RGBcolors[0]); randRGBcolors.append(RGBcolors[1]); randRGBcolors.append(RGBcolors[2])

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

        
