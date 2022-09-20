def Histogram(image, histtitle=False):
    """
    Plot a histogram from a given grey-level image

    Parameters
    ----------
        image : 3D numpy array

        histtitle : str, tittle of the figure (default = False)
    """
    
    from skimage.exposure import histogram
    fig, ax = plt.subplots(ncols=1, figsize=(5, 5))
    hist, hist_centers = skimage.exposure.histogram(image)
    ax.plot(hist_centers, hist, lw=2)
    ax.set_yscale('log')
    if histtitle != False:
        ax.set_title(histtitle)
