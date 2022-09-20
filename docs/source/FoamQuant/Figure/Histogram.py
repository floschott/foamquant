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
