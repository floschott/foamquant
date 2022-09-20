def Cut3D(image, histtitle=False):
    """
    This function plot an grey-level histogram of the image.

    Parameters
    -----------
        image : 3D numpy array
            This is the grey-value image

        histtitle : str (optional, default = None)
            Add a tittle to the hitogram

    Returns
    --------
        None
    """
    from skimage.exposure import histogram
    fig, ax = plt.subplots(ncols=1, figsize=(5, 5))
    hist, hist_centers = skimage.exposure.histogram(image)
    ax.plot(hist_centers, hist, lw=2)
    ax.set_yscale('log')
    if histtitle != False:
        ax.set_title(histtitle)
