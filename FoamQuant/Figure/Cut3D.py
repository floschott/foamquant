def Cut3D(image, zcut=False, ycut=False, xcut=False, histogram=False, showcuts=False, showaxes=False, showhistogram=False,
          histtitle=False, cmap='gray', interpolation=None, figblocksize=5, returnfig=False):
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