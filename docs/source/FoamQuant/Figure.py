# Image ----------------
def Cut3D(image, zcut=False, ycut=False, xcut=False, showcuts=False, nameaxes=None, cmap='gray', interpolation=None, figblocksize=5, returnfig=False, vmin=None, vmax=None, printminmax=False, colorbars=False, constrained_layout=True, thick=0, colorbarlab=None):
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
    :param cmap: Optional the color map used for the cuts, Default cmap = 'gray' 
    :type cmap: str or cmap type
    :param interpolation: Optional type of interpolation, Default interpolation = None 
    :type interpolation: str or None
    :param figblocksize: Optional size of the subfigure, Default figblocksize = 5 
    :type figblocksize: float
    :param returnfig: Optional, if should return the figure, if not returns None
    :type returnfig: Bool
    :param vmin: Optional, min value for the color range
    :type vmin: Bool
    :param vmax: Optional, max value for the color range
    :type vmax: Bool
    :param printminmax: Optional, print min and max for the whole image and the three projections
    :type printminmax: Bool
    :return: None or fig type
    """    
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    shapezyx = np.shape(image)
    if zcut == False:
        zcut = shapezyx[0]//2
    if ycut == False:
        ycut = shapezyx[1]//2
    if xcut == False:
        xcut = shapezyx[2]//2
    
    fig, ax = plt.subplots(ncols=3, figsize=(3*figblocksize, figblocksize), constrained_layout=constrained_layout)
        
    if thick == 0:
        if vmin!=None and vmax!=None:
            print('vmin =',vmin, 'vmax =',vmax) 
            neg1 = ax[0].imshow(image[zcut,:,:], cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
            neg2 = ax[1].imshow(image[:,ycut,:], cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
            neg3 = ax[2].imshow(image[:,:,xcut], cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
        else:    
            neg1 = ax[0].imshow(image[zcut,:,:], cmap=cmap, interpolation=interpolation)
            neg2 = ax[1].imshow(image[:,ycut,:], cmap=cmap, interpolation=interpolation)
            neg3 = ax[2].imshow(image[:,:,xcut], cmap=cmap, interpolation=interpolation)
    
    if thick > 0:
        if vmin!=None and vmax!=None:
            print('vmin =',vmin, 'vmax =',vmax) 
            neg1 = ax[0].imshow(np.nanmean(image[zcut-thick:zcut+thick,:,:],0), cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
            neg2 = ax[1].imshow(np.nanmean(image[:,ycut-thick:ycut+thick,:],1), cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
            neg3 = ax[2].imshow(np.nanmean(image[:,:,xcut-thick:xcut+thick],2), cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
        else:    
            neg1 = ax[0].imshow(np.nanmean(image[zcut-thick:zcut+thick,:,:],0), cmap=cmap, interpolation=interpolation)
            neg2 = ax[1].imshow(np.nanmean(image[:,ycut-thick:ycut+thick,:],1), cmap=cmap, interpolation=interpolation)
            neg3 = ax[2].imshow(np.nanmean(image[:,:,xcut-thick:xcut+thick],2), cmap=cmap, interpolation=interpolation)
    
    if showcuts:
        ax[0].plot([1,shapezyx[2]-1,shapezyx[2]-1,0,0],[shapezyx[1]-1,shapezyx[1]-1,0,0,shapezyx[1]-1],'r',linewidth=3) #zcut
        ax[1].plot([1,shapezyx[2]-1,shapezyx[2]-1,0,0],[shapezyx[0]-1,shapezyx[0]-1,0,0,shapezyx[0]-1],'b',linewidth=3) #ycut
        ax[2].plot([1,shapezyx[1]-1,shapezyx[1]-1,0,0],[shapezyx[0]-1,shapezyx[0]-1,0,0,shapezyx[0]-1],'g',linewidth=3) #xcut
        
        ax[0].plot([1,shapezyx[2]-1],[ycut,ycut],'b',linewidth=3) #ycut
        ax[0].plot([xcut,xcut],[1,shapezyx[1]-1],'g',linewidth=3) #xcut
        
        ax[1].plot([1,shapezyx[2]-1],[zcut,zcut],'r',linewidth=3) #zcut
        ax[1].plot([xcut,xcut],[1,shapezyx[0]-1],'g',linewidth=3) #xcut
        
        ax[2].plot([1,shapezyx[1]-1],[zcut,zcut],'r',linewidth=3) #zcut
        ax[2].plot([ycut,ycut],[1,shapezyx[0]-1],'b',linewidth=3) #ycut
    
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])
        
        
    plt.tight_layout()
    
    if printminmax:
        print('MIN:',np.nanmin(image),'MAX:',np.nanmax(image))
        print('Min0:',np.nanmin(image[zcut,:,:]),'Min0:',np.nanmax(image[zcut,:,:]))
        print('Min1:',np.nanmin(image[:,ycut,:]),'Max1',np.nanmax(image[:,ycut,:]))
        print('Min2:',np.nanmin(image[:,:,xcut]),'Max2:',np.nanmax(image[:,:,xcut]))
        
    if colorbars:
        fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.4, hspace=0.1)
        cb_ax = fig.add_axes([0.83, 0.1, 0.05, 0.8])
        cbar = fig.colorbar(neg1, cax=cb_ax)
        
    if colorbarlab != None:
        cbar.set_label(colorbarlab)
    
    if returnfig:
        return fig,ax,[neg1,neg2,neg3]
    
    
    
    
    
def Proj3D(image, nameaxes=None, cmap='gray', interpolation=None, figblocksize=5, returnfig=False, vmin=None, vmax=None, printminmax=False, colorbars=False, constrained_layout=True):
    """
    Plot a 3x1 figure showing three orthogonal projections of the 3D image.
    
    :param image: 3D image
    :type image: numpy array
    :param showaxes: Optional plot the axes
    :type showaxes: Bool
    :param cmap: Optional the color map used for the projections, Default cmap = 'gray' 
    :type cmap: str or cmap type
    :param interpolation: Optional the type of interpolation, Default interpolation = None 
    :type interpolation: str or None
    :param figblocksize: Optional size of the subfigure, Default figblocksize = 5 
    :type figblocksize: float
    :param returnfig: Optional if should return the figure, if not returns None
    :type returnfig: Bool
    :param vmin: Optional, min value for the color range
    :type vmin: Bool
    :param vmax: Optional, max value for the color range
    :type vmax: Bool
    :param printminmax: Optional, print min and max for the whole image and the three projections
    :type printminmax: Bool
    :return: None or fig
    """
    
    import numpy as np
    import matplotlib.pyplot as plt    
    
   
    fig, ax = plt.subplots(ncols=3, figsize=(3*figblocksize, figblocksize), constrained_layout=constrained_layout)
        
    if vmin!=None and vmax!=None:
        print('vmin =',vmin, 'vmax =',vmax) 
        neg1 = ax[0].imshow(np.nanmean(image,0), cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
        neg2 = ax[1].imshow(np.nanmean(image,1), cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
        neg3 = ax[2].imshow(np.nanmean(image,2), cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax)
    else:    
        neg1 = ax[0].imshow(np.nanmean(image,0), cmap=cmap, interpolation=interpolation)
        neg2 = ax[1].imshow(np.nanmean(image,1), cmap=cmap, interpolation=interpolation)
        neg3 = ax[2].imshow(np.nanmean(image,2), cmap=cmap, interpolation=interpolation)
    
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])
        
    plt.tight_layout()
    
    if printminmax:
        print('MIN:',np.nanmin(image),'MAX:',np.nanmax(image))
        print('Min0:',np.nanmin(np.nanmean(image,0)),'Min0:',np.nanmax(np.nanmean(image,0)))
        print('Min1:',np.nanmin(np.nanmean(image,1)),'Max1',np.nanmax(np.nanmean(image,1)))
        print('Min2:',np.nanmin(np.nanmean(image,2)),'Max2:',np.nanmax(np.nanmean(image,2)))
    
    if colorbars:
        fig.colorbar(neg1)
        fig.colorbar(neg2)
        fig.colorbar(neg3)
    
    if returnfig:
        return fig,ax,[neg1,neg2,neg3]

    
# Tensor----------------
def CutTensor3D(Grids, Coord, US, Count, zcut=False, ycut=False, xcut=False, scale_factor = 1,figblocksize=7, vmin=-1, vmax=1, Countmin=5, showscale=False, nameaxes=None, DeviatoricType=2,colorbarlab=None):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    shapezyx = np.shape(Coord)[:3]
    if zcut == False:
        zcut = shapezyx[0]//2
    if ycut == False:
        ycut = shapezyx[1]//2
    if xcut == False:
        xcut = shapezyx[2]//2
    
    fig, ax = plt.subplots(ncols=3, figsize=(3*figblocksize, figblocksize), constrained_layout=True)
    
    LLustxcut=np.zeros((shapezyx[0], shapezyx[1]))
    LLustzcut=np.zeros((shapezyx[1], shapezyx[2]))
    LLustycut=np.zeros((shapezyx[0], shapezyx[2]))
    LLxgridxcut=np.zeros((shapezyx[0], shapezyx[1]))
    LLygridxcut=np.zeros((shapezyx[0], shapezyx[1]))
    LLxgridycut=np.zeros((shapezyx[0], shapezyx[2]))
    LLygridycut=np.zeros((shapezyx[0], shapezyx[2]))
    LLxgridzcut=np.zeros((shapezyx[1], shapezyx[2]))
    LLygridzcut=np.zeros((shapezyx[1], shapezyx[2]))
    
    # xcut (only in z,y (0,1))
    for zi in range(shapezyx[0]):
        for yi in range(shapezyx[1]):
            
            if Count[zi][yi][xcut] > Countmin:
            
                z,y,x = np.diag(Coord[zi][yi][xcut])
                x,y = y,z
                #us = US[zi][yi][xcut]
                if DeviatoricType==3: # remove 1/2*log((a*b*c)**-3)
                    us = UMfromM(US[zi][yi][xcut])
                elif DeviatoricType==2: # remove log((a*b*c)**-3)
                    us = USfromS(US[zi][yi][xcut])
                elif DeviatoricType==1: # remove trace: a+b+c
                    us = SigdevfromSig(US[zi][yi][xcut])
                elif DeviatoricType==0: # keep as it is
                    us = US[zi][yi][xcut]
                else:
                    print('Error: provide DeviatoricType')
                    break
                us22 = us[0:2,0:2]
                ust = us[2][2]
                LLustxcut[zi][yi]=ust
                LLxgridxcut[zi][yi]=x
                LLygridxcut[zi][yi]=y

                Val, Vect = np.linalg.eig(us22)

                #ordering: min/max eigval
                mMV = np.zeros((2,2))
                mMv = np.sort(Val)
                for orderindex in range(2):
                    for vindex in range(2):
                        if mMv[orderindex] == Val[vindex]:
                            mMV[:,orderindex] = Vect[:,vindex]

                eigvm = mMv[0]
                eigvM = mMv[1]
                eigVMx = mMV[1,1]
                eigVmx = mMV[1,0]
                eigVMy = mMV[0,1]
                eigVmy = mMV[0,0]

                isnan = True
                if eigVMx > 0 and eigVMy >= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)); isnan = False
                elif eigVMx <= 0 and eigVMy > 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + np.pi/2; isnan = False
                elif eigVMx < 0 and eigVMy <= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)) + np.pi; isnan = False
                elif eigVMx >= 0 and eigVMy < 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + 3*np.pi/2; isnan = False
                elif eigVMx == 0:
                    beta = np.pi/2; isnan = False
                elif eigVMy == 0:
                    beta = np.pi/2; isnan = False

                if isnan == False:
                    alphas = np.linspace(0,2,30)*np.pi
                    a = scale_factor*abs(eigvM)
                    b = scale_factor*abs(eigvm)
                    xs = x + np.cos(beta)*np.cos(alphas)*a - np.sin(beta)*np.sin(alphas)*b
                    ys = y + np.sin(beta)*np.cos(alphas)*a + np.cos(beta)*np.sin(alphas)*b
                    LinePlusx = [x + np.cos(beta)*a, x - np.cos(beta)*a]
                    LinePlusy = [y + np.sin(beta)*a, y - np.sin(beta)*a]
                    LineMinusx = [x + np.cos(beta+np.pi/2)*b, x - np.cos(beta+np.pi/2)*b]
                    LineMinusy = [y + np.sin(beta+np.pi/2)*b, y - np.sin(beta+np.pi/2)*b]
                    #if np.abs(eigvM) > np.abs(eigvm):
                    #test for ploting the ellipse
                    if ust<0: # is compressed along the third direction
                        ax[2].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    if ust>0: # is elongated along the third direction
                        ax[2].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    #test for ploting the bars
                    if eigvM>0: #elongated along the major axis /else is compressed
                        ax[2].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if eigvm>0: #elongated along the second major axis /else is compressed
                        ax[2].plot(LineMinusx,LineMinusy, 'k-', alpha = 1, linewidth = 3)
            else:
                LLustxcut[zi][yi]=np.nan
                LLxgridxcut[zi][yi]=np.nan
                LLygridxcut[zi][yi]=np.nan
                
    # zcut (only in y,x (1,2))
    for yi in range(shapezyx[1]):
        for xi in range(shapezyx[2]):
            if Count[zcut][yi][xi] > Countmin:
                z,y,x = np.diag(Coord[zcut][yi][xi])
                x,y = x,y
                #us = US[zi][yi][xcut]
                if DeviatoricType==3: # remove 1/2*log((a*b*c)**-3)
                    us = UMfromM(US[zcut][yi][xi])
                elif DeviatoricType==2: # remove log((a*b*c)**-3)
                    us = USfromS(US[zcut][yi][xi])
                elif DeviatoricType==1: # remove trace: a+b+c
                    us = SigdevfromSig(US[zcut][yi][xi])
                elif DeviatoricType==0: # keep as it is
                    us = US[zcut][yi][xi]
                else:
                    print('Error: provide DeviatoricType')
                    break
                us22 = us[1:,1:]
                ust = us[0][0]
                LLustzcut[yi][xi]=ust
                LLxgridzcut[yi][xi]=x
                LLygridzcut[yi][xi]=y
                
                Val, Vect = np.linalg.eig(us22)
                
                #ordering: min/max eigval
                mMV = np.zeros((2,2))
                mMv = np.sort(Val)
                for orderindex in range(2):
                    for vindex in range(2):
                        if mMv[orderindex] == Val[vindex]:
                            mMV[:,orderindex] = Vect[:,vindex]

                eigvm = mMv[0]
                eigvM = mMv[1]
                eigVMx = mMV[1,1]
                eigVmx = mMV[1,0]
                eigVMy = mMV[0,1]
                eigVmy = mMV[0,0]

                isnan = True
                if eigVMx > 0 and eigVMy >= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)); isnan = False
                elif eigVMx <= 0 and eigVMy > 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + np.pi/2; isnan = False
                elif eigVMx < 0 and eigVMy <= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)) + np.pi; isnan = False
                elif eigVMx >= 0 and eigVMy < 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + 3*np.pi/2; isnan = False
                elif eigVMx == 0:
                    beta = np.pi/2; isnan = False
                elif eigVMy == 0:
                    beta = np.pi/2; isnan = False

                if isnan == False:
                    alphas = np.linspace(0,2,30)*np.pi
                    a = scale_factor*abs(eigvM)
                    b = scale_factor*abs(eigvm)
                    xs = x + np.cos(beta)*np.cos(alphas)*a - np.sin(beta)*np.sin(alphas)*b
                    ys = y + np.sin(beta)*np.cos(alphas)*a + np.cos(beta)*np.sin(alphas)*b
                    LinePlusx = [x + np.cos(beta)*a, x - np.cos(beta)*a]
                    LinePlusy = [y + np.sin(beta)*a, y - np.sin(beta)*a]
                    LineMinusx = [x + np.cos(beta+np.pi/2)*b, x - np.cos(beta+np.pi/2)*b]
                    LineMinusy = [y + np.sin(beta+np.pi/2)*b, y - np.sin(beta+np.pi/2)*b]
                    #if np.abs(eigvM) > np.abs(eigvm):                 
                    if ust<0:
                        ax[0].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    if ust>0:
                        ax[0].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    if eigvM>0:
                        ax[0].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if eigvm>0:
                        ax[0].plot(LineMinusx,LineMinusy, 'k-', alpha = 1, linewidth = 3)
            else:
                LLustzcut[yi][xi]=np.nan
                LLxgridzcut[yi][xi]=np.nan
                LLygridzcut[yi][xi]=np.nan

    # ycut (only in z,x (0,2))
    for zi in range(shapezyx[0]):
        for xi in range(shapezyx[2]):
            if Count[zi][ycut][xi] > Countmin:
                z,y,x = np.diag(Coord[zi][ycut][xi])
                x,y = x,z
                #us = US[zi][ycut][xi]
                if DeviatoricType==3: # remove 1/2*log((a*b*c)**-3)
                    us = UMfromM(US[zi][ycut][xi])
                if DeviatoricType==2: # remove log((a*b*c)**-3)
                    us = USfromS(US[zi][ycut][xi])
                elif DeviatoricType==1: # remove trace: a+b+c
                    us = SigdevfromSig(US[zi][ycut][xi])
                elif DeviatoricType==0: # keep as it is
                    us = US[zi][ycut][xi]
                else:
                    print('Error: provide DeviatoricType')
                    break
                us22 = np.zeros((2,2))
                us22[0,0] = us[0,0]
                us22[1,0] = us[2,0]
                us22[0,1] = us[0,2]
                us22[1,1] = us[2,2]
                ust = us[1,1]
                LLustycut[zi][xi]=ust
                
                LLxgridycut[zi][xi]=x
                LLygridycut[zi][xi]=y

                Val, Vect = np.linalg.eig(us22)

                #ordering: min/max eigval
                mMV = np.zeros((2,2))
                mMv = np.sort(Val)
                for orderindex in range(2):
                    for vindex in range(2):
                        if mMv[orderindex] == Val[vindex]:
                            mMV[:,orderindex] = Vect[:,vindex]

                eigvm = mMv[0]
                eigvM = mMv[1]
                eigVMx = mMV[1,1]
                eigVmx = mMV[1,0]
                eigVMy = mMV[0,1]
                eigVmy = mMV[0,0]

                isnan = True
                if eigVMx > 0 and eigVMy >= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)); isnan = False
                elif eigVMx <= 0 and eigVMy > 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + np.pi/2; isnan = False
                elif eigVMx < 0 and eigVMy <= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)) + np.pi; isnan = False
                elif eigVMx >= 0 and eigVMy < 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + 3*np.pi/2; isnan = False
                elif eigVMx == 0:
                    beta = np.pi/2; isnan = False
                elif eigVMy == 0:
                    beta = np.pi/2; isnan = False

                if isnan == False:
                    alphas = np.linspace(0,2,30)*np.pi
                    a = scale_factor*abs(eigvM)
                    b = scale_factor*abs(eigvm)
                    xs = x + np.cos(beta)*np.cos(alphas)*a - np.sin(beta)*np.sin(alphas)*b
                    ys = y + np.sin(beta)*np.cos(alphas)*a + np.cos(beta)*np.sin(alphas)*b
                    LinePlusx = [x + np.cos(beta)*a, x - np.cos(beta)*a]
                    LinePlusy = [y + np.sin(beta)*a, y - np.sin(beta)*a]
                    LineMinusx = [x + np.cos(beta+np.pi/2)*b, x - np.cos(beta+np.pi/2)*b]
                    LineMinusy = [y + np.sin(beta+np.pi/2)*b, y - np.sin(beta+np.pi/2)*b]
                    #if np.abs(eigvM) > np.abs(eigvm):
                    if ust<0:
                        ax[1].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    if ust>0:
                        ax[1].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    if eigvM>0:
                        ax[1].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if eigvm>0:
                        ax[1].plot(LineMinusx,LineMinusy, 'k-', alpha = 1, linewidth = 3)   
            else:
                LLustycut[zi][xi]=np.nan
                LLxgridycut[zi][xi]=np.nan
                LLygridycut[zi][xi]=np.nan
    ax[0].axis('equal')
    ax[1].axis('equal')
    ax[2].axis('equal')
    if len(np.shape(showscale))>0:
        for ra in showscale:
            ra= ra*scale_factor
            xs = np.cos(alphas)*ra
            ys = ra + np.sin(alphas)*ra
            ax[0].plot(xs,ys,'k')
            ax[1].plot(xs,ys,'k')
            ax[2].plot(xs,ys,'k')
    print('Normal min/max ax0',np.nanmin(LLustzcut), np.nanmax(LLustzcut))
    print('Normal min/max ax1',np.nanmin(LLustycut), np.nanmax(LLustycut))
    print('Normal min/max ax2',np.nanmin(LLustxcut), np.nanmax(LLustxcut))
    c=ax[2].pcolor(Grids[1], Grids[0], LLustxcut, vmin = vmin, vmax = vmax,cmap = 'seismic', shading='nearest', alpha=0.3)
    c=ax[1].pcolor(Grids[2], Grids[0], LLustycut, vmin = vmin, vmax = vmax,cmap = 'seismic', shading='nearest', alpha=0.3)
    c=ax[0].pcolor(Grids[2], Grids[1], LLustzcut, vmin = vmin, vmax = vmax,cmap = 'seismic', shading='nearest', alpha=0.3)
    cbar = fig.colorbar(c, ax=ax)
    
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])

    if colorbarlab != None:
        cbar.set_label(colorbarlab)
        
        
        
        
def ProjTensor3D(Grids, Coord, US, Count, scale_factor = 1,figblocksize=7, vmin=-1, vmax=1, Countmin=5, showscale=False, nameaxes=None, returnfig=None, dirsave=None, dpi=200, formt='jpeg', alpha=0.3, DeviatoricType=2,colorbarlab=None):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    shapezyx = np.shape(Coord)[:3]
    
    if returnfig==None:
        fig, ax = plt.subplots(ncols=3, figsize=(3*figblocksize, figblocksize), constrained_layout=True)
    else:
        fig = returnfig[0]
        ax = returnfig[1]
    
    LLustxcut=np.zeros((shapezyx[0], shapezyx[1]))
    LLustzcut=np.zeros((shapezyx[1], shapezyx[2]))
    LLustycut=np.zeros((shapezyx[0], shapezyx[2]))
    LLxgridxcut=np.zeros((shapezyx[0], shapezyx[1]))
    LLygridxcut=np.zeros((shapezyx[0], shapezyx[1]))
    LLxgridycut=np.zeros((shapezyx[0], shapezyx[2]))
    LLygridycut=np.zeros((shapezyx[0], shapezyx[2]))
    LLxgridzcut=np.zeros((shapezyx[1], shapezyx[2]))
    LLygridzcut=np.zeros((shapezyx[1], shapezyx[2]))
    
    # xcut (only in z,y (0,1))
    USzy = np.nanmean(US, 2)
    Countzy = np.nanmean(Count, 2)
    Coordzy = np.nanmean(Coord, 2)
    USyx = np.nanmean(US, 0)
    Countyx = np.nanmean(Count, 0)
    Coordyx = np.nanmean(Coord, 0)
    USzx = np.nanmean(US, 1)
    Countzx = np.nanmean(Count, 1)
    Coordzx = np.nanmean(Coord, 1)
    
    for zi in range(shapezyx[0]):
        for yi in range(shapezyx[1]):
            if Countzy[zi][yi] > Countmin:
                z,y,x = np.diag(Coordzy[zi][yi])
                x,y = y,z
                if DeviatoricType==3: # remove 1/2*log((a*b*c)**-3)
                    us = UMfromM(USzy[zi][yi])
                elif DeviatoricType==2: # remove log((a*b*c)**-3)
                    us = USfromS(USzy[zi][yi])
                elif DeviatoricType==1: # remove trace: a+b+c
                    us = SigdevfromSig(USzy[zi][yi])
                elif DeviatoricType==0: # keep as it is
                    us = USzy[zi][yi]
                else:
                    print('Error: provide DeviatoricType')
                    break
                us22 = us[0:2,0:2]
                ust = us[2][2]
                LLustxcut[zi][yi]=ust
                LLxgridxcut[zi][yi]=x
                LLygridxcut[zi][yi]=y
                
                Val, Vect = np.linalg.eig(us22)

                #ordering: min/max eigval
                mMV = np.zeros((2,2))
                mMv = np.sort(Val)
                for orderindex in range(2):
                    for vindex in range(2):
                        if mMv[orderindex] == Val[vindex]:
                            mMV[:,orderindex] = Vect[:,vindex]

                eigvm = mMv[0]
                eigvM = mMv[1]
                eigVMx = mMV[1,1]
                eigVmx = mMV[1,0]
                eigVMy = mMV[0,1]
                eigVmy = mMV[0,0]

                isnan = True
                if eigVMx > 0 and eigVMy >= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)); isnan = False
                elif eigVMx <= 0 and eigVMy > 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + np.pi/2; isnan = False
                elif eigVMx < 0 and eigVMy <= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)) + np.pi; isnan = False
                elif eigVMx >= 0 and eigVMy < 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + 3*np.pi/2; isnan = False
                elif eigVMx == 0:
                    beta = np.pi/2; isnan = False
                elif eigVMy == 0:
                    beta = np.pi/2; isnan = False

                if isnan == False:
                    alphas = np.linspace(0,2,30)*np.pi
                    a = scale_factor*abs(eigvM)
                    b = scale_factor*abs(eigvm)
                    xs = x + np.cos(beta)*np.cos(alphas)*a - np.sin(beta)*np.sin(alphas)*b
                    ys = y + np.sin(beta)*np.cos(alphas)*a + np.cos(beta)*np.sin(alphas)*b
                    LinePlusx = [x + np.cos(beta)*a, x - np.cos(beta)*a]
                    LinePlusy = [y + np.sin(beta)*a, y - np.sin(beta)*a]
                    LineMinusx = [x + np.cos(beta+np.pi/2)*b, x - np.cos(beta+np.pi/2)*b]
                    LineMinusy = [y + np.sin(beta+np.pi/2)*b, y - np.sin(beta+np.pi/2)*b]
                    #if np.abs(eigvM) > np.abs(eigvm):
                    #    ax[2].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    #else:
                    #    ax[2].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    #ax[2].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if ust<0:
                        ax[2].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    if ust>0:
                        ax[2].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    if eigvM>0:
                        ax[2].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if eigvm>0:
                        ax[2].plot(LineMinusx,LineMinusy, 'k-', alpha = 1, linewidth = 3)         
                else:
                    print('isnan')
            else:
                LLustxcut[zi][yi]=np.nan
                LLxgridxcut[zi][yi]=np.nan
                LLygridxcut[zi][yi]=np.nan
                
    # zcut (only in y,x (1,2))
    for yi in range(shapezyx[1]):
        for xi in range(shapezyx[2]):
            if Countyx[yi][xi] > Countmin:
                z,y,x = np.diag(Coordyx[yi][xi])
                x,y = x,y
                if DeviatoricType==3: # remove 1/2*log((a*b*c)**-3)
                    us = UMfromM(USyx[yi][xi])
                elif DeviatoricType==2: # remove log((a*b*c)**-3)
                    us = USfromS(USyx[yi][xi])
                elif DeviatoricType==1: # remove trace: a+b+c
                    us = SigdevfromSig(USyx[yi][xi])
                elif DeviatoricType==0: # keep as it is
                    us = USyx[yi][xi]
                else:
                    print('Error: provide DeviatoricType')
                    break
                us22 = us[1:,1:]
                ust = us[0][0]
                LLustzcut[yi][xi]=ust
                LLxgridzcut[yi][xi]=x
                LLygridzcut[yi][xi]=y
                
                Val, Vect = np.linalg.eig(us22)
                
                #ordering: min/max eigval
                mMV = np.zeros((2,2))
                mMv = np.sort(Val)
                for orderindex in range(2):
                    for vindex in range(2):
                        if mMv[orderindex] == Val[vindex]:
                            mMV[:,orderindex] = Vect[:,vindex]

                eigvm = mMv[0]
                eigvM = mMv[1]
                eigVMx = mMV[1,1]
                eigVmx = mMV[1,0]
                eigVMy = mMV[0,1]
                eigVmy = mMV[0,0]

                isnan = True
                if eigVMx > 0 and eigVMy >= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)); isnan = False
                elif eigVMx <= 0 and eigVMy > 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + np.pi/2; isnan = False
                elif eigVMx < 0 and eigVMy <= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)) + np.pi; isnan = False
                elif eigVMx >= 0 and eigVMy < 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + 3*np.pi/2; isnan = False
                elif eigVMx == 0:
                    beta = np.pi/2; isnan = False
                elif eigVMy == 0:
                    beta = np.pi/2; isnan = False

                if isnan == False:
                    alphas = np.linspace(0,2,30)*np.pi
                    a = scale_factor*abs(eigvM)
                    b = scale_factor*abs(eigvm)
                    xs = x + np.cos(beta)*np.cos(alphas)*a - np.sin(beta)*np.sin(alphas)*b
                    ys = y + np.sin(beta)*np.cos(alphas)*a + np.cos(beta)*np.sin(alphas)*b
                    LinePlusx = [x + np.cos(beta)*a, x - np.cos(beta)*a]
                    LinePlusy = [y + np.sin(beta)*a, y - np.sin(beta)*a]
                    LineMinusx = [x + np.cos(beta+np.pi/2)*b, x - np.cos(beta+np.pi/2)*b]
                    LineMinusy = [y + np.sin(beta+np.pi/2)*b, y - np.sin(beta+np.pi/2)*b]
                    #if np.abs(eigvM) > np.abs(eigvm):
                    #    ax[0].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    #else:
                    #    ax[0].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    #ax[0].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if ust<0:
                        ax[0].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    if ust>0:
                        ax[0].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    if eigvM>0:
                        ax[0].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if eigvm>0:
                        ax[0].plot(LineMinusx,LineMinusy, 'k-', alpha = 1, linewidth = 3)
                    
            else:
                LLustzcut[yi][xi]=np.nan
                LLxgridzcut[yi][xi]=np.nan
                LLygridzcut[yi][xi]=np.nan

    # ycut (only in z,x (0,2))
    for zi in range(shapezyx[0]):
        for xi in range(shapezyx[2]):
            if Countzx[zi][xi] > Countmin:
                z,y,x = np.diag(Coordzx[zi][xi])
                x,y = x,z
                if DeviatoricType==3: # remove 1/2*log((a*b*c)**-3)
                    us = UMfromM(USzx[zi][xi])
                elif DeviatoricType==2: # remove log((a*b*c)**-3)
                    us = USfromS(USzx[zi][xi])
                elif DeviatoricType==1: # remove trace: a+b+c
                    us = SigdevfromSig(USzx[zi][xi])
                elif DeviatoricType==0: # keep as it is
                    us = USzx[zi][yi]
                else:
                    print('Error: provide DeviatoricType')
                    break
                us22 = np.zeros((2,2))
                us22[0,0] = us[0,0]
                us22[1,0] = us[2,0]
                us22[0,1] = us[0,2]
                us22[1,1] = us[2,2]
                ust = us[1,1]
                LLustycut[zi][xi]=ust
                
                LLxgridycut[zi][xi]=x
                LLygridycut[zi][xi]=y

                Val, Vect = np.linalg.eig(us22)

                #ordering: min/max eigval
                mMV = np.zeros((2,2))
                mMv = np.sort(Val)
                for orderindex in range(2):
                    for vindex in range(2):
                        if mMv[orderindex] == Val[vindex]:
                            mMV[:,orderindex] = Vect[:,vindex]

                eigvm = mMv[0]
                eigvM = mMv[1]
                eigVMx = mMV[1,1]
                eigVmx = mMV[1,0]
                eigVMy = mMV[0,1]
                eigVmy = mMV[0,0]

                isnan = True
                if eigVMx > 0 and eigVMy >= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)); isnan = False
                elif eigVMx <= 0 and eigVMy > 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + np.pi/2; isnan = False
                elif eigVMx < 0 and eigVMy <= 0:
                    beta = np.arctan(np.abs(eigVMy/eigVMx)) + np.pi; isnan = False
                elif eigVMx >= 0 and eigVMy < 0:
                    beta = np.arctan(np.abs(eigVMx/eigVMy)) + 3*np.pi/2; isnan = False
                elif eigVMx == 0:
                    beta = np.pi/2; isnan = False
                elif eigVMy == 0:
                    beta = np.pi/2; isnan = False

                if isnan == False:
                    alphas = np.linspace(0,2,30)*np.pi
                    a = scale_factor*abs(eigvM)
                    b = scale_factor*abs(eigvm)
                    xs = x + np.cos(beta)*np.cos(alphas)*a - np.sin(beta)*np.sin(alphas)*b
                    ys = y + np.sin(beta)*np.cos(alphas)*a + np.cos(beta)*np.sin(alphas)*b
                    LinePlusx = [x + np.cos(beta)*a, x - np.cos(beta)*a]
                    LinePlusy = [y + np.sin(beta)*a, y - np.sin(beta)*a]
                    LineMinusx = [x + np.cos(beta+np.pi/2)*b, x - np.cos(beta+np.pi/2)*b]
                    LineMinusy = [y + np.sin(beta+np.pi/2)*b, y - np.sin(beta+np.pi/2)*b]
                    #if np.abs(eigvM) > np.abs(eigvm):
                    #    ax[1].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    #else:
                    #    ax[1].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    #ax[1].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if ust<0:
                        ax[1].plot(xs,ys, 'b-', alpha = 1, linewidth = 3)
                    if ust>0:
                        ax[1].plot(xs,ys, 'r-', alpha = 1, linewidth = 3)
                    if eigvM>0:
                        ax[1].plot(LinePlusx,LinePlusy, 'k-', alpha = 1, linewidth = 3)
                    if eigvm>0:
                        ax[1].plot(LineMinusx,LineMinusy, 'k-', alpha = 1, linewidth = 3)
            else:
                LLustycut[zi][xi]=np.nan
                LLxgridycut[zi][xi]=np.nan
                LLygridycut[zi][xi]=np.nan
    ax[0].axis('equal')
    ax[1].axis('equal')
    ax[2].axis('equal')
    if len(np.shape(showscale))>0:
        for ra in showscale:
            ra= ra*scale_factor
            xs = np.cos(alphas)*ra
            ys = ra + np.sin(alphas)*ra
            ax[0].plot(xs,ys,'k')
            ax[1].plot(xs,ys,'k')
            ax[2].plot(xs,ys,'k')
    print('Normal min/max ax0',np.nanmin(LLustzcut), np.nanmax(LLustzcut))
    print('Normal min/max ax1',np.nanmin(LLustycut), np.nanmax(LLustycut))
    print('Normal min/max ax2',np.nanmin(LLustxcut), np.nanmax(LLustxcut))
    c=ax[2].pcolor(Grids[1], Grids[0], LLustxcut, vmin = vmin, vmax = vmax,cmap = 'seismic', shading='nearest', alpha=alpha)
    c=ax[1].pcolor(Grids[2], Grids[0], LLustycut, vmin = vmin, vmax = vmax,cmap = 'seismic', shading='nearest', alpha=alpha)
    c=ax[0].pcolor(Grids[2], Grids[1], LLustzcut, vmin = vmin, vmax = vmax,cmap = 'seismic', shading='nearest', alpha=alpha)
    cbar = fig.colorbar(c, ax=ax)
    
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])
    
    if returnfig!=None:
        return fig
    if dirsave !=None:
        if 'eps' in formt:
            extent = ax[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(dirsave+'_21.eps', bbox_inches=extent.expanded(2.0, 2.0), format=formt)
            extent = ax[1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(dirsave+'_20.eps', bbox_inches=extent.expanded(2.0, 2.0), format=formt)
            extent = ax[2].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(dirsave+'_10.eps', bbox_inches=extent.expanded(2.0, 2.0), format=formt)
        else:
            fig.savefig(dirsave+'.jpeg', format=formt, dpi=dpi)

    if colorbarlab != None:
        cbar.set_label(colorbarlab)

def USfromS(S):
    import numpy as np
    
    Val, Vect = np.linalg.eig(S)
    a,b,c = Val
    req = np.power(a*b*c, 1/3)
    Ua = np.log(a/req)
    Ub = np.log(b/req)
    Uc = np.log(c/req)

    # Deformation tensor eigbasis
    US = [[Ua,0,0],
          [0,Ub,0],
          [0,0,Uc]]

    # Deformation tensor zyx
    US = Vect@US@np.transpose(Vect)
    return US

def UMfromM(M):
    import numpy as np
    
    Val, Vect = np.linalg.eig(M)
    a,b,c = Val
    req = np.power(a*b*c, 1/3)
    Ua = 1/2*np.log(a/req)
    Ub = 1/2*np.log(b/req)
    Uc = 1/2*np.log(c/req)

    # Deformation tensor eigbasis
    UM = [[Ua,0,0],
          [0,Ub,0],
          [0,0,Uc]]

    # Deformation tensor zyx
    UM = Vect@UM@np.transpose(Vect)
    return UM

def SigdevfromSig(B):
    import numpy as np
    
    Val, Vect = np.linalg.eig(B)
    a,b,c = Val
    trace = a+b+c
    adev = a-trace/3
    bdev = b-trace/3
    cdev = c-trace/3

    # Deformation tensor eigbasis
    Bdev = [[adev,0,0],
            [0,bdev,0],
            [0,0,cdev]]

    # Deformation tensor zyx
    Bdev = Vect@Bdev@np.transpose(Vect)
    return Bdev

# Vector----------------





# Scalar----------------
def CutScalar3D(Grids, Scalar, zcut=False, ycut=False, xcut=False, figblocksize=7, vmin=0, vmax=1, cmap='jet', alpha=1, nameaxes=None):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    lenz=len(Grids[0])
    leny=len(Grids[1])
    lenx=len(Grids[2])
    if zcut == False:
        zcut = lenz//2
    if ycut == False:
        ycut = leny//2
    if xcut == False:
        xcut = lenx//2
   
    fig, ax = plt.subplots(ncols=3, figsize=(3*figblocksize, figblocksize), constrained_layout=True)
   
    ax[0].axis('equal')
    ax[1].axis('equal')
    ax[2].axis('equal')
    
    print('SumProj min/max ax0',np.nanmin(Scalar[xcut]), np.nanmax(Scalar[xcut]))
    print('SumProj min/max ax1',np.nanmin(Scalar[ycut]), np.nanmax(Scalar[ycut]))
    print('SumProj min/max ax2',np.nanmin(Scalar[zcut]), np.nanmax(Scalar[zcut]))
    c=ax[2].pcolor(Grids[1], Grids[0], Scalar[:,:,xcut], vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
    c=ax[1].pcolor(Grids[2], Grids[0], Scalar[:,ycut,:], vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
    c=ax[0].pcolor(Grids[2], Grids[1], Scalar[zcut,:,:], vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
    fig.colorbar(c, ax=ax)
    
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])


def ProjScalar3D(Grids, Scalar, figblocksize=7, vmin=0, vmax=1, cmap='jet', alpha=1, sumproj=False, nameaxes=None):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(ncols=3, figsize=(3*figblocksize, figblocksize), constrained_layout=True)
   
    ax[0].axis('equal')
    ax[1].axis('equal')
    ax[2].axis('equal')
    
    if sumproj:
        print('SumProj min/max ax0',np.nanmin(np.nansum(Scalar,0)), np.nanmax(np.nansum(Scalar,0)))
        print('SumProj min/max ax1',np.nanmin(np.nansum(Scalar,1)), np.nanmax(np.nansum(Scalar,1)))
        print('SumProj min/max ax2',np.nanmin(np.nansum(Scalar,2)), np.nanmax(np.nansum(Scalar,2)))
        c=ax[2].pcolor(Grids[1], Grids[0], np.nansum(Scalar,2), vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
        c=ax[1].pcolor(Grids[2], Grids[0], np.nansum(Scalar,1), vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
        c=ax[0].pcolor(Grids[2], Grids[1], np.nansum(Scalar,0), vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
        fig.colorbar(c, ax=ax)
    else:
        print('MeanProj min/max ax0',np.nanmin(np.nanmean(Scalar,0)), np.nanmax(np.nanmean(Scalar,0)))
        print('MeanProj min/max ax1',np.nanmin(np.nanmean(Scalar,1)), np.nanmax(np.nanmean(Scalar,1)))
        print('MeanProj min/max ax2',np.nanmin(np.nanmean(Scalar,2)), np.nanmax(np.nanmean(Scalar,2)))
        c=ax[2].pcolor(Grids[1], Grids[0], np.nanmean(Scalar,2), vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
        c=ax[1].pcolor(Grids[2], Grids[0], np.nanmean(Scalar,1), vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
        c=ax[0].pcolor(Grids[2], Grids[1], np.nanmean(Scalar,0), vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest', alpha=alpha)
        fig.colorbar(c, ax=ax)
    
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])


# PLot contacts

def PlotContact(Ctc,color='k',marker='o',markersize=5,linestyle='-', linewidth=2,nameaxes=None, ax=None,figblocksize=5):
    import numpy as np
    if len(np.shape(ax))==0:
        fig, ax = plt.subplots(1,3, figsize = (figblocksize*3, figblocksize), constrained_layout=True)
    for i in range(len(Ctc['x1'])):
        if Ctc['x1'][i]>0 and Ctc['x2'][i]>0:
            ax[0].plot([Ctc['x1'][i],Ctc['x2'][i]],[Ctc['y1'][i],Ctc['y2'][i]], marker='o',linestyle='-', color=color)
            ax[1].plot([Ctc['x1'][i],Ctc['x2'][i]],[Ctc['z1'][i],Ctc['z2'][i]], marker='o',linestyle='-', color=color)
            ax[2].plot([Ctc['y1'][i],Ctc['y2'][i]],[Ctc['z1'][i],Ctc['z2'][i]], marker='o',linestyle='-', color=color)
    for i in range(3):
        ax[i].axis('equal')
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])

def PlotT1(Ctc,color=['r','g'],marker='o',markersize=5,linestyle='-', linewidth=2,nameaxes=None, ax=None,figblocksize=5,cat = 1):
    import numpy as np
    
    if len(np.shape(ax))==0:
        fig, ax = plt.subplots(1,3, figsize = (figblocksize*3, figblocksize), constrained_layout=True)
        
    if len(np.shape(color)) > 0:
        ln = Ctc[0]
        for i in range(len(ln['lostCoord'])):
            if len(ln['newCoord'][i]) == cat:
                z1,y1,x1 = ln['lostCoord'][i][0][0]
                z2,y2,x2 = ln['lostCoord'][i][0][1]
                ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-', color=color[0])
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color[0])
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color[0])
                z1,y1,x1 = ln['newCoord'][i][0][0]
                z2,y2,x2 = ln['newCoord'][i][0][1]
                ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-', color=color[1])
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color[1])
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color[1])
        nl = Ctc[1]
        for i in range(len(nl['lostCoord'])):
            if len(nl['lostCoord'][i]) == cat:
                z1,y1,x1 = nl['lostCoord'][i][0][0]
                z2,y2,x2 = nl['lostCoord'][i][0][1]
                ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-', color=color[0])
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color[0])
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color[0])
                z1,y1,x1 = nl['newCoord'][i][0][0]
                z2,y2,x2 = nl['newCoord'][i][0][1]
                ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-', color=color[1])
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color[1])
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color[1])

    else:
        ln = Ctc[0]
        for i in range(len(ln['lostCoord'])):
            if len(ln['newCoord'][i]) == cat:
                z1,y1,x1 = ln['lostCoord'][i][0][0]
                z2,y2,x2 = ln['lostCoord'][i][0][1]
                p = ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-')
                color = p[0].get_color()
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color)
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color)
                z1,y1,x1 = ln['newCoord'][i][0][0]
                z2,y2,x2 = ln['newCoord'][i][0][1]
                ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-', color=color)
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color)
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color)
        nl = Ctc[1]
        for i in range(len(nl['lostCoord'])):
            if len(nl['lostCoord'][i]) == cat:
                z1,y1,x1 = nl['lostCoord'][i][0][0]
                z2,y2,x2 = nl['lostCoord'][i][0][1]
                p = ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-')
                color = p[0].get_color()
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color)
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color)
                z1,y1,x1 = nl['newCoord'][i][0][0]
                z2,y2,x2 = nl['newCoord'][i][0][1]
                ax[0].plot([x1,x2],[y1,y2], marker='o',linestyle='-', color=color)
                ax[1].plot([x1,x2],[z1,z2], marker='o',linestyle='-', color=color)
                ax[2].plot([y1,y2],[z1,z2], marker='o',linestyle='-', color=color)
    
    for i in range(3):
        ax[i].axis('equal')
    if nameaxes!=None:
        ax[0].set_xlabel(nameaxes[2]); ax[0].set_ylabel(nameaxes[1])
        ax[1].set_xlabel(nameaxes[2]); ax[1].set_ylabel(nameaxes[0])
        ax[2].set_xlabel(nameaxes[1]); ax[2].set_ylabel(nameaxes[0])





# Histogram ----------------

def Histogram(image, histtitle=False):
    """
    Plot a 1x1 grey value histogram.
    
    :param image: 3D image.
    :type image: numpy array
    :return: None
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    import skimage.exposure
    
    fig, ax = plt.subplots(ncols=1, figsize=(5, 5))
    hist, hist_centers = skimage.exposure.histogram(image)
    ax.plot(hist_centers, hist, lw=2)
    ax.set_yscale('log')
    if histtitle != False:
        ax.set_title(histtitle)
        


# Colormaps ----------------   
# From: https://github.com/delestro/rand_cmap

def RandomCmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True, GiveRange=None):
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


    if type not in ('bright', 'soft','dark','special'):
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
        
        
    # Generate dark colors, by limiting the RGB spectrum
    if type == 'dark':
        low = 0.05
        high = 0.4
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)    
        
        
    # Generate ranged colors, by limiting the RGB spectrum
    if type == 'special' and len(GiveRange)>1:
        low = GiveRange[0]
        high = GiveRange[1]
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

def LinCmap(vmin=0,vmax=10, first_color="b", last_color="r", verbose=True):
    """
    Creates a linear colormap for matplotlib.
    
    :param vmin: min value for the color-range
    :type vmin: float
    :param vmax: max value for the color-range
    :type vmax: float
    :param first_color_black: first color
    :type first_color_black: str or matplotlib color
    :param last_color_black: last color
    :type last_color_black: str or matplotlib color
    :param verbose: If True, prints the number of labels and shows the colormap.
    :type verbose: Bool
    :return: matplotlib colormap
    """
    
    import matplotlib.colors as mcol
    import matplotlib.cm as cm
    
    CM = mcol.LinearSegmentedColormap.from_list("lincmap",[first_color,last_color])
    cnorm = mcol.Normalize(vmin=vmin,vmax=vmax)
    cpick = cm.ScalarMappable(norm=cnorm,cmap=CM)
    cpick.set_array([])
    
    return cpick



def CutSerie(directory, name, imrange, cutzyx = [0,None], nameaxes=['z','y','x'], Ncolumns=1, title='',figblocksize=5,endread='.tiff',n0=3, cmap='gray', interpolation=None, crop=[0,-1]):

    import numpy as np
    import matplotlib.pyplot as plt
    from tifffile import imread
    from FoamQuant.Helper import strindex
    
    Nrows = (len(imrange)-1) // Ncolumns + 1
    #print(len(imrange) // Ncolumns, len(imrange) % Ncolumns, Nrows)
    
    fig, ax = plt.subplots(Nrows,Ncolumns, figsize = (figblocksize*Ncolumns, figblocksize*Nrows), constrained_layout=True)

    for imi in range(len(imrange)):
        rowi = imi // Ncolumns
        columi = imi % Ncolumns

        image = imread(directory+name+strindex(imrange[imi], n0=n0)+endread)
        shape = np.shape(image)
        
        if cutzyx[0] == 0:
            if cutzyx[1] != None:
                cut = image[cutzyx[1]]
            else:
                zmid = shape[0]//2
                cut = image[zmid]
            if Ncolumns==1:
                ax[rowi].set_xlabel(nameaxes[2]); ax[rowi].set_ylabel(nameaxes[1])
            elif Nrows==1:
                ax[columi].set_xlabel(nameaxes[2]); ax[columi].set_ylabel(nameaxes[1])
            else:
                ax[rowi][columi].set_xlabel(nameaxes[2]); ax[rowi][columi].set_ylabel(nameaxes[1])
        if cutzyx[0] == 1:
            if cutzyx[1] != None:
                cut = image[:,cutzyx[1]]
            else:
                ymid = shape[1]//2
                cut = image[:,ymid]
            if Ncolumns==1:
                ax[rowi].set_xlabel(nameaxes[2]); ax[rowi].set_ylabel(nameaxes[0])
            elif Nrows==1:
                ax[columi].set_xlabel(nameaxes[2]); ax[columi].set_ylabel(nameaxes[0])
            else:
                ax[rowi][columi].set_xlabel(nameaxes[2]); ax[rowi][columi].set_ylabel(nameaxes[0])
                
        if cutzyx[0] == 2:
            if cutzyx[1] != None:
                cut = image[:,:,cutzyx[1]]
            else:
                xmid = shape[2]//2
                cut = image[:,:,xmid]
            if Ncolumns==1:
                ax[rowi].set_xlabel(nameaxes[1]); ax[rowi].set_ylabel(nameaxes[0])
            elif Nrows==1:
                ax[columi].set_xlabel(nameaxes[1]); ax[columi].set_ylabel(nameaxes[0])
            else:
                ax[rowi][columi].set_xlabel(nameaxes[1]); ax[rowi][columi].set_ylabel(nameaxes[0])
                
        if Ncolumns==1:
            ax[rowi].imshow(cut[crop[0]:crop[1]], cmap=cmap, interpolation=interpolation)
            ax[rowi].set_title(title+strindex(imrange[imi], n0=n0))
        elif Nrows==1:
            ax[columi].imshow(cut[crop[0]:crop[1]], cmap=cmap, interpolation=interpolation)
            ax[columi].set_title(title+strindex(imrange[imi], n0=n0))
        else:    
            ax[rowi][columi].imshow(cut[crop[0]:crop[1]], cmap=cmap, interpolation=interpolation)
            ax[rowi][columi].set_title(title+strindex(imrange[imi], n0=n0))


def ProjVec3D(LCoord, Lv, cmap = 'viridis', colorbarlab=None, nameaxes=['z','y','x'], scale = 1, vmin = None, vmax = None):

    from matplotlib.colors import Normalize
    import numpy as np
    import matplotlib.pyplot as plt
    
    z, y, x = np.transpose(LCoord)
    vz, vy, vx = np.transpose(Lv)
    
    speed = np.sqrt(vx**2 + vy**2 + vz**2)
    if vmin == None:
        vmin = np.nanmin(speed)
    if vmax == None:
        vmax = np.nanmax(speed)

    norm = Normalize(vmin=vmin, vmax=vmax)
    
    fig, ax = plt.subplots(1,3, figsize = (5*3, 5), constrained_layout=True)

    q0=ax[0].quiver(
        x, y,
        vx, vy,speed,norm=norm, 
        pivot='mid',cmap=cmap,scale=scale)
    ax[0].set_xlabel(nameaxes[2])
    ax[0].set_ylabel(nameaxes[1])
    
    q0=ax[1].quiver(
        x, z,
        vx, vz,speed,norm=norm, 
        pivot='mid',cmap=cmap,scale=scale)
    ax[1].set_xlabel(nameaxes[2])
    ax[1].set_ylabel(nameaxes[0])
    
    q0=ax[2].quiver(
        y, z,
        vy, vz,speed,norm=norm, 
        pivot='mid',cmap=cmap,scale=scale)
    ax[2].set_xlabel(nameaxes[1])
    ax[2].set_ylabel(nameaxes[0])
    
    if colorbarlab != None:
        plt.colorbar(q0, ax=ax[2], label=colorbarlab)        





    