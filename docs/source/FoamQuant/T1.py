def LostContact(Pairs1, Pairs_tsl2, Reg1):
    from tqdm import tqdm
    print('LostContact')
    
    #zip the pairs
    ZIPt0 = list(zip(Pairs1['lab1'], Pairs1['lab2']))
    ZIPt1 = list(zip(Pairs_tsl2['lab1'], Pairs_tsl2['lab2']))
    
    LCoupleLost=[]
    LCoupleLostX=[];LCoupleLostY=[];LCoupleLostZ=[]
    
    #Check if lost
    for i0 in tqdm(range(len(ZIPt0))):
        if ZIPt0[i0][0] in Pairs_tsl2['lab1'] or ZIPt0[i0][0] in Pairs_tsl2['lab2']: # if label1 from couple is in the next image
            if ZIPt0[i0][1] in Pairs_tsl2['lab1'] or ZIPt0[i0][1] in Pairs_tsl2['lab2']: # if label2 from couple is in the next image
                if (ZIPt0[i0][0],ZIPt0[i0][1]) not in ZIPt1 and (ZIPt0[i0][1],ZIPt0[i0][0]) not in ZIPt1: # Couple lost from Couple list in the next image
                    LCoupleLost.append([ZIPt0[i0][0],ZIPt0[i0][1]])
    
    # Retrieve the coordinates
    print('>>> Retrieve the coordinates')
    LReglab = Reg1['lab']
    for ci in tqdm(range(len(LCoupleLost))):
        for regi in range(len(LReglab)):
            if LReglab[regi] == LCoupleLost[ci][0]:
                x1 = Reg1['x'][regi]
                y1 = Reg1['y'][regi]
                z1 = Reg1['z'][regi]
            if LReglab[regi] == LCoupleLost[ci][1]:
                x2 = Reg1['x'][regi]
                y2 = Reg1['y'][regi]
                z2 = Reg1['z'][regi]
        LCoupleLostX.append([x1,x2])
        LCoupleLostY.append([y1,y2])
        LCoupleLostZ.append([z1,z2])
                    
    Pack = {'clab':LCoupleLost, 'cx':LCoupleLostX, 'cy':LCoupleLostY, 'cz':LCoupleLostZ}
    print(len(LCoupleLostX)/len(ZIPt0)*100, '%')
    
    return Pack

def NewContact(Pairs1, Pairs_tsl2):
    from tqdm import tqdm
    
    print('NewContact')
    
    #zip the pairs
    ZIPt0 = list(zip(Pairs1['lab1'], Pairs1['lab2']))
    ZIPt1 = list(zip(Pairs_tsl2['lab1'], Pairs_tsl2['lab2']))
    
    LCoupleNewLab=[];LCoupleNewX=[];LCoupleNewY=[];LCoupleNewZ=[]
    #Check if New
    for i0 in tqdm(range(len(ZIPt1))):
        if ZIPt1[i0][0] in Pairs1['lab1'] or ZIPt1[i0][0] in Pairs1['lab2']: # if label1 from couple is in the next image
            if ZIPt1[i0][1] in Pairs1['lab1'] or ZIPt1[i0][1] in Pairs1['lab2']: # if label2 from couple is in the next image
                if (ZIPt1[i0][0],ZIPt1[i0][1]) not in ZIPt0 and (ZIPt1[i0][1],ZIPt1[i0][0]) not in ZIPt0: # Couple New from Couple list in the next image
                
                    LCoupleNewX.append([Pairs_tsl2['x1'][i0],Pairs_tsl2['x2'][i0]])
                    LCoupleNewY.append([Pairs_tsl2['y1'][i0],Pairs_tsl2['y2'][i0]])
                    LCoupleNewZ.append([Pairs_tsl2['z1'][i0],Pairs_tsl2['z2'][i0]])
                    LCoupleNewLab.append([Pairs_tsl2['lab1'][i0],Pairs_tsl2['lab2'][i0]])
                    
    Pack = {'clab':LCoupleNewLab,'cx':LCoupleNewX, 'cy':LCoupleNewY, 'cz':LCoupleNewZ}
    print(len(LCoupleNewX)/len(ZIPt0)*100, '%')
    
    return Pack

def LostNewContact_Batch(pairsdirname, pairstrldirname, savedirnamelost, regdirname, savedirnamenew, imrange, verbose=True):   
    import pandas as pd
    from FoamQuant.Helper import strindex
    from FoamQuant.FromContact import Read_ContactPair
    from FoamQuant.FromLabelled import Read_RegionProp
    from FoamQuant.T1 import LostContact
    from FoamQuant.T1 import NewContact
    import csv
    import numpy as np
    
    for i in range(len(imrange)-1):
        ind1 = imrange[i]
        ind2 = imrange[i+1]
        imifordir1 = strindex(ind1, n0=3)
        imifordir2 = strindex(ind2, n0=3)

        Pairs1 = Read_ContactPair(pairsdirname[1], pairsdirname[0], [ind1], verbose=True)
        Pairs_tsl2 = Read_ContactPair(pairstrldirname[1], pairstrldirname[0], [ind2], verbose=True, extended=True)
        Reg1 = Read_RegionProp(regdirname[1], regdirname[0], [ind1], verbose=True)
        
        LostCoord = LostContact(Pairs1[0], Pairs_tsl2[0], Reg1)
        NewCoord = NewContact(Pairs1[0], Pairs_tsl2[0])
        
        with open(savedirnamelost[0] + savedirnamelost[1] + imifordir1+'_'+imifordir2 + '.tsv', 'w', newline='') as csvfile:        
            writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            #first line
            writer.writerow(['lab1','lab2','x1','x2','y1','y2','z1','z2'])
            #data line by line
            for i in range(len(LostCoord['cx'])):
                writer.writerow([LostCoord['clab'][i][0], LostCoord['clab'][i][1],
                                 LostCoord['cx'][i][0], LostCoord['cx'][i][1],
                                 LostCoord['cy'][i][0], LostCoord['cy'][i][1],
                                 LostCoord['cz'][i][0], LostCoord['cz'][i][1]])
                        
        with open(savedirnamenew[0] + savedirnamenew[1] + imifordir1+'_'+imifordir2 + '.tsv', 'w', newline='') as csvfile:        
            writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            #first line
            writer.writerow(['lab1','lab2','x1','x2','y1','y2','z1','z2'])
            #data line by line
            for i in range(len(NewCoord['cx'])):
                writer.writerow([NewCoord['clab'][i][0], NewCoord['clab'][i][1],
                                 NewCoord['cx'][i][0], NewCoord['cx'][i][1],
                                 NewCoord['cy'][i][0], NewCoord['cy'][i][1],
                                 NewCoord['cz'][i][0], NewCoord['cz'][i][1]])



def Read_lostnew(readdirname, imrange, verbose=True):
    import pandas as pd
    from FoamQuant.Helper import strindex
    import csv
    import numpy as np
    
    LPack=[]
    for i in range(len(imrange)):
        ind1 = imrange[i]
        ind2 = imrange[i]+1
        imifordir1 = strindex(ind1, n0=3)
        imifordir2 = strindex(ind2, n0=3)
        
        if verbose:
            print(readdirname[0]+readdirname[1]+imifordir1+'_'+imifordir2+'.tsv')
        
        Pairs = pd.read_csv(readdirname[0]+readdirname[1]+imifordir1+'_'+imifordir2+'.tsv', sep='\t',engine="python",  quoting=csv.QUOTE_NONE)
        x1 = np.asarray(Pairs['x1'])
        x2 = np.asarray(Pairs['x2'])
        y1 = np.asarray(Pairs['y1'])
        y2 = np.asarray(Pairs['y2'])
        z1 = np.asarray(Pairs['z1'])
        z2 = np.asarray(Pairs['z2'])
        lab1 = np.asarray(Pairs['lab1'])
        lab2 = np.asarray(Pairs['lab2'])
        
        Pack = {'lab1':lab1,'lab2':lab2,'x1':x1,'x2':x2,'y1':y1,'y2':y2,'z1':z1,'z2':z2}
        LPack.append(Pack)
    return LPack



def DetectT1_LostToNew(cont, lost, new, lcoord1,lcoord2, ncoord1,ncoord2):
    import numpy as np
    
    Llost_T1=[]
    Lnew_T1=[]
    Llost_T1_C=[]
    Lnew_T1_C=[]
    lostnoT1=[]
    lostnoT1coord=[]
    
    from tqdm import tqdm
    
    cont1 = cont['lab1']
    cont2 = cont['lab2']
    
    for li in tqdm(range(len(lost))): # lost loop
        
        lost_T1=[]
        new_T1=[]
        lost_T1_C=[]
        new_T1_C=[]
        
        l1 = lost[li][0]
        l2 = lost[li][1]
        
        neis1=[]; neis2=[]
        # get neighbours of the two lost labels ///
        for conti in range(len(cont1)):  # all contacts loop
            if cont1[conti] == l1 and cont2[conti] != l2: # all bbls in contact with l1, except bbl l2 (first bbl of the contact)
                neis1.append(cont2[conti])
            if cont2[conti] == l1 and cont1[conti] != l2: # all bbls in contact with l1, except bbl l2 (second bbl of the contact)
                neis1.append(cont1[conti])
            
            if cont1[conti] == l2 and cont2[conti] != l1: # all bbls in contact with l2, except bbl l1 (first bbl of the contact)
                neis2.append(cont2[conti])
            if cont2[conti] == l2 and cont1[conti] != l1: # all bbls in contact with l2, except bbl l1 (second bbl of the contact)
                neis2.append(cont1[conti])
                
        if len(neis1)==0 or len(neis2)==0:
            print('neighbours not found')
        
        # check if neighbours 1 are creating a new contact with neighbours 2 ///
        for ni in range(len(new)): # new loop
            for i1 in range(len(neis1)): # neighbours of l1 loop
                if neis1[i1] in new[ni]: # if l1 is included in a new contact
                    
                    if neis1[i1] == new[ni][0]: # is it the first label of the contact?
                        if new[ni][1] in neis2: # is the other label of the contact in the neighbours2 list?
                            if new[ni] not in new_T1: # if not already found: check new_contact list
                                lost_T1.append(lost[li])
                                new_T1.append(new[ni])
                                lost_T1_C.append([lcoord1[li],lcoord2[li]])
                                new_T1_C.append([ncoord1[ni],ncoord2[ni]])
                    else: # then it is the second label of the contact
                        if new[ni][0] in neis2: # is the other label of the contact in the neighbours2 list?
                            if new[ni] not in new_T1: # if not already found: check new_contact list
                                lost_T1.append(lost[li])
                                new_T1.append(new[ni])
                                lost_T1_C.append([lcoord1[li],lcoord2[li]])
                                new_T1_C.append([ncoord1[ni],ncoord2[ni]])
        Llost_T1.append(lost_T1)
        Lnew_T1.append(new_T1)
        Llost_T1_C.append(lost_T1_C)
        Lnew_T1_C.append(new_T1_C)
        if len(lost_T1)==0:
            lostnoT1.append(lost[li])
            lostnoT1coord.append([lcoord1[li],lcoord2[li]])
        
    return Llost_T1, Lnew_T1, Llost_T1_C, Lnew_T1_C, lostnoT1, lostnoT1coord




def DetectT1_NewToLost(cont, lost, new, lcoord1,lcoord2, ncoord1,ncoord2):
    import numpy as np
    
    Llost_T1=[]
    Lnew_T1=[]
    Llost_T1_C=[]
    Lnew_T1_C=[]
    newnoT1=[]
    newnoT1coord=[]
    
    from tqdm import tqdm
    
    cont1 = cont['lab1']
    cont2 = cont['lab2']
    
    for ni in tqdm(range(len(new))): # new loop
        
        lost_T1=[]
        new_T1=[]
        lost_T1_C=[]
        new_T1_C=[]
        
        n1 = new[ni][0]
        n2 = new[ni][1]
        
        neis1=[]; neis2=[]
        # get neighbours of the two new labels ///
        for conti in range(len(cont1)):  # all contacts loop
            if cont1[conti] == n1 and cont2[conti] != n2: # all bbls in contact with l1, except bbl l2 (first bbl of the contact)
                neis1.append(cont2[conti])
            if cont2[conti] == n1 and cont1[conti] != n2: # all bbls in contact with l1, except bbl l2 (second bbl of the contact)
                neis1.append(cont1[conti])
            
            if cont1[conti] == n2 and cont2[conti] != n1: # all bbls in contact with l2, except bbl l1 (first bbl of the contact)
                neis2.append(cont2[conti])
            if cont2[conti] == n2 and cont1[conti] != n1: # all bbls in contact with l2, except bbl l1 (second bbl of the contact)
                neis2.append(cont1[conti])
                
        if len(neis1)==0 or len(neis2)==0:
            print('neighbours not found')
        
        # check if neighbours 1 have lost a contact with neighbours 2 ///
        for li in range(len(lost)): # new loop
            for i1 in range(len(neis1)): # neighbours of l1 loop
                if neis1[i1] in lost[li]: # if l1 is included in a new contact
                    
                    if neis1[i1] == lost[li][0]: # is it the first label of the contact?
                        if lost[li][1] in neis2: # is the other label of the contact in the neighbours2 list?
                            if lost[li] not in lost_T1: # if not already found: check lost_contact list
                                lost_T1.append(lost[li])
                                new_T1.append(new[ni])
                                lost_T1_C.append([lcoord1[li],lcoord2[li]])
                                new_T1_C.append([ncoord1[ni],ncoord2[ni]])
                    else: # then it is the second label of the contact
                        if lost[li][0] in neis2: # is the other label of the contact in the neighbours2 list?
                            if lost[li] not in lost_T1: # if not already found: check lost_contact list
                                lost_T1.append(lost[li])
                                new_T1.append(new[ni])
                                lost_T1_C.append([lcoord1[li],lcoord2[li]])
                                new_T1_C.append([ncoord1[ni],ncoord2[ni]])
        Llost_T1.append(lost_T1)
        Lnew_T1.append(new_T1)
        Llost_T1_C.append(lost_T1_C)
        Lnew_T1_C.append(new_T1_C)
        if len(new_T1)==0:
            newnoT1.append(new[ni])
            newnoT1coord.append([ncoord1[ni],ncoord2[ni]])
        
    return Llost_T1, Lnew_T1, Llost_T1_C, Lnew_T1_C, newnoT1, newnoT1coord




                
                
                
def DetectT1_Batch(readdirnameall, readdirnamelost, readdirnamenew, namesave, dirsave, imrange, verbose=False, n0=3):
    
    print(imrange)
    
    import numpy as np 
    import pandas as pd
    import pickle as pkl
    from FoamQuant.Helper import strindex
    from FoamQuant.FromContact import Read_ContactPair
    import os
    
    #Check directory
    isExist = os.path.exists(dirsave)
    print('Path exist:', isExist)
    if not isExist:
        print('Saving path does not exist:\n', isExist)
        return
    
    #Batch loop
    for imi in imrange:
        
        # string index
        imifordir = strindex(imi, n0)
        
        # read lost and new contacts
        if verbose:
            print(readdirnamelost[0]+readdirnamelost[1]+strindex(imi, n0)+'_'+strindex(imi+1, n0))
            print(readdirnamenew[0]+readdirnamenew[1]+strindex(imi, n0)+'_'+strindex(imi+1, n0))
        lostcont = Read_lostnew(readdirnamelost, [imi], verbose=False)
        newcont = Read_lostnew(readdirnamenew, [imi], verbose=False)
        lostcont = lostcont[0]
        newcont = newcont[0]
        
        ncoord1=[]; ncoord2=[]; nc=[]
        for ci in range(len(newcont['lab2'])): #new contact loop
            nc.append([newcont['lab1'][ci],newcont['lab2'][ci]])
            ncoord1.append([newcont['z1'][ci],newcont['y1'][ci],newcont['x1'][ci]])
            ncoord2.append([newcont['z2'][ci],newcont['y2'][ci],newcont['x2'][ci]])
        lcoord1=[]; lcoord2=[]; lc=[]
        for ci in range(len(lostcont['lab2'])): #lost contact loop
            lc.append([lostcont['lab1'][ci],lostcont['lab2'][ci]])
            lcoord1.append([lostcont['z1'][ci],lostcont['y1'][ci],lostcont['x1'][ci]])
            lcoord2.append([lostcont['z2'][ci],lostcont['y2'][ci],lostcont['x2'][ci]])
            
        # read all contacts
        if verbose:
            print(readdirnameall[0]+readdirnameall[1]+strindex(imi, n0))
        cont = Read_ContactPair(readdirnameall[1], readdirnameall[0], [imi], verbose=False, endread='.tsv', n0=3, extended=False)
        cont = cont[0]
        
        # from lost contacts, find new contacts
        lost_T1, new_T1, lost_T1_C, new_T1_C, lostnoT1, lostnoT1coord = DetectT1_LostToNew(cont, lc, nc, lcoord1,lcoord2, ncoord1,ncoord2)
        Pack = {"lost": lost_T1, "new": new_T1, "lostCoord": lost_T1_C, "newCoord": new_T1_C, "lostnoT1":lostnoT1, "lostnoT1coord":lostnoT1coord}
        #Save as pickle
        with open(dirsave + namesave +imifordir+ '_LostNew.pkl','wb') as file:
            pkl.dump(Pack, file, pkl.HIGHEST_PROTOCOL)
        # from new contacts, find lost contacts
        lost_T1, new_T1, lost_T1_C, new_T1_C, newnoT1, newnoT1coord = DetectT1_NewToLost(cont, lc, nc, lcoord1,lcoord2, ncoord1,ncoord2)
        Pack = {"lost": lost_T1, "new": new_T1, "lostCoord": lost_T1_C, "newCoord": new_T1_C, "newnoT1":newnoT1, "newnoT1coord":newnoT1coord}
        #Save as pickle
        with open(dirsave + namesave +imifordir+ '_NewLost.pkl','wb') as file:
            pkl.dump(Pack, file, pkl.HIGHEST_PROTOCOL)
            
        if verbose:
                print(namesave +imifordir+': done')


def ReadT1(readdirname,imrange, verbose=False, n0=3):
    
    import numpy as np
    import pickle as pkl
    from FoamQuant.Helper import strindex
    
    import pandas as pd
    
    LT1_nl=[]
    for i in range(len(imrange)):
        ind = imrange[i]
        imifordir = strindex(ind, n0=3)
        
        if verbose:
            print(readdirname[0]+readdirname[1]+imifordir+'_NewLost'+'.tsv')
        
        with open(readdirname[0]+readdirname[1]+imifordir+'_NewLost'+'.pkl','rb') as f:
                LT1_nl.append(pkl.load(f))
                
    LT1_ln=[]
    for i in range(len(imrange)):
        ind = imrange[i]
        imifordir = strindex(ind, n0=3)
        
        if verbose:
            print(readdirname[0]+readdirname[1]+imifordir+'_LostNew'+'.tsv')
        
        with open(readdirname[0]+readdirname[1]+imifordir+'_LostNew'+'.pkl','rb') as f:
                LT1_ln.append(pkl.load(f))
    
    return [LT1_nl,LT1_ln]