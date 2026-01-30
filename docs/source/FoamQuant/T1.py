def LostContact(Pairs1, Pairs_tsl2, Reg1):
    """
    LostContact
    -----------
    Detects contacts that disappear between two consecutive time steps
    while both labels are still present.

    Parameters
    ----------
    Pairs1 : dict
        Contact pairs at time t (keys: 'lab1', 'lab2').
    Pairs_tsl2 : dict
        Contact pairs at time t+1 after tracking (keys: 'lab1', 'lab2').
    Reg1 : dict
        Region properties at time t (keys: 'lab', 'x', 'y', 'z').

    Outputs
    -------
    Pack : dict
        Dictionary containing:
        - 'clab' : list of lost contact label pairs
        - 'cx'   : x coordinates of lost contacts
        - 'cy'   : y coordinates of lost contacts
        - 'cz'   : z coordinates of lost contacts
    """

    from tqdm import tqdm
    print('LostContact')

    ZIPt0 = list(zip(Pairs1['lab1'], Pairs1['lab2']))
    ZIPt1 = list(zip(Pairs_tsl2['lab1'], Pairs_tsl2['lab2']))

    LCoupleLost = []
    LCoupleLostX = []
    LCoupleLostY = []
    LCoupleLostZ = []

    for i0 in tqdm(range(len(ZIPt0))):
        if ZIPt0[i0][0] in Pairs_tsl2['lab1'] or ZIPt0[i0][0] in Pairs_tsl2['lab2']:
            if ZIPt0[i0][1] in Pairs_tsl2['lab1'] or ZIPt0[i0][1] in Pairs_tsl2['lab2']:
                if (ZIPt0[i0][0], ZIPt0[i0][1]) not in ZIPt1 and \
                   (ZIPt0[i0][1], ZIPt0[i0][0]) not in ZIPt1:
                    LCoupleLost.append([ZIPt0[i0][0], ZIPt0[i0][1]])

    print('>>> Retrieve the coordinates')
    LReglab = Reg1['lab']
    for ci in tqdm(range(len(LCoupleLost))):
        for regi in range(len(LReglab)):
            if LReglab[regi] == LCoupleLost[ci][0]:
                x1, y1, z1 = Reg1['x'][regi], Reg1['y'][regi], Reg1['z'][regi]
            if LReglab[regi] == LCoupleLost[ci][1]:
                x2, y2, z2 = Reg1['x'][regi], Reg1['y'][regi], Reg1['z'][regi]

        LCoupleLostX.append([x1, x2])
        LCoupleLostY.append([y1, y2])
        LCoupleLostZ.append([z1, z2])

    Pack = {
        'clab': LCoupleLost,
        'cx': LCoupleLostX,
        'cy': LCoupleLostY,
        'cz': LCoupleLostZ
    }

    print(len(LCoupleLostX) / len(ZIPt0) * 100, '%')
    return Pack

def NewContact(Pairs1, Pairs_tsl2):
    """
    NewContact
    ----------
    Detects contacts that appear between two consecutive time steps.

    Parameters
    ----------
    Pairs1 : dict
        Contact pairs at time t.
    Pairs_tsl2 : dict
        Contact pairs at time t+1 after tracking (extended format).

    Outputs
    -------
    Pack : dict
        Dictionary containing:
        - 'clab' : list of new contact label pairs
        - 'cx'   : x coordinates
        - 'cy'   : y coordinates
        - 'cz'   : z coordinates
    """

    from tqdm import tqdm
    print('NewContact')

    ZIPt0 = list(zip(Pairs1['lab1'], Pairs1['lab2']))
    ZIPt1 = list(zip(Pairs_tsl2['lab1'], Pairs_tsl2['lab2']))

    LCoupleNewLab = []
    LCoupleNewX = []
    LCoupleNewY = []
    LCoupleNewZ = []

    for i0 in tqdm(range(len(ZIPt1))):
        if ZIPt1[i0][0] in Pairs1['lab1'] or ZIPt1[i0][0] in Pairs1['lab2']:
            if ZIPt1[i0][1] in Pairs1['lab1'] or ZIPt1[i0][1] in Pairs1['lab2']:
                if (ZIPt1[i0][0], ZIPt1[i0][1]) not in ZIPt0 and \
                   (ZIPt1[i0][1], ZIPt1[i0][0]) not in ZIPt0:
                    LCoupleNewX.append([Pairs_tsl2['x1'][i0], Pairs_tsl2['x2'][i0]])
                    LCoupleNewY.append([Pairs_tsl2['y1'][i0], Pairs_tsl2['y2'][i0]])
                    LCoupleNewZ.append([Pairs_tsl2['z1'][i0], Pairs_tsl2['z2'][i0]])
                    LCoupleNewLab.append([Pairs_tsl2['lab1'][i0], Pairs_tsl2['lab2'][i0]])

    Pack = {
        'clab': LCoupleNewLab,
        'cx': LCoupleNewX,
        'cy': LCoupleNewY,
        'cz': LCoupleNewZ
    }

    print(len(LCoupleNewX) / len(ZIPt0) * 100, '%')
    return Pack

def LostNewContact_Batch(pairsdirname, pairstrldirname,
                         savedirnamelost, regdirname,
                         savedirnamenew, imrange, verbose=True):
    """
    LostNewContact_Batch
    -------------------
    Batch computation of lost and new contacts between consecutive images
    and saves them as TSV files.

    Parameters
    ----------
    pairsdirname : tuple
        Directory info for original contact pairs.
    pairstrldirname : tuple
        Directory info for tracked contact pairs.
    savedirnamelost : tuple
        Output directory for lost contacts.
    regdirname : tuple
        Directory info for region properties.
    savedirnamenew : tuple
        Output directory for new contacts.
    imrange : iterable of int
        Image indices to process.
    verbose : bool, optional
        If True, prints progress.

    Outputs
    -------
    None
        Writes TSV files for lost and new contacts.
    """

    import csv
    import numpy as np
    from FoamQuant.Helper import strindex
    from FoamQuant.FromContact import Read_ContactPair
    from FoamQuant.FromLabelled import Read_RegionProp

    for i in range(len(imrange) - 1):
        ind1, ind2 = imrange[i], imrange[i + 1]
        imifordir1 = strindex(ind1, n0=3)
        imifordir2 = strindex(ind2, n0=3)

        Pairs1 = Read_ContactPair(pairsdirname[1], pairsdirname[0], [ind1], verbose=True)
        Pairs_tsl2 = Read_ContactPair(pairstrldirname[1], pairstrldirname[0],
                                      [ind2], verbose=True, extended=True)
        Reg1 = Read_RegionProp(regdirname[1], regdirname[0], [ind1], verbose=True)

        LostCoord = LostContact(Pairs1[0], Pairs_tsl2[0], Reg1)
        NewCoord = NewContact(Pairs1[0], Pairs_tsl2[0])

        for Pack, savedir in zip([LostCoord, NewCoord],
                                 [savedirnamelost, savedirnamenew]):

            with open(savedir[0] + savedir[1] +
                      imifordir1 + '_' + imifordir2 + '.tsv',
                      'w', newline='') as csvfile:

                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow(['lab1', 'lab2', 'x1', 'x2', 'y1', 'y2', 'z1', 'z2'])

                for i in range(len(Pack['cx'])):
                    writer.writerow([
                        Pack['clab'][i][0], Pack['clab'][i][1],
                        Pack['cx'][i][0], Pack['cx'][i][1],
                        Pack['cy'][i][0], Pack['cy'][i][1],
                        Pack['cz'][i][0], Pack['cz'][i][1]
                    ])

def Read_lostnew(readdirname, imrange, verbose=True):
    """
    Read_lostnew
    ------------
    Reads lost or new contact TSV files into Python dictionaries.

    Parameters
    ----------
    readdirname : tuple
        Directory path information.
    imrange : iterable of int
        Image indices.
    verbose : bool, optional
        If True, prints file paths.

    Outputs
    -------
    LPack : list of dict
        One dictionary per image containing labels and coordinates.
    """

    import pandas as pd
    import csv
    import numpy as np
    from FoamQuant.Helper import strindex

    LPack = []
    for ind1 in imrange:
        ind2 = ind1 + 1
        imifordir1 = strindex(ind1, n0=3)
        imifordir2 = strindex(ind2, n0=3)

        if verbose:
            print(readdirname[0] + readdirname[1] +
                  imifordir1 + '_' + imifordir2 + '.tsv')

        Pairs = pd.read_csv(
            readdirname[0] + readdirname[1] +
            imifordir1 + '_' + imifordir2 + '.tsv',
            sep='\t', quoting=csv.QUOTE_NONE
        )

        Pack = {k: np.asarray(Pairs[k]) for k in
                ['lab1', 'lab2', 'x1', 'x2', 'y1', 'y2', 'z1', 'z2']}

        LPack.append(Pack)

    return LPack
