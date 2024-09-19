Elementary plastic rearrangement (T1) detection
===============================================

In this jupyternotebook we are going to detect T1 events between each
individual bubble in a set of images.

.. code:: ipython3

    # Import Libraries
    import numpy as np
    import matplotlib.pyplot as plt
    from tifffile import imread, imsave
    import skimage.measure 
    import pickle as pkl
    import os
    import spam
    
    # Import FoamQuant library
    from FoamQuant import *
    
    # Set matplotlib default font size
    plt.rc('font', size=25) 

.. code:: ipython3

    # Create the processing pipeline
    ProcessPipeline = ['P4_BubbleSegmented','P5_BubbleNoEdge','Q7_RegProp','Q8_Tracking', 
                       'Q9_Contacts','Q10_Translated','Q11_LostNewContacts','Q12_T1']
    
    for Pi in ProcessPipeline:
        if  os.path.exists(Pi):
            print('path already exist:',Pi)
        else:
            print('Created:',Pi)
            os.mkdir(Pi)


.. parsed-literal::

    path already exist: P4_BubbleSegmented
    path already exist: P5_BubbleNoEdge
    path already exist: Q7_RegProp
    path already exist: Q8_Tracking
    path already exist: Q9_Contacts
    path already exist: Q10_Translated
    path already exist: Q11_LostNewContacts
    path already exist: Q12_T1


A) Type of imput data
---------------------

The images we are extracting the individual bubble stress tensor are
bubble-segmented images with removed labels on the edge of the images.

.. code:: ipython3

    # Read/Save image names and directories
    dirread = ProcessPipeline[1]+'/'
    nameread = 'BubbleNoEdge_'
    imrange = [1,2,3,4,5]

.. code:: ipython3

    # Read the first image of the series
    Lab = imread(dirread+nameread+strindex(imrange[0], 3)+'.tiff')
    
    # Create a random colormap
    rcmap = RandomCmap(5000)
    
    # Show a 3D-cut view of the volume
    Cut3D(Lab, 
          showcuts=True,
          cmap=rcmap, 
          interpolation='nearest', 
          figblocksize=7,           # tune this parrameter to change the figure size
          zcut=30,                  # tune this parrameter to change the orthogonal z-cut position
          ycut=False,               # tune this parrameter to change the orthogonal y-cut position
          xcut=False)               # tune this parrameter to change the orthogonal x-cut position


.. parsed-literal::

    Number of labels: 5000



.. image:: Jupy_FoamQuant_T1_detection_files/Jupy_FoamQuant_T1_detection_7_1.png



.. image:: Jupy_FoamQuant_T1_detection_files/Jupy_FoamQuant_T1_detection_7_2.png


B) Track the bubbles
--------------------

.. code:: ipython3

    # Read/Save names and directories
    nameread = 'BubbleNoEdge_'
    namesave = 'Props_'
    dirread = ProcessPipeline[1]+'/'
    dirsave = ProcessPipeline[2]+'/'
    # Images indexes
    imrange = [1,2,3,4,5]

.. code:: ipython3

    # Region properties
    RegionProp_Batch(nameread,
                     namesave,
                     dirread,
                     dirsave,
                     imrange,
                     verbose=True,
                     field=[40,220,40,220,40,220], # tune this parrameter if you wish
                     endread='.tiff', 
                     endsave='.tsv')


.. parsed-literal::

    Path exist: True
    Props_001: done
    Props_002: done
    Props_003: done
    Props_004: done
    Props_005: done


C) Bubble tracking
------------------

.. code:: ipython3

    # Read/Save image names and directories
    nameread = 'Props_'
    namesave = 'Tracking_'
    dirread = ProcessPipeline[2]+'/'
    dirsave = ProcessPipeline[3]+'/'
    # Images indexes
    imrange = [1,2,3,4,5]

.. code:: ipython3

    # Tracking
    LabelTracking_Batch(nameread, 
                        namesave, 
                        dirread, 
                        dirsave,
                        imrange,
                        verbose=False,
                        endread='.tsv',
                        endsave='.tsv',
                        n0=3,
                        searchbox=[-10,10,-10,10,-10,10],   # tune this parrameter if you wish
                        Volpercent=0.05)              # tune this parrameter if you wish


.. parsed-literal::

    Path exist: True


.. parsed-literal::

    100%|██████████| 504/504 [00:00<00:00, 812.91it/s]


.. parsed-literal::

    Lost tracking: 13 2.579365079365079 %


.. parsed-literal::

    100%|██████████| 506/506 [00:00<00:00, 792.46it/s]


.. parsed-literal::

    Lost tracking: 9 1.7786561264822136 %


.. parsed-literal::

    100%|██████████| 511/511 [00:00<00:00, 803.22it/s]


.. parsed-literal::

    Lost tracking: 7 1.36986301369863 %


.. parsed-literal::

    100%|██████████| 514/514 [00:00<00:00, 736.20it/s]

.. parsed-literal::

    Lost tracking: 16 3.11284046692607 %


.. parsed-literal::

    


D) Contact topology
-------------------

.. code:: ipython3

    # Read/Save image names and directories
    nameread = 'BubbleSeg_'
    nameread_noedge = 'BubbleNoEdge_'
    namesave = 'Contact_'
    dirread = ProcessPipeline[0]+'/'
    dirread_noedge = ProcessPipeline[1]+'/'
    dirsave = ProcessPipeline[4]+'/'
    # Images indexes
    imrange = [1,2,3,4,5]

.. code:: ipython3

    GetContacts_Batch(nameread,
                      nameread_noedge,
                      namesave,
                      dirread,
                      dirread_noedge,
                      dirsave,
                      imrange,
                      verbose=True)


.. parsed-literal::

    Path exist: True
    Contact_001: done
    Contact_002: done
    Contact_003: done
    Contact_004: done
    Contact_005: done


E) Translate contact pairs
--------------------------

.. code:: ipython3

    # Read/Save image names and directories
    nameread = 'Contact_pair_'
    nameread_track = 'Tracking_'
    namesave = 'Translated_pair_'
    dirread = ProcessPipeline[4]+'/'
    dirread_track = ProcessPipeline[3]+'/'
    dirsave = ProcessPipeline[5]+'/'
    # Images indexes
    imrange = [1,2,3,4,5]

.. code:: ipython3

    Translate_Pairs_Batch(nameread, 
                          namesave, 
                          nameread_track, 
                          dirread_track, 
                          dirread, 
                          dirsave, 
                          imrange, 
                          endsave='.tsv', 
                          n0=3)


.. parsed-literal::

    Tracking_001_002 : done
    Contact_pair_002: done


.. parsed-literal::

    100%|██████████| 5026/5026 [00:00<00:00, 5483.67it/s]


.. parsed-literal::

    Tracking_002_003 : done
    Contact_pair_003: done


.. parsed-literal::

    100%|██████████| 5033/5033 [00:00<00:00, 5498.18it/s]


.. parsed-literal::

    Tracking_003_004 : done
    Contact_pair_004: done


.. parsed-literal::

    100%|██████████| 5005/5005 [00:00<00:00, 5525.60it/s]


.. parsed-literal::

    Tracking_004_005 : done
    Contact_pair_005: done


.. parsed-literal::

    100%|██████████| 5037/5037 [00:00<00:00, 5455.60it/s]


F) Detect lost and new contacts
-------------------------------

.. code:: ipython3

    # Read/Save image names and directories
    pairsdirname = [ProcessPipeline[4]+'/','Contact_pair_']
    pairstrldirname = [ProcessPipeline[5]+'/','Translated_pair_']
    regdirname = [ProcessPipeline[2]+'/','Props_']
    
    savedirnamelost = [ProcessPipeline[6]+'/','Lost_']
    savedirnamenew = [ProcessPipeline[6]+'/','New_']
    
    # Images indexes
    imrange = [1,2,3,4,5]

.. code:: ipython3

    LostNewContact_Batch(pairsdirname, 
                         pairstrldirname, 
                         savedirnamelost, 
                         regdirname, 
                         savedirnamenew, 
                         imrange, 
                         verbose=True)


.. parsed-literal::

    Contact_pair_001: done
    Translated_pair_002: done
    Props_001: done
    LostContact


.. parsed-literal::

    100%|██████████| 4960/4960 [00:00<00:00, 16217.84it/s]


.. parsed-literal::

    >>> Retrieve the coordinates


.. parsed-literal::

    100%|██████████| 66/66 [00:00<00:00, 5437.30it/s]


.. parsed-literal::

    1.3306451612903225 %
    NewContact


.. parsed-literal::

    100%|██████████| 5026/5026 [00:00<00:00, 16435.79it/s]


.. parsed-literal::

    1.2096774193548387 %
    Contact_pair_002: done
    Translated_pair_003: done
    Props_002: done
    LostContact


.. parsed-literal::

    100%|██████████| 5026/5026 [00:00<00:00, 16132.90it/s]


.. parsed-literal::

    >>> Retrieve the coordinates


.. parsed-literal::

    100%|██████████| 68/68 [00:00<00:00, 5437.28it/s]


.. parsed-literal::

    1.3529645841623557 %
    NewContact


.. parsed-literal::

    100%|██████████| 5033/5033 [00:00<00:00, 16098.69it/s]


.. parsed-literal::

    1.1142061281337048 %
    Contact_pair_003: done
    Translated_pair_004: done
    Props_003: done
    LostContact


.. parsed-literal::

    100%|██████████| 5033/5033 [00:00<00:00, 15592.16it/s]


.. parsed-literal::

    >>> Retrieve the coordinates


.. parsed-literal::

    100%|██████████| 83/83 [00:00<00:00, 5058.00it/s]


.. parsed-literal::

    1.6491158354857938 %
    NewContact


.. parsed-literal::

    100%|██████████| 5005/5005 [00:00<00:00, 15237.07it/s]


.. parsed-literal::

    1.7683290284124777 %
    Contact_pair_004: done
    Translated_pair_005: done
    Props_004: done
    LostContact


.. parsed-literal::

    100%|██████████| 5005/5005 [00:00<00:00, 16167.55it/s]


.. parsed-literal::

    >>> Retrieve the coordinates


.. parsed-literal::

    100%|██████████| 46/46 [00:00<00:00, 5329.19it/s]


.. parsed-literal::

    0.919080919080919 %
    NewContact


.. parsed-literal::

    100%|██████████| 5037/5037 [00:00<00:00, 16093.94it/s]

.. parsed-literal::

    1.118881118881119 %


.. parsed-literal::

    


.. code:: ipython3

    Lost = Read_lostnew(savedirnamelost, 
                        imrange[:-1], 
                        verbose=True)
    
    New = Read_lostnew(savedirnamenew, 
                        imrange[:-1], 
                        verbose=True)


.. parsed-literal::

    Q11_LostNewContacts/Lost_001_002.tsv
    Q11_LostNewContacts/Lost_002_003.tsv
    Q11_LostNewContacts/Lost_003_004.tsv
    Q11_LostNewContacts/Lost_004_005.tsv
    Q11_LostNewContacts/New_001_002.tsv
    Q11_LostNewContacts/New_002_003.tsv
    Q11_LostNewContacts/New_003_004.tsv
    Q11_LostNewContacts/New_004_005.tsv


.. code:: ipython3

    # Show all the lost and new contact at the first time step
    fig, ax = plt.subplots(1,3, figsize = (5*3, 5), constrained_layout=True)
    PlotContact(Lost[0],color='r',ax=ax, nameaxes=['z','y','x'])
    PlotContact(New[0],color='g',ax=ax, nameaxes=['z','y','x'])



.. image:: Jupy_FoamQuant_T1_detection_files/Jupy_FoamQuant_T1_detection_24_0.png


G) Detect T1 events
-------------------

.. code:: ipython3

    # Read/Save image names and directories
    pairsdirname = [ProcessPipeline[4]+'/','Contact_pair_']
    readdirnamelost = [ProcessPipeline[6]+'/','Lost_']
    readdirnamenew = [ProcessPipeline[6]+'/','New_']
    
    namesave = 'T1_'
    dirsave = ProcessPipeline[7]+'/'
    
    # Images indexes
    imrange = [1,2,3,4,5]

.. code:: ipython3

    DetectT1_Batch(pairsdirname,
                   readdirnamelost,
                   readdirnamenew,
                   namesave,
                   dirsave,
                   imrange[:-1],
                   verbose=True,
                   n0=3)


.. parsed-literal::

    [1, 2, 3, 4]
    Path exist: True
    Q11_LostNewContacts/Lost_001_002
    Q11_LostNewContacts/New_001_002
    Q9_Contacts/Contact_pair_001


.. parsed-literal::

    100%|██████████| 66/66 [00:00<00:00, 365.20it/s]
    100%|██████████| 60/60 [00:00<00:00, 363.57it/s]


.. parsed-literal::

    T1_001: done
    Q11_LostNewContacts/Lost_002_003
    Q11_LostNewContacts/New_002_003
    Q9_Contacts/Contact_pair_002


.. parsed-literal::

    100%|██████████| 68/68 [00:00<00:00, 360.37it/s]
    100%|██████████| 56/56 [00:00<00:00, 359.14it/s]


.. parsed-literal::

    T1_002: done
    Q11_LostNewContacts/Lost_003_004
    Q11_LostNewContacts/New_003_004
    Q9_Contacts/Contact_pair_003


.. parsed-literal::

    100%|██████████| 83/83 [00:00<00:00, 351.87it/s]
    100%|██████████| 89/89 [00:00<00:00, 350.14it/s]


.. parsed-literal::

    T1_003: done
    Q11_LostNewContacts/Lost_004_005
    Q11_LostNewContacts/New_004_005
    Q9_Contacts/Contact_pair_004


.. parsed-literal::

    100%|██████████| 46/46 [00:00<00:00, 362.46it/s]
    100%|██████████| 56/56 [00:00<00:00, 364.62it/s]

.. parsed-literal::

    T1_004: done


.. parsed-literal::

    


.. code:: ipython3

    # Read/Save image names and directories
    namesave = 'T1_'
    dirsave = ProcessPipeline[7]+'/'
    # Read the T1
    T1_NewLost = ReadT1([dirsave,namesave],imrange[:-1], verbose=False, n0=3)

.. code:: ipython3

    # show the T1 at the first time step
    fig, ax = plt.subplots(1,3, figsize = (5*3, 5), constrained_layout=True)
    PlotT1([T1_NewLost[0][0],T1_NewLost[1][0]], ax=ax, color=None, nameaxes=['z','y','x'])



.. image:: Jupy_FoamQuant_T1_detection_files/Jupy_FoamQuant_T1_detection_29_0.png

