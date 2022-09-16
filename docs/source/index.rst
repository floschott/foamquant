Welcome to foamquant's documentation! (project in progress)
===================================

**foamquant** is a toolbox specifically created for processing time series of 3D images of evolving liquid foams by using open source libraries Scikit-Image [vanderWalt2014]_ and SPAM [Stamani2020]_. 

The **objective** is a greater accessibility of time-resolved liquid foam 3D analysis tools for the foam-physicis and food scientists communitites.

We propose a readthedocs documentation with three sections: **API reference**, **Jupyter Notebooks examples** and  **Tutorials**


Tools
-----------------

.. figure:: Diagram.PNG
   :scale: 50%
   
Diagram showing the organisation of the python functions in the different sub-modules. The functions in red are not included in the package yet.



Elastic and plastic fields quantification tools
-----------------

1) Elastic deformation field:

   - Shape field, defined in [Graner2008]_ and first used in [Raufaste2015]_

   - Texture field, defined in [Graner2008]_ 
   
.. figure:: shape_texture_3d.PNG
   :scale: 50%

Shape and texture in 3D

   - Elastic internal strain field, defined in [Graner2008]_ 

2) Detection of created / disappeared films (contacts)

3) Plastic deformation field (T1)

Label traking 
-----------------
The tracking method was inspired by ID-track presented in [Ando2013]_ and uses DVC SPAM-ddic [Hall2010]_ [Ando2013]_.

1) Bubble traking

.. figure:: tracking_3d.PNG


Tracking of five bubbles, showing various tracked properties: elastic internal strain, contact, velocity, and volume.

2) Film traking


Jupyter Notebooks
-----------------
Example of extracted properties:

1) Liquid fraction and Plateau border radius

2) Bubble volume distribution and Volume image

3) Coordination distribution and Coordination image

4) Foam flow field (thanks to ID-track in [Ando2013]_ and DVC SPAM-ddic [Hall2010]_ [Ando2013]_)

5) Bubble tracking and time resolved visualization in ParaView


References
============
.. [vanderWalt2014] S. van der Walt et al., scikit-image: Image processing in Python. PeerJ 2:e453 (2014) https://doi.org/10.7717/peerj.453 

.. [Stamani2020] Stamati et al., (2020). spam: Software for Practical Analysis of Materials. Journal of Open Source Software, 5(51), 2286, https://doi.org/10.21105/joss.02286

.. [Ando2013] Andò,E. et al., Experimental micromechanics: grain-scale observation of sand deformation, Géotechnique Letters 2, 107–112, (2012) https://doi.org/10.1680/geolett.12.00027

.. [Hall2010] S. A. Hall et al., Discrete and continuum analysis of localised deformation in sand using X-ray μCT and volumetric digital image correlation. Géotechnique, 60(5), 315-322, (2010) https://doi.org/10.1680/geot.2010.60.5.315

.. [Graner2008] F. Graner et al., Discrete rearranging disordered patterns, part I: Robust statistical tools in two or three dimensions, Eur. Phys. J. E 25, 349–369 (2008) https://doi.org/10.1140/epje/i2007-10298-8

.. [Raufaste2015] Raufaste, C. et al., Three-dimensional foam flow resolved by fast X-ray tomographic microscopy, EPL, 111, 38004, (2015) https://doi.org/10.1209/0295-5075/111/38004




.. note::

   This project is under development. The Jupyter notebook are not uploaded yet, neither are all the tools.


.. toctree::
   :glob:
