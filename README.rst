Welcome to foamquant's documentation!
===================================

**foamquant** is a toolbox specifically created for processing time series of 3D images of evolving liquid foams. The toolbox is greatly based on open source libraries: Scikit-Image [vanderWalt2014]_ and SPAM [Stamani2020]_. 

The **objective** is a greater accessibility of time-resolved liquid foam 3D analysis tools. The tools focus mainly on flow studies, but can also be used for coarsening and drainage studies.

We propose three sections: **Jupyter Notebooks**, **Elastic and plastic deformations quantification tools** and **Label traking**.

Jupyter Notebooks
-----------------
How we can use Scikit-Image and SPAM tools for extracting:

1) Liquid fraction and Plateau border size

2) Bubble volume distribution

3) Coordination distribution (uing SPAM [Stamani2020]_)

4) Foam flow field (inspired by ID-track in [Ando2013]_ or DVC SPAM-ddic [Hall2010]_ [Ando2013]_)


Elastic and plastic deformations quantification tools
-----------------

1) Elastic deformation field:

   - Shape field [Graner2008]_ [Raufaste2015]_

   - Texture field [Graner2008]_ 

   - Elastic internal strain field [Raufaste2015]_

2) Detection of created / disappeared films (contacts)

3) Plastic deformation field (T1)

Label traking 
-----------------
The tracking tool was inspired by ID-track presented in [Ando2013]_ and uses DVC SPAM-ddic [Hall2010]_ [Ando2013]_.

1) Bubble traking

2) Film traking


References
============
.. [vanderWalt2014] S. van der Walt et al., scikit-image: Image processing in Python. PeerJ 2:e453 (2014) https://doi.org/10.7717/peerj.453 

.. [Stamani2020] Stamati et al., (2020). spam: Software for Practical Analysis of Materials. Journal of Open Source Software, 5(51), 2286, https://doi.org/10.21105/joss.02286

.. [Ando2013] Andò,E. et al., Experimental micromechanics: grain-scale observation of sand deformation, Géotechnique Letters 2, 107–112, (2012) https://doi.org/10.1680/geolett.12.00027

.. [Hall2010] S. A. Hall et al., Discrete and continuum analysis of localised deformation in sand using X-ray μCT and volumetric digital image correlation. Géotechnique, 60(5), 315-322, (2010) https://doi.org/10.1680/geot.2010.60.5.315

.. [Graner2008] F. Graner et al., Discrete rearranging disordered patterns, part I: Robust statistical tools in two or three dimensions, Eur. Phys. J. E 25, 349–369 (2008) https://doi.org/10.1140/epje/i2007-10298-8

.. [Raufaste2015] Raufaste, C. et al., Three-dimensional foam flow resolved by fast X-ray tomographic microscopy, EPL, 111, 38004, (2015) https://doi.org/10.1209/0295-5075/111/38004




.. note::

   This project is under development. The Jupyter notebook are not uploaded yet, neither is the tool library.

