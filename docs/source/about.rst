Modules overview
=======

**foamquant** is composed of 11 modules: 

.. figure:: Diagram.png
   :alt: FoamQuant modules
   :scale: 20%

Typical analysis pipeline
-----------------
.. figure:: typicalpipeline.png
   :alt: Processing pipeline
   :width: 50%
   :align: center


Process
-----------------

The ``Process`` module focuses on the first image processing steps, aiming to produce phase-segmented and bubble-segmented images for the subsequent quantification. The provided functions are wrappers around existing tools from ``scikit-image`` [scikit-image]_, ``SPAM`` [Stamati2020]_, and ``PoreSpy`` [porespy_gostick_2019]_, enabling the processing of the images in a sequence (batch-wise). 

As illustrated in :numref:`fig-pipeline`, these include subsequent background removal, phase segmentation, masking, speckle removal, bubble segmentation, and edge bubble removal. An illustration of reconstructed, phase-segmented and bubble-segmented tomographic images is given in :numref:`fig-processing` (a--c). 

For more details, a complete example of a processing pipeline is available at `foamquant.readthedocs.io <https://foamquant.readthedocs.io>`_.

Note that although used on foams, this generic module is also applicable to analyze a wider range of porous materials, e.g., bread crumb during baking [schott2023structural]_.


.. figure:: Figure_segmentation.png
   :alt: Processing pipeline
   :width: 80%
   :align: center
   
(a) **raw reconstructed image**, (b) **phase segmented image** and (c) **bubble segmented image**.

FromBinary
-----------------

Functions to quantify the liquid fraction from a batch of phase segmented images.

.. figure:: fromliqfrac.png
   :scale: 40%
   
The liquid fraction along a cartesian mesh can be returned **structured** or **unstructured**.

FromLabelled
-----------------

Functions to quantify the bubbles regions properties from a batch of labelled images.

.. figure:: fromlab.png
   :scale: 40%
   
The function save the regions properties in a **.csv**

Tracking
-----------------

Functions to track the bubbles and their properties from a batch of labelled images.

.. figure:: tracking.png
   :scale: 40%
   
The color (from green to black) indicates the time index. The red points are the lost tracking positions.

Passage and Average
-----------------

Functions to convert scalar, vectorial or tensorial properties from cartesian to cylindrical and spherical, and perform time/space averages.

.. figure:: passage_average.png
   :scale: 40%
   
In this example the displacement field is first expressed in a cylindrical basic and then averaged.



Two ways of measuring the internal strain field
-----------------

   - Shape field, defined in [Graner2008]_ and first used in [Raufaste2015]_

   - Texture field, defined in [Graner2008]_ 
   
.. figure:: shape_texture_3d.PNG
   :scale: 50%

Label traking
-----------------
The tracking method was inspired by ID-track presented in [Ando2013]_.

.. figure:: tracking_3d.PNG
   :scale: 70%
   
Tracking of five bubbles, showing various tracked properties: elastic internal strain, number of neighbours, velocity, and volume.


References
============
.. [vanderWalt2014] S. van der Walt et al., scikit-image: Image processing in Python. PeerJ 2:e453 (2014) https://doi.org/10.7717/peerj.453

.. [stamati2020] Stamati et al., (2020). spam: Software for Practical Analysis of Materials. Journal of Open Source Software, 5(51), 2286, https://doi.org/10.21105/joss.02286

.. [Ando2013] Andò,E. et al., Experimental micromechanics: grain-scale observation of sand deformation, Géotechnique Letters 2, 107–112, (2012) https://doi.org/10.1680/geolett.12.00027

.. [Hall2010] S. A. Hall et al., Discrete and continuum analysis of localised deformation in sand using X-ray μCT and volumetric digital image correlation. Géotechnique, 60(5), 315-322, (2010) https://doi.org/10.1680/geot.2010.60.5.315

.. [Graner2008] F. Graner et al., Discrete rearranging disordered patterns, part I: Robust statistical tools in two or three dimensions, Eur. Phys. J. E 25, 349–369 (2008) https://doi.org/10.1140/epje/i2007-10298-8

.. [Raufaste2015] Raufaste, C. et al., Three-dimensional foam flow resolved by fast X-ray tomographic microscopy, EPL, 111, 38004, (2015) https://doi.org/10.1209/0295-5075/111/38004

.. [Schott2023] F. Schott et al., Three-dimensional liquid foam flow through a hopper resolved by fast X-ray microtomography, Soft Matter, (2023) https://doi.org/10.1039/d2sm01299e

.. [Schott20232] F. Schott et al., Structural formation during bread baking in a combined microwave-convective oven determined by sub-second in-situ synchrotron X-ray microtomography, Food Research International, (2023) https://doi.org/10.1016/j.foodres.2023.113283
