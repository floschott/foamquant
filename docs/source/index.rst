Welcome to foamquant's documentation! (beginning of it)
===================================

**foamquant** FoamQuant is a toolbox specifically created for processing time series of large 3D images of evolving liquid foams. The toolbox is greatly based on open source libraries: scikit-image [1] and SPAM [2]. 

The objective is a greater accessibility of time-resolved liquid foam 3D analysis tools. It focus mainly on flow studies, but can also be used in a certain extend for coarsening and drainage studies. 


"Basic" tools:

1) Liquid fraction

2) Volume and dispersity

3) Coordination (uing SPAM [2])


"Advanced" tools:

1) Flow field (using ID-track [3] or DVC SPAM-ddic [4])

2) Elastic deformation field

- Shape field [5][6]

- Texture field [5]

- Elastic internal strain field

5) Lost/New contact detection (uing SPAM [2])

6) Plastic deformation field (topological transitions, uing SPAM [2])


Tracking tools:

1) Bubble traking (using ID-track [3])

2) Films traking (using ID-track [3])


References:

[1] Van der Walt et al, doi: https://doi.org/10.7717/peerj.453

[2] Stamati et al, doi: https://doi.org/10.21105/joss.02286

[3] Stamati et al, doi: https://doi.org/10.21105/joss.02286

[4] Stamati et al, doi: https://doi.org/10.21105/joss.02286

[5] Stamati et al, doi: https://doi.org/10.21105/joss.02286

[6] Stamati et al, doi: https://doi.org/10.21105/joss.02286


Check the :doc:`functions` section for further information.

.. note::

   This project is under development.

Contents
--------

.. toctree::

   usage
   api
   process
   quantif
   figure
