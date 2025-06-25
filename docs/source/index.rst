Welcome to foamquant documentation
===================================

The ``pip`` installable Python library, **foamquant**, was primarily developed to study time-resolved X-ray tomograms of evolving liquid foams and their fundamental local physical behavior. However, we believe this library can be adapted to the broader class of cellular materials, using other imaging methods, and put to work on other scientific and engineering applications.

The toolbox has four main dependencies: `numpy <https://numpy.org/doc/stable/>`_, `scikit-image <https://scikit-image.org/>`_ [scikit-image]_, `SPAM <https://ttk.gricad-pages.univ-grenoble-alpes.fr/spam/>`_ [Stamati2020]_, and `porespy <https://porespy.org/>`_ [porespy_gostick_2019]_.

Main objectives:

- Integrating standard image processing methods from ``scikit-image``, ``SPAM``, and ``porespy``.
- Including new mechanical quantification methods specifically dedicated to studying evolving liquid foam structure: local elasticity and plasticity.

This website provides an **API reference** and ready-to-use **jupyter notebooks** image analysis pipelines.

.. code::
   # Requirements
   numpy==1.17.2
   matplotlib==3.3.4
   matplotlib-inline==0.1.6
   scikit-image==0.18.3
   tifffile==2021.11.2
   pandas==0.25.1
   spam==0.6.0.3

.. toctree::
   :maxdepth: 2

   about
   api
   examples

   
References
-----------------------------------

.. [scikit-image] van der Walt, S., et al. (2014) scikit-image: image processing in Python. PeerJ 2, e453.

.. [Stamati2020] Stamati, H., et al. (2020). SPAM: Software for Practical Analysis of Materials, 5, 2286.

.. [porespy_gostick_2019] Gostick, J. et al. (2019). PoreSpy: A Python toolkit for quantitative analysis of porous media images. Journal of Open Source Software, 4(37), 1296.

.. note::

   FoamQuant requires Python 3.8 or higher.

   - Make sure to install dependencies: ``numpy``, ``scikit-image``, ``SPAM``, and ``porespy``.
   - Consider using a virtual environment to avoid conflicts.
   - For more information or questions regarding this project please contact: `Florian Schott <florian.schott@solid.lth.se>`_ or `Rajmund Mokso <rajmo@dtu.dk>`_.
