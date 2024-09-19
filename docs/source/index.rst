Welcome to foamquant documentation
===================================

.. important::
   This project is under construction. New functionalities are constantly added to this package.
   The project is directly related to the study [Schott2023]_.

=======

**foamquant** is a toolbox specifically created for processing time series of 3D images of evolving liquid foams by using open source libraries Scikit-Image [vanderWalt2014]_ and SPAM [Stamani2020]_. 

The **objective** is a greater accessibility of time-resolved liquid foam 3D analysis tools for the foam-physicis and food scientists communitites. 
We propose the following documentation: **API** and **Jupyter Notebooks**.

Installation: neither ``conda``, nor ``pip`` installations are available yet.

Main dependencies: `SPAM <https://ttk.gricad-pages.univ-grenoble-alpes.fr/spam/>`_ and  `scikit-image <https://scikit-image.org/>`_

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

   
.. note::
   Contacts: `Florian Schott <florian.schott@solid.lth.se>`_ and `Rajmund Mokso <rajmo@dtu.dk>`_
