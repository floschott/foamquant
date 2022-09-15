Usage
=====

.. _installation:

Installation
------------

To use foamquant, first install it using pip:

.. code-block:: console

   (.venv) $ pip install foamquant

Creating recipes
----------------

Use the ``FoamQuant.Figure.Cut3D()`` function:

.. autofunction:: FoamQuant.Figure.Histogram

The ``blabla`` parameter should be either ``"blabla"``, ``"blabla"``,
or ``"blabla"``. Otherwise, :py:func:`FoamQuant.Figure.Histogram`
will raise an exception.

.. autoexception:: FoamQuant.Figure.Histogram

For example:

>>> import foamquant
>>> foamquant.Figure.Histogram()


.. toctree::
   
   index
   usage
   api
