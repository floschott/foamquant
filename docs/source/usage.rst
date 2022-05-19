Usage
=====

.. _installation:

Installation
------------

To use foamquant, first install it using pip:

.. code-block:: console

   (.venv) $ pip install foamquant


Using foamquant
----------------

To create a random color map that can be loaded in Paraview or Tomviz,
you can use the ``foamquant.json_rand_dictionary()`` function:

.. autofunction:: foamquant.json_rand_dictionary

The ``Ncolors`` parameter should be ``"int"``. Otherwise, :py:func:`foamquant.json_rand_dictionary`
will raise an error.

For example:

>>> import foamquant
>>> foamquant.json_rand_dictionary(5000, 'Random_ColorMap_5000colors',first_color_black=True)


