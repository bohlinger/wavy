Importing multiple satellite products
#####################################

The multisat module is a shortcut to import multiple satellite products at once. It uses the consolidate module [:ref:`consolidate-label`] under the hood. In short it is used like:

.. code-block:: python3

    >>> from wavy.multisat_module import multisat_class as ms
    >>> mso = ms(sd='2022-02-01', ed='2022-02-03', name=["s3a", "s3b"])
