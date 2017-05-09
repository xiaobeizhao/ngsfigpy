
============================
NGSfig
============================


Example: Legoplot in NGSfig
----------------------------------

.. code-block:: python

  import warnings
  from pkg_resources import get_distribution, resource_stream, Requirement
  from ngsfig.graphics import trinucleotide

  if get_distribution("setuptools").version < "33.1.1":
    warnings.warn("The code below has only been tested via setuptools (=33.1.1)")
  
  with resource_stream(Requirement.parse("ngsfig"), "data/trinucleotide_demo.txt") as input:
    trinucleotide.demo(input)

Run this script or paste it into a Python console.


.. class:: no-web

  .. image:: https://raw.githubusercontent.com/xiaobeizhao/ngsfigpy/master/data/trinucleotide_demo.png
    :alt: NGSfig | Legoplot
    :width: 100%
    :align: center

.. class:: no-web no-pdf

           
Download and Install
--------------------
* `GitHub <https://github.com/xiaobeizhao/xmiscpy>`
* `GitHub <https://github.com/xiaobeizhao/ngsfigpy>`

  
License
-------
Code and documentation are available according to the GNU LGPL License.


