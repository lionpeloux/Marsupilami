Welcome to Marsupilami's documentation !
=======================================

Marsupilami est un solver pour le formfinding et le calcul des structures du Génie Civil basé sur des méthodes de dynamique explicite.

.. image:: img/marsupilami.png
    :width: 125px
    :align: center
    :alt: marsupilami

Marsupilami permet de faire aussi bien de la dynamique (avec schéma d'intégration explicite) que du calcul statique (par dynamique explicite amortie, telle la relaxation dynamique avec amortissement cinétique ou visqueux).

Positionnement :

 - rigueur scientifique
 - open source
 - contributif
 - hebergé sur GitHub
 - documenté sur Read The Docs
 - intégration ciblée pour Rhino et Grasshopper


The main documentation for the site is organized into a couple sections:


* :ref:`arch-doc`
* :ref:`user-doc`
* :ref:`developer-doc`

.. _arch-doc:

.. toctree::
   :maxdepth: 1
   :caption: Architecture

   arch_overview.rst
   arch_element.rst
   arch_constraint.rst

.. _user-doc:

.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   usr_overview.rst
   usr_tutorial.rst
   usr_theory.rst

.. _developer-doc:

.. toctree::
    :maxdepth: 2
    :caption: Developer Documentation

    dev_overview.rst
    dev_htmldoc.rst
    dev_autoapi.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
