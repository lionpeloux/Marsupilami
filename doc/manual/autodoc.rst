.. Marsupilami documentation master file, created by
   sphinx-quickstart on Sun Oct 18 13:44:27 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*******************************
Autodoc SetUp : WIP
*******************************

Introduction
====================

Plusieurs solutions ont été testées pour réaliser une doc automatique façon "API" à partir de ``sphinx``. En effet, ``sphinx`` dispose d'un module ``autodoc`` mais celui-ci est reservé à python. Cependant, ``sphinx`` permet l'ajout de ``domain`` pour cibler d'autres langages. Pour .NET, une tentative existe mais n'est qu'en verison alpha à ce jour :

* `sphinxcontrib-dotnetdomain <https://github.com/rtfd/sphinxcontrib-dotnetdomain>`_ : un domain pour sphinx prenant en charge .Net (édité par RTD)

* `sphinx-autoapi <https://github.com/rtfd/sphinx-autoapi>`_ : une extension sphinx pour générer de la doc façon "API" automatiquement, avec prise en compte .Net (édité par RTD)

* sphinx-autoapi needs `docfx <http://vicancy.github.io/docascode/#/tutorial/docfx_getting_started.md>`_

* `breathe <https://breathe.readthedocs.org/en/latest/>`_ : une extension sphinx pour faire le lien entre doxygen et sphinx (actuellement en place, mais orienté C++, ne semble pas très bien adapté pour le C#)

* `doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ : utilitaire avec interface graphique qui permet de générer un autodoc façon "API" à partir d'un code Csharp. Fonctionne très bien. Pas d'intégration en standard avec sphinx.

* Je devrais jetter un coup d'oeil également à `MKDocs <http://www.mkdocs.org>`_

Tout celà fonctionne avec python. Pour s'affranchir de difficultées liées à de multiples installations python, un ``virtualenv`` a été mis en place au sein du dossier ``manual``.

.. note:: Un ``virtualenv`` n'est autre qu'une distribution locale et **isolée** de python dans laquelle on peut se placer pour travailler. Idéal pour les problèmes liés aux liens symbolics lorsqu'on a de multiples versions de python sur la même machine.

Environnement Virtuel (python)
==============================

Installer ``virtualenv`` pour la distribution ``python`` utilisée par défaut pour le système avec pip (la commande sudo permet d'executer la commande avec les privilèges admin si nécessaire)::

  $ [sudo] pip install virtualenv

Créer un environnement virtuel dans le dossier. ``../Marsupilami/doc/manual/_virtualenv`` avec la commande suivante::

  $ virtualenv ../Marsupilami/doc/manual/_virtualenv

Se placer dans cet environnement virtuel avec la commande suivante::

  $ source ../Marsupilami/doc/manual/_virtualenv/bin/activate

.. note:: L'invite de commande du terminal change alors de nom pour montrer que les appels pythons se feront désormais à partir de cet environnement virtuel.

Dans cet envrionnement virtuel, installer les packages requis::

  $ pip install sphinx
  $ pip install breathe
  $ pip install sphinxcontrib-dotnetdomain
  $ pip install sphinx-autoapi

Vérifier les modules installés en appelant la commande::

  $ pip list


Générer la doc
==============

La doc est générée :

* à partir des fichiers ``.rst`` du dossier ``../Marsupilami/doc/manual/``

* à partir des fichiers xml générés par ``doxygen`` et situés dans le dossier ``../Marsupilami/doc/manual/_doxygen/``

Un mémo `REST <http://rest-sphinx-memo.readthedocs.org/en/latest/index.html>`_ est disponible sur RTD.

La génération de la XML doc (et d'un HTML façon "API") avec ``doxygen`` est sans soucis.

La configuration de `breathe` se fait dans le fichier ``conf.py`` du projet sphinx::

  extensions = ['sphinx.ext.autodoc','autoapi.extension','sphinxcontrib.dotnetdomain','breathe']
  breathe_projects = { "Marsupilami": "../manual/_doxygen/xml/" }
  breathe_default_members = ('members', 'undoc-members','protected-members','private-members')
  breathe_default_project = 'Marsupilami'

Une fois placé dans l'environnement python virtuel, se placer dans le dossier ``../Marsupilami/doc/manual/`` et invoquer la commande (selon les cas)::

  $ make html
  $ make clean html

Autre point d'aide pour breathe : `tuto <https://github.com/Cruel/readthedocs-breathe/tree/master/docs>`_
