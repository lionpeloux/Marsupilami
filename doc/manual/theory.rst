.. Marsupilami documentation master file, created by
   sphinx-quickstart on Sun Oct 18 13:44:27 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Theory
======

Marsupilami is a non linear solver for networks of beams, bars and cables.
It is based on a dynamic relaxation algorithm.


Overview
--------

Darboux vector
Those equations can be formulated with the \emph{Darboux vector} of the chosen material frame,
which represents the rotational velocity of the frame along :math:`\boldsymbol{x}(s)` :

Darboux Vector
--------------

.. math::
  \boldsymbol{d}'_{i}(s) = \boldsymbol{\Omega}_m(s) \times \boldsymbol{d}_i(s)
	\quad,\quad
	\boldsymbol{\Omega}_m(s)
	=
	\begin{bmatrix}
		\tau(s) \\
		\kappa_{1}(s) \\
		\kappa_{2}(s)
	\end{bmatrix}



Dynamic relaxation
------------------



Elements
----------
