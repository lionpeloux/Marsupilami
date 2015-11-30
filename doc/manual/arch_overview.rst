.. Marsupilami documentation master file, created by
   sphinx-quickstart on Sun Oct 18 13:44:27 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview
========

Intro
-----

Conceptuellement, le système étudié est modélisé comme une "soupe" de points matériels auxquels sont affectés des masses. L'état du système est suivi au cours du temps. A chaque instant, on connaît pour chaque noeud du système :

- M[i] : sa masse
- X[i] : sa position
- V[i] : sa vitesse
- A[i] : son accélération
- R[i] : la résultante des forces s'appliquant sur ce noeud
- F[i] : les efforts extérieurs s'appliquant sur ce noeud

**ToTHINK** : pour la torsion, il faut considérer un noeud orienté et un moment de torsion

Ces points matériels, appelés également noeuds, sont liés les uns aux autres par des lois d'interactions ou relations qui définissent entre autre :

- des **éléments**, agissant sur `M` et `R` (Beam, Cable, Chain, Bar, Tie, Strut, ...)
- des **contraintes cinématiques**, agissant sur `X` et `V`
- des **contraintes méchaniques**, agissant sur `M` et `R`
- des **liaisons**, agissant sur `X`, `V`, `M`, `R`


Modèle
-----------

Faut il définir l'abstraction d'un **model** ? Un modèle est un ensemble de masses ponctuelles caractérisé par :

- un état initial
- des relations entre les masses ponctuelles

Solver
-----------

Un **solver** est un code qui permet de faire évoluer un modèle dans le temps, vers un état d'équilibre ou non.
element

Elément
-------

Un élément représente une interaction méchanique entre un groupe de noeuds. Cette relation donne lieu à la définition d'efforts internes qui n'ont de sens que pour l'élément lui même (`N`, `T`, `M`, `Q`). Du point de vue du système, seul la résultante des efforts en chaque noeud n'a de sens.

.. table:: Listing of Element Types

+---------------+------+------+---------+---------+-------------+----------+
| Element       | Node | Loop | Torsion | Flexion | Compression | Traction |
+===============+======+======+=========+=========+=============+==========+
| Torsion Beam  | N    | yes  | x       | x       | x           | x        |
+---------------+------+------+---------+---------+-------------+----------+
| Flexion Beam  | N    | yes  |         | x       | x           | x        |
+---------------+------+------+---------+---------+-------------+----------+
| Cable         | N    | yes  |         |         |             | x        |
+---------------+------+------+---------+---------+-------------+----------+
| Sliding Cable | N    | no   |         |         |             | x        |
+---------------+------+------+---------+---------+-------------+----------+
| Bar           | 2    | no   |         |         | x           | x        |
+---------------+------+------+---------+---------+-------------+----------+
| Tie           | 2    | no   |         |         |             | x        |
+---------------+------+------+---------+---------+-------------+----------+
| Strut         | 2    | no   |         |         | x           |          |
+---------------+------+------+---------+---------+-------------+----------+
Certain éléments peuvent s'utiliser dans une configuration bouclée. D'autres types d'éléments peuvent entrer dans cette liste, comme des mailles de fillet ou des éléments de nexorade.


Contrainte
----------

On regroupe sous l'appelation constraint les autres types de relations entre noeuds.

.. table:: Listing of Constraint Types

+------------+-------+---+---+---+---+---+
| Constraint | Nodes | M | X | V | R | F |
+============+=======+===+===+===+===+===+
| Link       | 1-N   | X | X | X | X |   |
+------------+-------+---+---+---+---+---+
| Kinematic  | 1-N   |   | X | X |   |   |
+------------+-------+---+---+---+---+---+
| Mechanic   | 1-N   | X |   |   | X |   |
+------------+-------+---+---+---+---+---+
| Force      | 1-N   |   |   |   |   | X |
+------------+-------+---+---+---+---+---+

On peut considérer des actions exterieures qui soient dépendantes de la géométrie (forces suiveuses comme pression / amplitude de chargement fonction d'une déformée) ou du temps (sollicitations dynamiques)
