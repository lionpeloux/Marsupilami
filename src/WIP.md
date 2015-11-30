# WIP

Marsupilami est un solver pour le formfinding et le calcul des structures du Génie Civil basé sur des méthodes de dynamique explicite.

Il permet de faire aussi bien de la dynamique (avec schéma d'intégration explicite) que du calcul statique (par dynamique explicite amortie, telle la relaxation dynamique avec amortissement cinétique ou visqueux).

Positionnement

- rigueur scientifique
- open source
- contributif
- hebergé sur GitHub
- documenté sur Read The Docs
- intégration ciblée pour Rhino et Grasshopper

## Crédit

**Copyright** : Laboratoire Navier
**Authors** : Lionel du Peloux, Cyril Douthe, Frédérice Tayeb, Baptiste Lefevre
**Licence** : MIT

**Citation** : "Marsupilami, a dynalmic explicit solver for the formfinding and design of ligthweight structures in Civil Engineering, Laboratoire Navier, Paris, France, 2015"


## Doc

La documentation est rédigée très simplement en [ReStructuredText](https://fr.wikipedia.org/wiki/ReStructuredText) avec des inclusions [LaTeX](https://fr.wikipedia.org/wiki/LaTeX) selon les besoins. Elle est générée automatiquement avec [sphinx](http://sphinx-doc.org/). Elle est hébergée sur [Read The Docs](https://readthedocs.org/)

La documentation comprend :

- une présentation générale de l'outil
- une explication des fondements théoriques
- des exemples
- des conseils aux développeurs / contributeurs

## Dépendences

Le coeur de calcul est rendu indépendant de Rhino / Grasshopper

- [Math.NET Numerics](http://numerics.mathdotnet.com/) : vecteurs, matrices, algèbre linéaire
- [Rhinocommon.dll](http://developer.rhino3d.com/guides/rhinocommon/what_is_rhinocommon/) : OpenNurbs / intersections / distances
- [Grasshopper.dll](http://www.grasshopper3d.com/) : couche de composants grasshopper

# Architecture

Conceptuellement, le système étudié est modélisé comme une "soupe" de points matériels auxquels sont affectés des masses. L'état du système est suivi au cours du temps. A chaque instant, on connaît pour chaque noeud du système :

- `M[i]` : sa masse
- `X[i]` : sa position
- `V[i]` : sa vitesse
- `A[i]` : son accélération
- `R[i]` : la résultante des forces s'appliquant sur ce noeud
- `F[i]` : les efforts extérieurs s'appliquant sur ce noeud

**ToTHINK** : pour la torsion, il faut considérer un noeud orienté et un moment de torsion

Ces points matériels, appelés également noeuds, sont liés les uns aux autres par des lois d'interactions ou relations qui définissent entre autre :

- des **éléments**, agissant sur `M` et `R` (Beam, Cable, Chain, Bar, Tie, Strut, ...)
- des **contraintes cinématiques**, agissant sur `X` et `V`
- des **contraintes méchaniques**, agissant sur `M` et `R`
- des **liaisons**, agissant sur `X`, `V`, `M`, `R`


## model

Faut il définir l'abstraction d'un **model** ? Un modèle est un ensemble de masses ponctuelles caractérisé par :

- un état initial
- des relations entre les masses ponctuelles

## solver

un **solver** est un code qui permet de faire évoluer un modèle dans le temps, vers un état d'équilibre ou non.
element

## element

Un élément représente une interaction méchanique entre un groupe de noeuds. Cette relation donne lieu à la définition d'efforts internes qui n'ont de sens que pour l'élément lui même (`N`, `T`, `M`, `Q`). Du point de vue du système, seul la résultante des efforts en chaque noeud n'a de sens.

| Element       | Nodes |	Loop | Torsion | Flexion | Compression | Traction |
| ------------- | :---: | :--: | :-----: | :-----: | :---------: | :------: |
| Torsion Beam  | N     | yes  | x       | x       | x           | x        |
| Flexion Beam  | N     | yes  |         | x       | x           | x        |
| Cable         | N     | yes  |         |         |             | x        |
| Sliding Cable | N     | no   |         |         |             | x        |
| Bar           | 2     | no   |         |         | x           | x        |
| Tie           | 2     | no   |         |         |             | x        |
| Strut         | 2     | no   |         |         | x           |          ||

Certain éléments peuvent s'utiliser dans une configuration bouclée. D'autres types d'éléments peuvent entrer dans cette liste, comme des mailles de fillet ou des éléments de nexorade.

## constraint

On regroupe sous l'appelation constraint les autres types de relations entre noeuds.

| Constraint    |  Nodes |  M  |  X  |  V  |  R  |  F  |
| ------------- | :----: | :-: | :-: | :-: | :-: | :-: |
| Link          | 1-N    | x   | x   | x   | x   |     |
| Kinematic     | 1-N    |     | x   | x   |     |     |
| Mechanic      | 1-N    | x   |     |     | x   |     |
| Force         | 1-N    |     |     |     |     | x   ||

On peut considérer des actions exterieures qui soient dépendantes de la géométrie (forces suiveuses comme pression / amplitude de chargement fonction d'une déformée) ou du temps (sollicitations dynamiques)

- [x] sdsd

- closed beam
- presstress (N or %L)
- applied displacement
- fabric element
