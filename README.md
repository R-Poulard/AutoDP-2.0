# AutoDP - Stage M2 AMIIB à l'Université Paris-Saclay

Ce dépôt Git correspond au travail effectué lors de mon stage de deuxième semestre de M2 AMIIB au sein de l'équipe AMIBIO du LIX (Laboratoire d'informatique de l'École polytechnique).

## Objectif du Stage

L'objectif principal de ce stage était de finaliser l'implémentation d'AutoDP, un pipeline d'automatisation basé sur les travaux de B. Marchand et al. ([Publication](https://almob.biomedcentral.com/articles/10.1186/s13015-023-00229-z#Sec3)). L'amélioration porte sur la génération du code C ainsi que du fichier binaire correspondant, en sortie du pipeline. Cette version d'AutoDP permet l'implémentation d'un modèle de Turner standard avec des pseudonœuds additionnels, en utilisant le ViennaRNA Package ([ViennaRNA](https://www.tbi.univie.ac.at/RNA/)) et une fonctionnalité unique d'ajout de grammaire auxiliaire.
La version modifiée du ViennaRNA Package, permettant l'ajout de la grammaire auxiliaire dans le schéma d'équation de programmation dynamique n'est pas encore disponibles sur le site officiel de ViennaRNA. Merci de les citer et de prendre contact avec les développeurs pour toute utilisation de ce paquet ou demande de cette version encore en developpement.

## Structure du Dépôt

- **`auto-dp/ressources/stack-Turner/`**  
  Contient les modèles utilisés pour la génération du code.

- **`auto-dp/workflow/script/produce_code2_Stack_Turner.py`**  
  Script Python permettant de générer les fichiers `.C` nécessaires.

- **`auto-dp/workflow/script/produce_code2_Stack_Turner.py`**  
  Un autre script Python pour la génération de code.

- **`test_sequence/`**  
  Ce répertoire contient les séquences de test pour le modèle de "base-pairing".

- **`sequence_H.txt` et `sequence_shuffled.txt`**  
  Fichiers contenant les séquences utilisées pour tester le pouvoir discriminant du programme.

- **`simple_rfam_data/`**  
  Répertoire pour récupérer les séquences RFAM, en grande partie inspiré de ce [repository](https://github.com/bmarchand/simple_rfam_structure_extraction).

- **`test_sequence/`**  
  Ce répertoire contient les séquences de test pour le modèle de "base-pairing".

- **`Test/`**
  Ce répertoire contient les script nécessaire pour tester les séquences comme décrit dans le rapport
  
## Informations Supplémentaires

Pour plus d'informations sur le pipeline AutoDP et son utilisation, veuillez consulter le [repository GitHub originel de B. Marchand](https://github.com/bmarchand/auto-dp/tree/main).
