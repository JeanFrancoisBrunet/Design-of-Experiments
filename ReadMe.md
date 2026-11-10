# DOE — Plans d'expériences

Application de bureau (Python / Tkinter) pour **construire et analyser des plans d'expériences**. L'interface est en Python ; tous les calculs statistiques et les graphiques sont réalisés par **R** (R de base + `jsonlite`, aucun autre package).

## Fonctions

| Onglet | Famille de plans | Contenu |
|---|---|---|
| 1. Criblage | Plackett-Burman (8 à 32 essais) ou fractions 2^(k−p) | Repérer les facteurs influents avec un minimum d'essais |
| 2. Factoriel | Plans complets et fractionnaires à 2 niveaux | Générateurs, résolution, structure d'alias, effets et interactions, réduction du modèle, prévision, optimisation |
| 3. Taguchi | Tableaux orthogonaux L4 à L32 (2, 3 ou 5 niveaux) | Interactions, niveau fictif, rapport signal/bruit (4 types), ANOVA, contributions, niveaux optimaux, intervalle de confirmation, prévision |
| 4. Surface | Composite centré (circonscrit, à faces centrées, inscrit), Box-Behnken | Modèle du second degré, point stationnaire, lignes de niveau, surface 3D, optimisation |
| 5. Mélange | Réseau simplexe, centroïde, sommets extrêmes | Modèles de Scheffé (linéaire à cubique), pseudo-composants, diagramme ternaire, optimisation |

Chaque famille suit les mêmes sous-onglets :

1. **Plan** : définition des facteurs, options, génération de la matrice d'essais (ordre aléatoire reproductible), saisie des réponses dans le tableau, collage depuis un tableur, export CSV.
2. **Analyse** : ANOVA, effets, coefficients, équations (unités codées et réelles), R² / R² ajusté / R² prévu, observations inhabituelles, normalité des résidus, graphiques (Pareto, Daniel, effets principaux, interactions, cube, résidus).
3. **Graphiques du modèle** (factoriel, surface, mélange).
4. **Prévision et optimisation** : valeurs prévues avec intervalles de confiance et de prévision ; optimisation multi-réponses par fonction de désirabilité (Derringer et Suich).

L'onglet **Rapport** rassemble les résultats calculés : on choisit ceux à inclure et on exporte une page HTML autonome (plans, textes et graphiques). Le **projet** complet (plans, réponses, résultats) s'enregistre dans un fichier `.doeproj` que l'on peut rouvrir.

## Installation

```bash
sudo apt install -y r-base r-cran-jsonlite python3-tk python3-pil.imagetk
```

`python3-pil.imagetk` est facultatif (affichage optimisé des graphiques).

## Lancement

```bash
cd ~/Projects/DOE
python3 doe.py
```

## Fichiers

| Fichier | Rôle |
|---|---|
| `doe.py` | Interface Tkinter, gestion des plans, du rapport et du projet |
| `cli_doe.R` | Passerelle : lit `params.json` et `plan.csv`, écrit `resultat.json` et les graphiques PNG |
| `plans_doe.R` | Construction des plans (criblage, factoriel, Taguchi, surface, mélange) |
| `moteur_doe.R` | Moteur de modélisation, prévision, optimisation |
| `taguchi_doe.R` | Analyse des plans de Taguchi |
| `graphiques_doe.R` | Tous les graphiques |
| `telegram_bot_doe.py` | Bot Telegram (voir ci-dessous) |
| `icons/DOE.png` | Icône |

Les cinq fichiers R doivent rester dans le même dossier que `doe.py`.

## Utilisation

- **Nombres** : la virgule et le point sont acceptés comme séparateur décimal.
- **Facteurs** : niveaux numériques (bas / haut) ou, pour les plans de criblage et factoriels, deux niveaux qualitatifs. Dans un plan de Taguchi, les niveaux d'un facteur sont séparés par des points-virgules. Dans un plan de mélange, on saisit les bornes de chaque composant (entre 0 et 1) ; si des bornes supérieures contraignent le domaine, choisir le plan à sommets extrêmes.
- **Réponses** : double-clic sur une cellule marquée ✎ du tableau, ou « Saisir ou coller des réponses… » (une colonne par réponse, ordre d'exécution ou ordre standard). Plusieurs réponses peuvent être ajoutées, renommées ou supprimées.
- **Risque α** : réglé en haut de la fenêtre, il s'applique aux tests et aux intervalles.
- **Réduction du modèle** : élimination descendante des termes non significatifs, en respectant la hiérarchie des termes.
- **Plans sans degré de liberté pour l'erreur** : les effets sont alors évalués par la méthode de Lenth, et les intervalles ne sont pas calculables ; le programme le signale.
- **Optimisation** : pour chaque réponse, choisir le but (maximiser, minimiser, atteindre une cible) et, si besoin, les limites et l'importance. Le bouton « Utiliser l'optimum » reporte les réglages trouvés dans l'onglet de prévision.

## Bot Telegram

`telegram_bot_doe.py` donne accès aux mêmes plans et analyses depuis Telegram. Il importe `doe.py` (même moteur R) et doit rester dans le même dossier.

**Configuration** : le jeton se lit dans la variable `TELEGRAM_TOKEN`, sinon dans `~/.telegram_config`, section `[telegram]`, clé `token_doe`. Les identifiants autorisés viennent de `ALLOWED_IDS` ou de la clé `allowed_ids` ; le `chat_id` du fichier est toujours autorisé et reçoit le message de démarrage. Un utilisateur refusé reçoit son identifiant Telegram.

**Dépendance** : `requests` (`sudo apt install -y python3-requests`).

**Lancement** : `python3 telegram_bot_doe.py`

**Démarche**
1. Créer un plan : `/criblage`, `/factoriel`, `/taguchi`, `/surface` ou `/melange`. La commande affiche la syntaxe : un facteur par ligne (`Nom : niveau_bas niveau_haut`), puis des options (`type=…`, `rep=…`, `graine=…`).
2. `/reponses` puis envoyer les mesures (une ligne par essai, une colonne par réponse ; un fichier csv ou txt est aussi accepté).
3. `/analyse`, `/graphiques`, `/prevoir`, `/optimiser`. Les options du modèle utilisées sont reprises d'une commande à l'autre.
4. `/famille` change de plan actif, `/plan` l'affiche, `/export` le renvoie en CSV, `/alpha` règle le risque.

Les plans restent en mémoire du bot, le temps de la conversation.

## Principes de calcul

- Facteurs codés −1 / +1 ; pseudo-composants pour les mélanges.
- Sélection séquentielle des colonnes indépendantes, sommes des carrés ajustées, terme de courbure en présence de points centraux.
- Fractions définies par générateurs ; résolution et alias affichés.
- Taguchi : tableaux construits sur les corps de Galois, rapport signal/bruit, regroupement automatique des termes non significatifs dans l'erreur, intervalle de confirmation.
- Optimisation : échantillonnage aléatoire puis affinage local (L-BFGS-B) ; pour les mélanges, combinaisons convexes des sommets puis transferts entre composants.

## Dépannage

- « R est introuvable » ou « jsonlite absent » : relancer la commande d'installation ci-dessus.
- Les graphiques n'apparaissent pas : vérifier que R dispose du périphérique `cairo` (inclus dans `r-base` standard).
