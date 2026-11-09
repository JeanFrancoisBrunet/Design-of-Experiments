#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
#  DOE — Plans d'expériences — Raspberry Pi 5
#  Application de bureau (Tkinter). Les calculs sont délégués à R (cli_doe.R).
#   1. Plan de criblage (Plackett-Burman, fractions) : repérer les facteurs influents
#   2. Plans factoriels complets et fractionnaires : matrice, analyse, graphiques, interactions,
#      optimisation du modèle, prévision
#   3. Plans de Taguchi (tableaux orthogonaux, rapport signal/bruit) : mêmes fonctions
#   4. Plans de surface de réponse (composite centré, Box-Behnken)
#   5. Plans de mélange (réseau simplexe, centroïde, sommets extrêmes)
#
#  Dépendances : Python 3 + Tkinter (python3-tk) ; R + jsonlite (r-base, r-cran-jsonlite)
#                Pillow (facultatif, affichage des graphiques)
#  Lancement   : python3 doe.py
#  Auteur : Jean-François BRUNET – JFBConseils – Novembre 2026
#  Rev.   :
# =============================================================================

import atexit
import base64
import csv
import datetime
import html
import json
import math
import os
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import zipfile
from pathlib import Path

DOSSIER = Path(__file__).resolve().parent
ICONE = DOSSIER / "icons" / "DOE.png"
SCRIPT_R = DOSSIER / "cli_doe.R"
FICHIERS_R = ("cli_doe.R", "plans_doe.R", "moteur_doe.R", "graphiques_doe.R", "taguchi_doe.R")
TITRE = "DOE — Plans d'expériences"
EXT_PROJET = ".doeproj"
FORMAT_PROJET = "DOE"
VERSION_PROJET = 1
MAX_ESSAIS = 5000
SYSTEME = ("Ordre_essai", "Ordre_standard", "Type_point")      # colonnes gérées par le programme
RESERVES = {c.lower() for c in SYSTEME}


class ErreurDOE(Exception):
    """Erreur destinée à l'utilisateur (message lisible)."""


# =============================================================================
#  Nombres et mise en forme
# =============================================================================
def lire_nombre(texte, nom="valeur"):
    """Nombre saisi (virgule ou point décimal) ; chaîne vide -> None."""
    if texte is None:
        return None
    t = str(texte).strip().replace(" ", "").replace(" ", "")
    if t == "":
        return None
    try:
        v = float(t.replace(",", "."))
    except ValueError:
        raise ErreurDOE(f"Valeur numérique invalide pour {nom} : « {texte} ».")
    if math.isnan(v) or math.isinf(v):
        raise ErreurDOE(f"Valeur numérique invalide pour {nom} : « {texte} ».")
    return v


def est_nombre(texte):
    try:
        return lire_nombre(texte) is not None
    except ErreurDOE:
        return False


def lire_alpha(texte):
    t = str(texte or "").strip()
    pct = t.endswith("%")
    v = lire_nombre(t.rstrip("%").strip(), "le risque α")
    if v is None:
        raise ErreurDOE("Saisir le risque α.")
    if pct or v >= 1:
        v /= 100
    if not (0 < v < 0.5):
        raise ErreurDOE("Le risque α doit être compris entre 0 et 0,5.")
    return v


def fmt_nombre(v, chiffres=5):
    if v is None:
        return "—"
    if isinstance(v, bool):
        return "oui" if v else "non"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if math.isnan(v):
            return "—"
        if v == int(v) and abs(v) < 1e9:
            return str(int(v))
        return f"{v:.{chiffres}g}".replace(".", ",")
    return str(v)


LIBELLES = {"p_value": "p-value", "Contribution_pct": "Contribution (%)", "Contribution_pure_pct": "Contribution pure (%)",
            "IC_inf": "IC inf.", "IC_sup": "IC sup.", "Valeur_codée": "Valeur codée", "Valeur_réelle": "Valeur réelle",
            "RSB": "RSB (dB)", "t_Lenth": "t de Lenth", "Somme_carrés": "Somme des carrés", "Carré_moyen": "Carré moyen",
            "Erreur_type": "Erreur type", "Ecart_type": "Écart-type", "Écart_type": "Écart-type"}


def libelle(c):
    return LIBELLES.get(c, str(c).replace("_", " "))


def fmt_cellule(nom, v, libre=False):
    if (not libre and isinstance(v, (int, float)) and not isinstance(v, bool) and str(nom).lower().startswith("p_")):
        return "< 0,0001" if v < 1e-4 else f"{v:.4f}".replace(".", ",")
    return fmt_nombre(v)


def table_texte(cols, rows, libre=False):
    cellules = [[fmt_cellule(c, v, libre) for c, v in zip(cols, r)] for r in rows]
    ent = [str(c) if libre else libelle(c) for c in cols]
    larg = [max(len(ent[i]), *(len(r[i]) for r in cellules)) if cellules else len(ent[i]) for i in range(len(cols))]

    def ligne(vals):
        return "  ".join(v.ljust(larg[i]) if i == 0 else v.rjust(larg[i]) for i, v in enumerate(vals))
    return "\n".join([ligne(ent), "  ".join("-" * w for w in larg)] + [ligne(r) for r in cellules])


def _liste(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return list(v)
    if isinstance(v, dict):
        return list(v.values())
    return [v]


# =============================================================================
#  Fractions factorielles : générateurs et résolution
# =============================================================================
# {nombre de facteurs: [générateurs, ...]} — la fraction de résolution la plus élevée pour chaque nombre d'essais
FRACTIONS = {
    3: [["C=AB"]],
    4: [["D=ABC"]],
    5: [["E=ABCD"], ["D=AB", "E=AC"]],
    6: [["F=ABCDE"], ["E=ABC", "F=BCD"], ["D=AB", "E=AC", "F=BC"]],
    7: [["G=ABCDEF"], ["F=ABCD", "G=ABDE"], ["E=ABC", "F=BCD", "G=ACD"], ["D=AB", "E=AC", "F=BC", "G=ABC"]],
    8: [["G=ABCD", "H=ABEF"], ["F=ABC", "G=ABD", "H=BCDE"], ["E=BCD", "F=ACD", "G=ABC", "H=ABD"]],
    9: [["G=ABCD", "H=ACEF", "I=CDEF"], ["F=BCDE", "G=ACDE", "H=ABDE", "I=ABCE"], ["E=ABC", "F=BCD", "G=ACD", "H=ABD", "I=ABCD"]],
    10: [["G=BCDF", "H=ACDF", "I=ABDE", "J=ABCE"], ["F=ABCD", "G=ABCE", "H=ABDE", "I=ACDE", "J=BCDE"],
         ["E=ABC", "F=BCD", "G=ACD", "H=ABD", "I=ABCD", "J=AB"]],
    11: [["G=CDE", "H=ABCD", "I=ABF", "J=BDEF", "K=ADEF"], ["F=ABC", "G=BCD", "H=CDE", "I=ACD", "J=ADE", "K=BDE"],
         ["E=ABC", "F=BCD", "G=ACD", "H=ABD", "I=ABCD", "J=AB", "K=AC"]],
}
ROMAIN = {3: "III", 4: "IV", 5: "V", 6: "VI", 7: "VII", 8: "VIII"}


def resolution_fraction(k, gens):
    """Résolution d'une fraction 2^(k-p) d'après ses générateurs (None si invalide)."""
    try:
        mots = []
        for g in gens:
            cible, rhs = g.replace(" ", "").upper().split("=")
            rhs = rhs.lstrip("+-")
            m = 1 << (ord(cible) - 65)
            for ch in rhs:
                m ^= 1 << (ord(ch) - 65)
            mots.append(m)
        mini = 99
        for s in range(1, 2 ** len(mots)):
            m = 0
            for j, w in enumerate(mots):
                if s >> j & 1:
                    m ^= w
            mini = min(mini, bin(m).count("1"))
        return mini
    except Exception:
        return None


def options_fractions(k, res_max=99):
    """Liste de (libellé, générateurs) pour k facteurs."""
    out = []
    for gens in FRACTIONS.get(k, []):
        p = len(gens)
        r = resolution_fraction(k, gens)
        if r is None or r > res_max:
            continue
        out.append((f"2^({k}−{p}) — {2 ** (k - p)} essais — résolution {ROMAIN.get(r, r)}", gens))
    return out


# =============================================================================
#  Plan d'expériences
# =============================================================================
class Plan:
    """Un plan : colonnes (ordre d'exécution), facteurs, réponses saisies, informations de construction."""

    def __init__(self, famille, facteurs, cols, meta, titre="", lignes=None, tables=None):
        self.famille = famille
        self.facteurs = facteurs                  # liste de dict (nom, type, bas, haut | niveaux)
        self.cols = cols                          # nom de colonne -> liste de valeurs
        self.meta = meta or {}
        self.titre = titre
        self.lignes = list(lignes or [])
        self.tables = tables or {}                # tableaux descriptifs (alias, affectation des colonnes)
        self.reponses = []

    @property
    def n(self):
        return len(next(iter(self.cols.values()))) if self.cols else 0

    @property
    def noms_facteurs(self):
        return [f["nom"] for f in self.facteurs]

    def ordre_colonnes(self):
        return [c for c in SYSTEME if c in self.cols] + self.noms_facteurs + self.reponses

    def nom_reponse_libre(self, base="Y"):
        i = len(self.reponses) + 1
        while f"{base}{i}" in self.cols:
            i += 1
        return f"{base}{i}"

    def verifier_nom_reponse(self, nom):
        nom = (nom or "").strip()
        if not nom:
            raise ErreurDOE("Le nom de la réponse est vide.")
        if nom.lower() in RESERVES or nom in self.noms_facteurs:
            raise ErreurDOE("Ce nom est déjà utilisé par une colonne du plan.")
        if nom in self.reponses:
            raise ErreurDOE("Il existe déjà une réponse de ce nom.")
        if any(c in nom for c in ',;"\n\t'):
            raise ErreurDOE("Le nom d'une réponse ne peut contenir ni virgule, ni point-virgule, ni guillemet.")
        return nom

    def ajouter_reponse(self, nom, valeurs=None):
        nom = self.verifier_nom_reponse(nom)
        self.cols[nom] = list(valeurs) if valeurs is not None else [None] * self.n
        self.reponses.append(nom)
        return nom

    def supprimer_reponse(self, nom):
        if nom in self.reponses:
            self.reponses.remove(nom)
            del self.cols[nom]

    def renommer_reponse(self, ancien, nouveau):
        if ancien == nouveau:
            return
        nouveau = self.verifier_nom_reponse(nouveau)
        self.cols = {(nouveau if c == ancien else c): v for c, v in self.cols.items()}
        self.reponses = [nouveau if r == ancien else r for r in self.reponses]

    def nb_renseignes(self, nom):
        return sum(1 for v in self.cols.get(nom, []) if v is not None)

    def meta_r(self):
        return {**self.meta, "famille": self.famille, "facteurs": self.facteurs, "reponses": list(self.reponses)}

    def ecrire_csv_r(self, chemin):
        cols = self.ordre_colonnes()
        with open(chemin, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            for i in range(self.n):
                w.writerow(["" if self.cols[c][i] is None else repr(self.cols[c][i]) if isinstance(self.cols[c][i], float) else self.cols[c][i] for c in cols])

    def exporter_csv(self, chemin):
        """CSV lisible par un tableur français (point-virgule, virgule décimale, UTF-8 avec marque)."""
        cols = self.ordre_colonnes()
        with open(chemin, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(cols)
            for i in range(self.n):
                w.writerow(["" if self.cols[c][i] is None else (fmt_export(self.cols[c][i])) for c in cols])

    def vers_dict(self):
        return {"famille": self.famille, "facteurs": self.facteurs, "cols": self.cols, "meta": self.meta, "titre": self.titre,
                "lignes": self.lignes, "tables": self.tables, "reponses": self.reponses}

    @staticmethod
    def depuis_dict(d):
        if not isinstance(d, dict) or not isinstance(d.get("cols"), dict) or not isinstance(d.get("facteurs"), list):
            raise ErreurDOE("Plan invalide dans le projet.")
        cols = {str(k): list(v) for k, v in d["cols"].items()}
        longueurs = {len(v) for v in cols.values()}
        if len(longueurs) > 1:
            raise ErreurDOE("Plan invalide dans le projet (colonnes de longueurs différentes).")
        p = Plan(str(d.get("famille", "")), d["facteurs"], cols, d.get("meta") or {}, str(d.get("titre", "")), d.get("lignes"), d.get("tables"))
        p.reponses = [r for r in d.get("reponses", []) if r in cols]
        return p


def fmt_export(v):
    if isinstance(v, float):
        return repr(v).replace(".", ",")
    return str(v)


def convertir_reponse(texte):
    """Texte saisi -> float ou None (vide)."""
    t = str(texte).strip().replace(" ", "").replace(" ", "")
    if t == "" or t.lower() in ("na", "nan", "n/a", "—", "-"):
        return None
    try:
        v = float(t.replace(",", "."))
    except ValueError:
        raise ErreurDOE(f"« {texte} » n'est pas un nombre.")
    if math.isnan(v) or math.isinf(v):
        raise ErreurDOE(f"« {texte} » n'est pas un nombre valide.")
    return v


def lire_tableau_texte(texte):
    """Texte collé ou fichier -> (noms ou None, lignes de nombres). Séparateurs : tabulation, point-virgule, espace, virgule."""
    lignes = [l for l in texte.replace("\r", "").split("\n") if l.strip()]
    if not lignes:
        raise ErreurDOE("Aucune donnée.")
    premiere = lignes[0]
    if "\t" in premiere:
        sep = "\t"
    elif ";" in premiere:
        sep = ";"
    elif "," in premiere and "." in premiere:
        sep = ","                                  # virgule séparatrice, point décimal
    else:
        sep = None                                 # espaces ; la virgule est alors la décimale

    def decouper(l):
        if sep:
            return [c.strip() for c in l.split(sep)]
        return l.split()
    rangees = [decouper(l) for l in lignes]
    larg = max(len(r) for r in rangees)
    rangees = [r + [""] * (larg - len(r)) for r in rangees]
    entete = None
    try:
        for c in rangees[0]:
            if c.strip():
                convertir_reponse(c)
    except ErreurDOE:
        entete = [c.strip() for c in rangees[0]]
        rangees = rangees[1:]
    vals = []
    for r in rangees:
        vals.append([convertir_reponse(c) for c in r])
    if not vals:
        raise ErreurDOE("Aucune valeur numérique.")
    return entete, vals


# =============================================================================
#  Résultats de R
# =============================================================================
class Resultat:
    """Résultat d'une analyse R : titre, lignes de texte, tableaux, figures PNG, remarques."""

    def __init__(self, brut, dossier):
        self.brut = brut
        self.dossier = Path(dossier)
        self.commande = brut.get("commande", "")
        self.titre = str(brut.get("titre", ""))
        self.conclusion = str(brut.get("conclusion") or "")
        self.lignes = [str(x) for x in _liste(brut.get("lignes"))]
        self.avertissements = [str(x) for x in _liste(brut.get("avertissements"))]
        self.optimum = brut.get("optimum") if isinstance(brut.get("optimum"), dict) else None
        self.figures = [self.dossier / str(f) for f in _liste(brut.get("figures"))
                        if Path(str(f)).name == str(f) and (self.dossier / str(f)).is_file()]
        self.tables = {}
        self.libres = {str(x) for x in _liste(brut.get("tables_libres"))}
        tabs = brut.get("tables")
        if isinstance(tabs, dict):
            for nom, t in tabs.items():
                if not isinstance(t, dict) or not t:
                    continue
                cols = list(t.keys())
                vals = [_liste(t[c]) for c in cols]
                self.tables[nom] = (cols, [list(r) for r in zip(*vals)])

    def texte(self, titres=True):
        parts = []
        if titres and self.titre:
            parts.append(self.titre + "\n" + "=" * len(self.titre))
        if self.lignes:
            parts.append("\n".join(self.lignes))
        for nom, (cols, rows) in self.tables.items():
            parts.append(libelle(nom) + "\n" + table_texte(cols, rows, libre=nom in self.libres))
        if self.avertissements:
            parts.append("\n".join("Remarque : " + a for a in self.avertissements))
        return "\n\n".join(parts)

    def liberer(self):
        shutil.rmtree(self.dossier, ignore_errors=True)


_dossiers_temp = []


@atexit.register
def _nettoyer_temp():
    for d in _dossiers_temp:
        shutil.rmtree(d, ignore_errors=True)


def verifier_r():
    rs = shutil.which("Rscript")
    if not rs:
        return False, "R est introuvable. Installer : sudo apt install -y r-base r-cran-jsonlite"
    manquants = [f for f in FICHIERS_R if not (DOSSIER / f).is_file()]
    if manquants:
        return False, "Fichier(s) R manquant(s) dans le dossier du programme : " + ", ".join(manquants)
    try:
        r = subprocess.run([rs, "-e", "library(jsonlite)"], capture_output=True, text=True, timeout=60)
    except Exception as e:
        return False, f"R inutilisable : {e}"
    if r.returncode != 0:
        return False, "Le package R « jsonlite » est absent. Installer : sudo apt install -y r-cran-jsonlite"
    return True, "R prêt."


def appeler_r(commande, params, plan=None, timeout=240, conserver=True):
    """Lance cli_doe.R ; retourne (dict du résultat, dossier). Le dossier est supprimé si conserver=False."""
    rs = shutil.which("Rscript")
    if not rs:
        raise ErreurDOE("R est introuvable. Installer : sudo apt install -y r-base r-cran-jsonlite")
    dossier = Path(tempfile.mkdtemp(prefix="doe_"))
    _dossiers_temp.append(dossier)
    try:
        corps = {"commande": commande, **params}
        if plan is not None:
            plan.ecrire_csv_r(dossier / "plan.csv")
            corps["meta"] = plan.meta_r()
        (dossier / "params.json").write_text(json.dumps(corps, ensure_ascii=False), encoding="utf-8")
        env = dict(os.environ, LC_ALL="C.UTF-8", LANG="C.UTF-8")
        try:
            p = subprocess.run([rs, str(SCRIPT_R), str(dossier)], capture_output=True, text=True, timeout=timeout, env=env)
        except subprocess.TimeoutExpired:
            raise ErreurDOE("Le calcul R a dépassé le délai autorisé.")
        sortie = dossier / "resultat.json"
        if not sortie.is_file():
            detail = (p.stderr or p.stdout or "").strip().splitlines()
            raise ErreurDOE("R n'a produit aucun résultat. " + (detail[-1] if detail else ""))
        brut = json.loads(sortie.read_text(encoding="utf-8"))
        if not brut.get("ok"):
            raise ErreurDOE(str(brut.get("erreur", "Erreur de calcul.")))
        if not conserver:
            shutil.rmtree(dossier, ignore_errors=True)
        return brut, dossier
    except Exception:
        shutil.rmtree(dossier, ignore_errors=True)
        raise


def executer_r(plan, commande, params, echelle=1.0, alpha=0.05):
    brut, dossier = appeler_r(commande, {**params, "echelle": echelle, "alpha_test": alpha}, plan)
    return Resultat(brut, dossier)


def generer_plan(famille, commande, params, facteurs):
    """Construit un plan avec R. facteurs : définition envoyée à R (conservée dans le plan)."""
    brut, dossier = appeler_r(commande, params, None, conserver=False)
    cols = {str(k): _liste(v) for k, v in (brut.get("plan") or {}).items()}
    if not cols:
        raise ErreurDOE("R n'a produit aucun plan.")
    n = len(next(iter(cols.values())))
    if n > MAX_ESSAIS:
        raise ErreurDOE(f"Le plan compte {n} essais : le maximum est {MAX_ESSAIS}.")
    tables = {}
    for nom, t in (brut.get("tables") or {}).items():
        if isinstance(t, dict):
            tables[nom] = {c: _liste(v) for c, v in t.items()}
    meta = brut.get("meta") or {}
    return Plan(famille, facteurs, cols, meta, str(brut.get("titre", "")), [str(x) for x in _liste(brut.get("lignes"))], tables)


# =============================================================================
#  Projet : enregistrement et reprise
# =============================================================================
RE_FIGURE_PROJET = re.compile(r"^figures/(\d{1,5})/([A-Za-z0-9_.\-]{1,80}\.png)$")
LIMITE_JSON_PROJET = 200_000_000
LIMITE_PNG_PROJET = 30_000_000


def enregistrer_projet(chemin, plans, entrees, meta=None):
    """plans : {famille: Plan} ; entrees : liste de (famille, heure, Resultat, inclus)."""
    projet = {"format": FORMAT_PROJET, "version": VERSION_PROJET, "enregistre": datetime.datetime.now().isoformat(timespec="seconds"),
              "meta": meta or {}, "plans": {k: p.vers_dict() for k, p in plans.items() if p is not None}, "journal": []}
    fichiers = []
    for i, (fam, heure, res, inclus) in enumerate(entrees):
        for f in res.figures:
            fichiers.append((f"figures/{i}/{f.name}", f))
        projet["journal"].append({"famille": fam, "heure": heure, "inclus": bool(inclus), "brut": res.brut})
    chemin = Path(chemin)
    tmp = chemin.with_name(chemin.name + ".tmp")
    try:
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("projet.json", json.dumps(projet, ensure_ascii=False))
            for nom, f in fichiers:
                z.write(f, nom)
        os.replace(tmp, chemin)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def lire_projet(chemin):
    invalide = ErreurDOE("Ce fichier n'est pas un projet valide du programme DOE.")
    try:
        z = zipfile.ZipFile(chemin)
    except (zipfile.BadZipFile, OSError):
        raise invalide
    with z:
        infos = {i.filename: i for i in z.infolist()}
        if "projet.json" not in infos or infos["projet.json"].file_size > LIMITE_JSON_PROJET or len(infos) > 20_000:
            raise invalide
        try:
            projet = json.loads(z.read("projet.json").decode("utf-8"))
        except (ValueError, UnicodeDecodeError):
            raise invalide
        if not isinstance(projet, dict) or projet.get("format") != FORMAT_PROJET:
            raise invalide
        version = projet.get("version")
        if not isinstance(version, int) or version > VERSION_PROJET:
            raise ErreurDOE("Ce projet a été créé par une version plus récente du programme.")
        plans = {}
        for fam, d in (projet.get("plans") or {}).items():
            plans[fam] = Plan.depuis_dict(d)
        entrees = []
        journal = projet.get("journal")
        for i, e in enumerate(journal if isinstance(journal, list) else []):
            if not isinstance(e, dict) or not isinstance(e.get("brut"), dict):
                raise invalide
            dossier = Path(tempfile.mkdtemp(prefix="doe_"))
            _dossiers_temp.append(dossier)
            for nom, info in infos.items():
                m = RE_FIGURE_PROJET.match(nom)
                if m and int(m.group(1)) == i and info.file_size <= LIMITE_PNG_PROJET:
                    (dossier / m.group(2)).write_bytes(z.read(nom))
            entrees.append((str(e.get("famille", "")), str(e.get("heure", "")), Resultat(e["brut"], dossier), bool(e.get("inclus", True))))
    meta = projet.get("meta") if isinstance(projet.get("meta"), dict) else {}
    return meta, plans, entrees


# =============================================================================
#  Interface — éléments communs
# =============================================================================
def charger_image(chemin, largeur_max, hauteur_max=None, agrandir=False):
    """PhotoImage ajustée à largeur_max × hauteur_max (proportions conservées)."""
    import tkinter as tk
    try:
        from PIL import Image, ImageTk
        im = Image.open(chemin)
        ratio = largeur_max / im.width
        if hauteur_max:
            ratio = min(ratio, hauteur_max / im.height)
        if not agrandir:
            ratio = min(ratio, 1.0)
        if abs(ratio - 1.0) > 0.01:
            im = im.resize((max(1, int(im.width * ratio)), max(1, int(im.height * ratio))), Image.LANCZOS)
        return ImageTk.PhotoImage(im)
    except Exception:
        pass
    photo = tk.PhotoImage(file=str(chemin))
    ratio = largeur_max / photo.width()
    if hauteur_max:
        ratio = min(ratio, hauteur_max / photo.height())
    if ratio < 0.99:
        photo = photo.subsample(math.ceil(1 / ratio))
    return photo


def echelle_ecran(root):
    return max(1.0, min(2.0, root.winfo_screenwidth() * 0.96 / 900, root.winfo_screenheight() * 0.88 / 600))


def _placer_fenetre(win, parent, largeur_min=460):
    win.update_idletasks()
    w, h = max(win.winfo_reqwidth(), largeur_min), win.winfo_reqheight() + 8
    win.minsize(w, h)
    try:
        x = parent.winfo_rootx() + (parent.winfo_width() - w) // 2
        y = parent.winfo_rooty() + (parent.winfo_height() - h) // 3
    except Exception:
        x = y = 100
    win.geometry(f"{w}x{h}+{max(x, 0)}+{max(y, 0)}")


def boite_dialogue(parent, titre, message, boutons=("OK",), defaut=0, details=None):
    """Boîte de dialogue à boutons personnalisés ; retourne l'indice du bouton choisi (le dernier si fermée)."""
    import tkinter as tk
    from tkinter import ttk
    win = tk.Toplevel(parent)
    win.title(titre)
    win.transient(parent)
    choix = {"i": len(boutons) - 1}

    def fin(i):
        choix["i"] = i
        win.destroy()
    cadre = ttk.Frame(win, padding=(16, 14, 16, 12))
    cadre.pack(fill="both", expand=True)
    barre = ttk.Frame(cadre)
    barre.pack(side="bottom", anchor="e", pady=(14, 0))
    boutons_tk = []
    for i, texte in enumerate(boutons):
        b = tk.Button(barre, text=texte, command=lambda i=i: fin(i), padx=14, pady=3, default="active" if i == defaut else "normal")
        b.pack(side="left", padx=4)
        boutons_tk.append(b)
    if message:
        ttk.Label(cadre, text=message, wraplength=520, justify="left", font=("TkDefaultFont", 10, "bold")).pack(anchor="w", fill="x", pady=(0, 6))
    if details:
        grille = ttk.Frame(cadre)
        grille.pack(anchor="w", fill="x")
        for r, (lib, val) in enumerate(details):
            ttk.Label(grille, text=f"{lib} :").grid(row=r, column=0, sticky="nw", padx=(0, 6), pady=2)
            ttk.Label(grille, text=val, wraplength=440, justify="left").grid(row=r, column=1, sticky="nw", pady=2)
    win.bind("<Escape>", lambda e: fin(len(boutons) - 1))
    win.bind("<Return>", lambda e: fin(defaut))
    win.protocol("WM_DELETE_WINDOW", lambda: fin(len(boutons) - 1))
    _placer_fenetre(win, parent)
    win.update()
    try:
        win.grab_set()
        boutons_tk[defaut].focus_set()
    except Exception:
        pass
    parent.wait_window(win)
    return choix["i"]


def choisir_fichier_ecriture(parent, **options):
    from tkinter import filedialog
    while True:
        try:
            chemin = filedialog.asksaveasfilename(parent=parent, confirmoverwrite=False, **options)
        except Exception:
            chemin = filedialog.asksaveasfilename(parent=parent, **options)
        if not chemin:
            return ""
        p = Path(chemin)
        if not p.exists():
            return chemin
        rep = boite_dialogue(parent, TITRE, "Ce fichier existe déjà. Voulez-vous le remplacer ?", ("Remplacer", "Choisir un autre nom", "Annuler"),
                             defaut=1, details=[("Dossier", str(p.parent)), ("Fichier", p.name)])
        if rep == 0:
            return chemin
        if rep == 2:
            return ""
        options["initialdir"], options["initialfile"] = str(p.parent), p.name


class CadreDefilant:
    """Zone à défilement vertical ; les widgets se placent dans .interieur."""

    def __init__(self, parent):
        import tkinter as tk
        from tkinter import ttk
        self.frame = ttk.Frame(parent)
        self.canvas = tk.Canvas(self.frame, highlightthickness=0, borderwidth=0)
        sb = ttk.Scrollbar(self.frame, orient="vertical", command=self.canvas.yview)
        self.interieur = ttk.Frame(self.canvas)
        self.interieur.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self._fen = self.canvas.create_window((0, 0), window=self.interieur, anchor="nw")
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self._fen, width=e.width))
        self.canvas.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.bind("<Enter>", self._lier)
        self.canvas.bind("<Leave>", self._delier)
        self.interieur.bind("<Enter>", self._lier)

    def _lier(self, _e=None):
        self.canvas.bind_all("<MouseWheel>", self._roulette)
        self.canvas.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-2, "units"))
        self.canvas.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(2, "units"))

    def _delier(self, _e=None):
        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            self.canvas.unbind_all(seq)

    def _roulette(self, e):
        self.canvas.yview_scroll(-1 if e.delta > 0 else 1, "units")


class PanneauGraphiques:
    """Affiche des graphiques PNG ajustés à la place disponible ; ◀ ▶ pour passer de l'un à l'autre."""

    def __init__(self, parent, plein_ecran=None, fermer=None):
        import tkinter as tk
        from tkinter import ttk
        self.frame = ttk.Frame(parent)
        self.figures, self.i, self.photo, self._job, self._taille = [], 0, None, None, (0, 0)
        barre = ttk.Frame(self.frame)
        barre.pack(fill="x", pady=(0, 4))
        ttk.Button(barre, text="◀", width=3, command=lambda: self.aller(-1)).pack(side="left")
        self.lbl = ttk.Label(barre, text="", width=16, anchor="center")
        self.lbl.pack(side="left")
        ttk.Button(barre, text="▶", width=3, command=lambda: self.aller(1)).pack(side="left")
        ttk.Button(barre, text="Enregistrer le graphique…", command=self.enregistrer).pack(side="left", padx=10)
        if plein_ecran:
            ttk.Button(barre, text="Plein écran", command=plein_ecran).pack(side="left")
        if fermer:
            ttk.Button(barre, text="Fermer", command=fermer).pack(side="right")
        self.zone = tk.Canvas(self.frame, background="white", highlightthickness=1, highlightbackground="#999999")
        self.zone.pack(fill="both", expand=True)
        self.zone.bind("<Configure>", self._reajuster)

    def definir(self, figures):
        self.figures = list(figures)
        self.i = 0
        self.afficher()

    def aller(self, pas):
        if len(self.figures) > 1:
            self.i = (self.i + pas) % len(self.figures)
            self.afficher()

    def afficher(self):
        self.zone.delete("all")
        self.photo = None
        w, h = self.zone.winfo_width(), self.zone.winfo_height()
        self._taille = (w, h)
        if not self.figures:
            self.lbl.configure(text="")
            return
        self.lbl.configure(text=f"Graphique {self.i + 1} / {len(self.figures)}")
        try:
            self.photo = charger_image(self.figures[self.i], max(200, w - 8), max(150, h - 8), agrandir=True)
            self.zone.create_image(w // 2, h // 2, image=self.photo)
        except Exception as e:
            self.zone.create_text(10, 10, anchor="nw", text=f"Affichage impossible : {e}")

    def _reajuster(self, _e=None):
        if self._job:
            self.frame.after_cancel(self._job)
        self._job = self.frame.after(200, self._maj)

    def _maj(self):
        self._job = None
        w, h = self.zone.winfo_width(), self.zone.winfo_height()
        if self.figures and (abs(w - self._taille[0]) > 20 or abs(h - self._taille[1]) > 20):
            self.afficher()

    def enregistrer(self):
        if not self.figures:
            return
        src = self.figures[self.i]
        chemin = choisir_fichier_ecriture(self.frame.winfo_toplevel(), title="Enregistrer le graphique", defaultextension=".png", initialfile=src.name,
                                          filetypes=[("Image PNG", "*.png")], initialdir=str(Path.home()))
        if chemin:
            shutil.copy(src, chemin)


class Visionneuse:
    def __init__(self, parent, figures):
        import tkinter as tk
        self.win = tk.Toplevel(parent)
        self.win.title(f"{TITRE} — graphiques")
        sw, sh = self.win.winfo_screenwidth(), self.win.winfo_screenheight()
        self.win.geometry(f"{int(sw * 0.96)}x{int(sh * 0.90)}+{int(sw * 0.02)}+0")
        self.panneau = PanneauGraphiques(self.win, fermer=self.win.destroy)
        self.panneau.frame.pack(fill="both", expand=True, padx=8, pady=8)
        self.win.bind("<Left>", lambda e: self.panneau.aller(-1))
        self.win.bind("<Right>", lambda e: self.panneau.aller(1))
        self.win.bind("<Escape>", lambda e: self.win.destroy())
        self.win.update_idletasks()
        self.panneau.definir(figures)


class EntreeJournal:
    def __init__(self, famille, zone, res, heure=None, inclus=True):
        self.famille, self.zone, self.res, self.inclus = famille, zone, res, inclus
        self.heure = heure or datetime.datetime.now().isoformat(timespec="seconds")


class ZoneResultat:
    """Texte à chasse fixe à gauche, graphiques à droite (séparateur déplaçable)."""

    def __init__(self, parent, app, famille, nom):
        import tkinter as tk
        from tkinter import ttk, font as tkfont
        self.app, self.famille, self.nom = app, famille, nom
        self.resultat = None
        self.frame = ttk.Frame(parent)
        barre = ttk.Frame(self.frame)
        barre.pack(fill="x", pady=(6, 2))
        ttk.Button(barre, text="Enregistrer le résultat…", command=self.enregistrer).pack(side="left")
        ttk.Button(barre, text="Graphiques en plein écran…", command=self.agrandir).pack(side="left", padx=8)
        self.pan = ttk.PanedWindow(self.frame, orient="horizontal")
        self.pan.pack(fill="both", expand=True)
        gauche = ttk.Frame(self.pan)
        self.texte = tk.Text(gauche, wrap="none", height=10, width=40, state="disabled", padx=8, pady=6,
                             font=tkfont.nametofont("TkFixedFont"), background="white", relief="solid", borderwidth=1)
        sy = ttk.Scrollbar(gauche, orient="vertical", command=self.texte.yview)
        sx = ttk.Scrollbar(gauche, orient="horizontal", command=self.texte.xview)
        self.texte.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sx.pack(side="bottom", fill="x")
        sy.pack(side="right", fill="y")
        self.texte.pack(side="left", fill="both", expand=True)
        self.pan.add(gauche, weight=2)
        self.graphes = PanneauGraphiques(self.pan, plein_ecran=self.agrandir)
        self._graphes_visibles = False
        gras = tkfont.nametofont("TkFixedFont").copy()
        gras.configure(weight="bold")
        self.texte.tag_configure("titre", font=gras, foreground="#1f4e8c")
        self.texte.tag_configure("remarque", foreground="#8a5a00")
        self.texte.tag_configure("erreur", foreground="#b00020")
        self.texte.tag_configure("section", font=gras)

    def _montrer_graphes(self, oui):
        if oui and not self._graphes_visibles:
            self.pan.add(self.graphes.frame, weight=3)
            self._graphes_visibles = True
            self.frame.after(80, self._placer_separateur)
            self.frame.after(500, lambda: self._placer_separateur(1))
        elif not oui and self._graphes_visibles:
            self.pan.forget(self.graphes.frame)
            self._graphes_visibles = False

    def _placer_separateur(self, essai=0):
        try:
            larg = self.pan.winfo_width()
            if larg < 200:
                raise ValueError("zone non affichée")
            pos = self.pan.sashpos(0)
            if essai == 0 or pos < 0.15 * larg or pos > 0.85 * larg:
                self.pan.sashpos(0, int(larg * 0.42))
        except Exception:
            if essai < 8:
                self.frame.after(150, lambda: self._placer_separateur(essai + 1))

    def agrandir(self):
        from tkinter import messagebox
        if self.resultat is None or not self.resultat.figures:
            messagebox.showinfo(TITRE, "Aucun graphique à agrandir.")
            return
        Visionneuse(self.app.root, self.resultat.figures)

    def _ecrire(self, contenu, tag=None):
        self.texte.insert("end", contenu, tag or ())

    def vider(self):
        self.texte.configure(state="normal")
        self.texte.delete("1.0", "end")
        self.texte.configure(state="disabled")

    def reinitialiser(self):
        self.resultat = None
        self._montrer_graphes(False)
        self.graphes.definir([])
        self.vider()

    def afficher_erreur(self, message):
        self.reinitialiser()
        self.texte.configure(state="normal")
        self._ecrire("Erreur : " + message + "\n", "erreur")
        self.texte.configure(state="disabled")

    def afficher(self, res):
        self.resultat = res
        self.vider()
        t = self.texte
        t.configure(state="normal")
        self._ecrire(res.titre + "\n", "titre")
        self._ecrire("=" * len(res.titre) + "\n\n")
        if res.lignes:
            self._ecrire("\n".join(res.lignes) + "\n\n")
        for nom, (cols, rows) in res.tables.items():
            self._ecrire(libelle(nom) + "\n", "section")
            self._ecrire(table_texte(cols, rows, libre=nom in res.libres) + "\n\n")
        for a in res.avertissements:
            self._ecrire("Remarque : " + a + "\n", "remarque")
        t.configure(state="disabled")
        t.yview_moveto(0)
        self._montrer_graphes(bool(res.figures))
        self.graphes.definir(res.figures)

    def enregistrer(self):
        from tkinter import messagebox
        if self.resultat is None:
            messagebox.showinfo(TITRE, "Aucun résultat à enregistrer.")
            return
        chemin = choisir_fichier_ecriture(self.app.root, title="Enregistrer le résultat", defaultextension=".txt",
                                          filetypes=[("Texte", "*.txt")], initialdir=str(Path.home()))
        if not chemin:
            return
        p = Path(chemin)
        p.write_text(self.resultat.texte() + "\n", encoding="utf-8")
        for f in self.resultat.figures:
            shutil.copy(f, p.with_name(f"{p.stem}_{f.name}"))
        self.app.statut(f"Résultat enregistré : {p.name}" + (" (figures PNG à côté)" if self.resultat.figures else ""))


class DialogueCollage:
    """Saisie de réponses : collage de colonnes, ou lecture d'un fichier texte / CSV."""

    def __init__(self, famille):
        import tkinter as tk
        from tkinter import ttk, font as tkfont
        self.fam = famille
        self.win = tk.Toplevel(famille.app.root)
        self.win.title(f"{TITRE} — saisir des réponses")
        self.win.transient(famille.app.root)
        cadre = ttk.Frame(self.win, padding=12)
        cadre.pack(fill="both", expand=True)
        ttk.Label(cadre, wraplength=560, justify="left",
                  text="Coller ici une colonne par réponse (un essai par ligne). Une première ligne de noms est facultative : "
                       "sans noms, les valeurs vont aux réponses du plan, dans l'ordre ; avec des noms, une réponse de même nom est remplacée "
                       "et toute autre est créée.").pack(anchor="w")
        self.v_ordre = tk.StringVar(value="execution")
        ligne = ttk.Frame(cadre)
        ligne.pack(anchor="w", pady=6)
        ttk.Label(ligne, text="Ordre des lignes :").pack(side="left")
        ttk.Radiobutton(ligne, text="ordre d'exécution", variable=self.v_ordre, value="execution").pack(side="left", padx=6)
        ttk.Radiobutton(ligne, text="ordre standard", variable=self.v_ordre, value="standard").pack(side="left", padx=6)
        self.texte = tk.Text(cadre, width=70, height=14, font=tkfont.nametofont("TkFixedFont"), undo=True)
        self.texte.pack(fill="both", expand=True)
        barre = ttk.Frame(cadre)
        barre.pack(fill="x", pady=(8, 0))
        ttk.Button(barre, text="Ouvrir un fichier…", command=self.ouvrir).pack(side="left")
        ttk.Button(barre, text="Valider", command=self.valider).pack(side="right")
        ttk.Button(barre, text="Annuler", command=self.win.destroy).pack(side="right", padx=6)
        self.win.bind("<Escape>", lambda e: self.win.destroy())
        _placer_fenetre(self.win, famille.app.root, 620)
        self.texte.focus_set()

    def ouvrir(self):
        from tkinter import filedialog, messagebox
        chemin = filedialog.askopenfilename(parent=self.win, title="Ouvrir un fichier de réponses",
                                            filetypes=[("Texte et CSV", "*.txt *.csv *.tsv *.dat"), ("Tous les fichiers", "*.*")])
        if not chemin:
            return
        try:
            contenu = Path(chemin).read_bytes()
            try:
                texte = contenu.decode("utf-8-sig")
            except UnicodeDecodeError:
                texte = contenu.decode("latin-1")
        except OSError as e:
            messagebox.showerror(TITRE, str(e), parent=self.win)
            return
        self.texte.delete("1.0", "end")
        self.texte.insert("1.0", texte)

    def valider(self):
        from tkinter import messagebox
        try:
            self.fam.integrer_reponses(self.texte.get("1.0", "end"), self.v_ordre.get() == "standard")
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e), parent=self.win)
            return
        self.win.destroy()

# =============================================================================
#  Interface — une famille de plans (criblage, factoriel, Taguchi, surface, mélange)
# =============================================================================
FAMILLES = {
    "criblage": {"titre": "1. Criblage", "nmin": 2, "nmax": 31, "defaut": 5, "mot_n": "Nombre de facteurs", "mot_f": "Facteurs", "nom": "Facteur"},
    "factoriel": {"titre": "2. Factoriel", "nmin": 2, "nmax": 11, "defaut": 3, "mot_n": "Nombre de facteurs", "mot_f": "Facteurs", "nom": "Facteur"},
    "taguchi": {"titre": "3. Taguchi", "nmin": 2, "nmax": 31, "defaut": 4, "mot_n": "Nombre de facteurs", "mot_f": "Facteurs", "nom": "Facteur"},
    "surface": {"titre": "4. Surface", "nmin": 2, "nmax": 7, "defaut": 3, "mot_n": "Nombre de facteurs", "mot_f": "Facteurs", "nom": "Facteur"},
    "melange": {"titre": "5. Mélange", "nmin": 2, "nmax": 10, "defaut": 3, "mot_n": "Nombre de composants", "mot_f": "Composants", "nom": "Composant"},
}
MODELES = {
    "criblage": [("Effets principaux", "ordre1"), ("Effets principaux et interactions à 2 facteurs", "ordre2")],
    "factoriel": [("Effets principaux", "ordre1"), ("Jusqu'aux interactions à 2 facteurs", "ordre2"),
                  ("Jusqu'aux interactions à 3 facteurs", "ordre3"), ("Jusqu'aux interactions à 4 facteurs", "ordre4")],
    "surface": [("Linéaire", "lineaire"), ("Linéaire et interactions à 2 facteurs", "interactions"), ("Quadratique complet", "quadratique")],
    "melange": [("Linéaire", "lineaire"), ("Quadratique", "quadratique"), ("Cubique spécial", "cubique_special"), ("Cubique complet", "cubique")],
}
MODELE_DEFAUT = {"criblage": "ordre1", "factoriel": "ordre2", "surface": "quadratique", "melange": "quadratique"}
TYPES_RSB = [("Plus grand est meilleur", "plus_grand"), ("Plus petit est meilleur", "plus_petit"),
             ("Cible nominale (moyenne et dispersion)", "nominal_1"), ("Cible nominale (dispersion seule)", "nominal_2")]
BUTS = [("Maximiser", "max"), ("Minimiser", "min"), ("Atteindre une cible", "cible")]
TABLEAUX_TAGUCHI = ["L4", "L8", "L9", "L12", "L16", "L18", "L25", "L27", "L32"]
TYPES_SURFACE = [("Composite centré circonscrit (CCC)", "ccc"), ("Composite centré à faces centrées (CCF)", "ccf"),
                 ("Composite centré inscrit (CCI)", "cci"), ("Box-Behnken (3 à 5 facteurs)", "bb")]
ALPHAS = [("Rotatable", "rotatable"), ("Orthogonal", "orthogonal"), ("Valeur personnalisée", "personnalise")]
TYPES_MELANGE = [("Réseau simplexe", "reseau"), ("Centroïde simplexe", "centroide"), ("Sommets extrêmes (bornes sur les composants)", "sommets")]
TYPES_GRAPHE = [("Lignes de niveau", "contour"), ("Surface en 3D", "surface"), ("Les deux", "les_deux")]


def lettre(i):
    return chr(65 + i) if i < 26 else str(i + 1)


def code_de(libelle_, liste):
    return dict(liste).get(libelle_, liste[0][1])


def libelle_de(code, liste):
    return next((l for l, c in liste if c == code), liste[0][0])


class Famille:
    def __init__(self, app, cle):
        tk = app.tk
        self.app, self.cle, self.cfg = app, cle, FAMILLES[cle]
        self.plan = None
        self.zones = {}
        self.lignes_fact = []
        self.interactions = []
        self.dernier_optimum = None
        self.v_nb = tk.StringVar(value=str(self.cfg["defaut"]))
        self.opt = {"rep": tk.StringVar(value="1"), "centres": tk.StringVar(value="0"), "alea": tk.BooleanVar(value=True), "graine": tk.StringVar(),
                    "methode": tk.StringVar(), "nb_essais": tk.StringVar(value="Automatique"), "type": tk.StringVar(), "fraction": tk.StringVar(),
                    "gens": tk.StringVar(), "tableau": tk.StringVar(value="Automatique"), "nrep": tk.StringVar(value="3"),
                    "alpha_mode": tk.StringVar(value=ALPHAS[0][0]), "alpha": tk.StringVar(), "degre": tk.StringVar(value="2"),
                    "axiaux": tk.BooleanVar(value=False), "aretes": tk.BooleanVar(value=False), "faces": tk.BooleanVar(value=False)}
        if cle == "criblage":
            self.opt["methode"].set("Plackett-Burman")
        if cle == "factoriel":
            self.opt["type"].set("Plan complet")
        if cle == "surface":
            self.opt["type"].set(TYPES_SURFACE[0][0])
            self.opt["centres"].set("")
        if cle == "melange":
            self.opt["type"].set(TYPES_MELANGE[0][0])
            self.opt["centres"].set("1")
        self.v_rep = tk.StringVar()
        modeles = MODELES.get(cle)
        self.v_modele = tk.StringVar(value=libelle_de(MODELE_DEFAUT.get(cle), modeles) if modeles else "")
        self.v_reduire = tk.BooleanVar(value=False)
        self.v_aout = tk.StringVar(value="0,05")
        self.v_pseudo = tk.BooleanVar(value=True)
        self.v_rsb = tk.StringVar(value=TYPES_RSB[0][0])
        self.v_regroup = tk.BooleanVar(value=True)
        self.v_nconf = tk.StringVar(value="3")
        self.v_edit = tk.StringVar()
        self.v_ordre = tk.StringVar(value="Ordre d'exécution")
        self.v_ax, self.v_ay, self.v_graphe = tk.StringVar(), tk.StringVar(), tk.StringVar(value=TYPES_GRAPHE[0][0])
        self.v_ca, self.v_cb, self.v_cc = tk.StringVar(), tk.StringVar(), tk.StringVar()
        self.fixes_vars, self.prev_vars, self.optim_lignes = {}, {}, {}
        self.frac_options = []

    # ------------------------------------------------------------------ construction
    def construire(self, parent):
        ttk = self.app.ttk
        self.nb = ttk.Notebook(parent)
        self._page_plan(self.nb)
        self._page_analyse(self.nb)
        if self.cle in ("factoriel", "surface", "melange"):
            self._page_surfaces(self.nb)
        self._page_prevision(self.nb)
        self._maj_nb()
        return self.nb

    def _zone(self, parent, nom):
        z = ZoneResultat(parent, self.app, self, nom)
        z.frame.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        self.zones[nom] = z
        self.app.zones.append(z)
        return z

    # ------------------------------------------------------------------ page 1 : plan
    def _page_plan(self, nb):
        ttk = self.app.ttk
        f = ttk.Frame(nb)
        nb.add(f, text="1. Plan")
        pan = ttk.PanedWindow(f, orient="horizontal")
        pan.pack(fill="both", expand=True, padx=6, pady=6)
        gauche, droite = ttk.Frame(pan), ttk.Frame(pan)
        pan.add(gauche, weight=2)
        pan.add(droite, weight=5)
        self._definition(gauche)
        self._tableau(droite)

    def _definition(self, parent):
        tk, ttk = self.app.tk, self.app.ttk
        c = self.cfg
        top = ttk.Frame(parent)
        top.pack(fill="x", pady=(0, 4))
        ttk.Label(top, text=c["mot_n"]).pack(side="left")
        sp = ttk.Spinbox(top, from_=c["nmin"], to=c["nmax"], width=4, textvariable=self.v_nb, command=self._maj_nb)
        sp.pack(side="left", padx=6)
        sp.bind("<Return>", lambda e: self._maj_nb())
        sp.bind("<FocusOut>", lambda e: self._maj_nb())
        cadre = ttk.LabelFrame(parent, text=c["mot_f"])
        cadre.pack(fill="both", expand=True)
        self.defil = CadreDefilant(cadre)
        self.defil.frame.pack(fill="both", expand=True, padx=4, pady=4)
        opts = ttk.LabelFrame(parent, text="Options du plan")
        opts.pack(fill="x", pady=(6, 0))
        getattr(self, "_options_" + self.cle)(opts)
        barre = ttk.Frame(parent)
        barre.pack(fill="x", pady=(8, 0))
        ttk.Button(barre, text="Générer le plan", command=self.generer).pack(side="left")
        ttk.Button(barre, text="Effacer le plan", command=self.effacer_plan).pack(side="left", padx=6)

    def _lg(self, parent, r, texte, widget, col=0, span=1):
        ttk = self.app.ttk
        ttk.Label(parent, text=texte).grid(row=r, column=col, sticky="e", padx=(6, 4), pady=2)
        widget.grid(row=r, column=col + 1, columnspan=span, sticky="w", pady=2)

    def _commun(self, parent, r, centres=True, repetitions=True):
        ttk = self.app.ttk
        if repetitions:
            self._lg(parent, r, "Répétitions du plan", ttk.Entry(parent, textvariable=self.opt["rep"], width=6))
            r += 1
        if centres:
            self._lg(parent, r, "Points centraux", ttk.Entry(parent, textvariable=self.opt["centres"], width=6))
            r += 1
        ttk.Checkbutton(parent, text="Ordre d'exécution aléatoire", variable=self.opt["alea"]).grid(row=r, column=0, columnspan=3, sticky="w", padx=6, pady=2)
        self._lg(parent, r + 1, "Valeur d'initialisation (facultatif)", ttk.Entry(parent, textvariable=self.opt["graine"], width=8))
        return r + 2

    def _options_criblage(self, p):
        ttk = self.app.ttk
        cb = ttk.Combobox(p, textvariable=self.opt["methode"], state="readonly", width=24, values=["Plackett-Burman", "Fractionnaire 2^(k−p)"])
        self._lg(p, 0, "Méthode", cb, span=2)
        cb.bind("<<ComboboxSelected>>", lambda e: self._etat_options())
        self.cb_essais = ttk.Combobox(p, textvariable=self.opt["nb_essais"], state="readonly", width=12, values=["Automatique", "8", "12", "16", "20", "24", "32"])
        self._lg(p, 1, "Essais (Plackett-Burman)", self.cb_essais)
        self.cb_frac = ttk.Combobox(p, textvariable=self.opt["fraction"], state="readonly", width=38)
        self._lg(p, 2, "Fraction", self.cb_frac, span=2)
        self.cb_frac.bind("<<ComboboxSelected>>", lambda e: self._choisir_fraction())
        self.ent_gens = ttk.Entry(p, textvariable=self.opt["gens"], width=38)
        self._lg(p, 3, "Générateurs", self.ent_gens, span=2)
        self._commun(p, 4)

    def _options_factoriel(self, p):
        ttk = self.app.ttk
        cb = ttk.Combobox(p, textvariable=self.opt["type"], state="readonly", width=24, values=["Plan complet", "Plan fractionnaire"])
        self._lg(p, 0, "Type", cb, span=2)
        cb.bind("<<ComboboxSelected>>", lambda e: self._etat_options())
        self.cb_frac = ttk.Combobox(p, textvariable=self.opt["fraction"], state="readonly", width=38)
        self._lg(p, 1, "Fraction", self.cb_frac, span=2)
        self.cb_frac.bind("<<ComboboxSelected>>", lambda e: self._choisir_fraction())
        self.ent_gens = ttk.Entry(p, textvariable=self.opt["gens"], width=38)
        self._lg(p, 2, "Générateurs", self.ent_gens, span=2)
        self._commun(p, 3)

    def _options_taguchi(self, p):
        ttk = self.app.ttk
        self.cb_tab = ttk.Combobox(p, textvariable=self.opt["tableau"], state="readonly", width=14, values=["Automatique"] + TABLEAUX_TAGUCHI)
        self._lg(p, 0, "Tableau orthogonal", self.cb_tab)
        self._lg(p, 1, "Colonnes de réponses", ttk.Spinbox(p, from_=1, to=20, width=4, textvariable=self.opt["nrep"]))
        ttk.Label(p, text="(répétitions ou conditions de bruit)", foreground="#666666").grid(row=1, column=2, sticky="w")
        cad = ttk.LabelFrame(p, text="Interactions à étudier")
        cad.grid(row=2, column=0, columnspan=3, sticky="ew", padx=6, pady=4)
        self.v_ia, self.v_ib = self.app.tk.StringVar(), self.app.tk.StringVar()
        self.cb_ia = ttk.Combobox(cad, textvariable=self.v_ia, state="readonly", width=14, postcommand=lambda: self._maj_noms_inter())
        self.cb_ib = ttk.Combobox(cad, textvariable=self.v_ib, state="readonly", width=14, postcommand=lambda: self._maj_noms_inter())
        self.cb_ia.grid(row=0, column=0, padx=3, pady=3)
        ttk.Label(cad, text="×").grid(row=0, column=1)
        self.cb_ib.grid(row=0, column=2, padx=3)
        ttk.Button(cad, text="Ajouter", command=self.ajouter_interaction).grid(row=0, column=3, padx=3)
        self.liste_inter = self.app.tk.Listbox(cad, height=3, width=34, exportselection=False)
        self.liste_inter.grid(row=1, column=0, columnspan=3, sticky="ew", padx=3, pady=3)
        ttk.Button(cad, text="Retirer", command=self.retirer_interaction).grid(row=1, column=3, padx=3)
        self._commun(p, 3, centres=False, repetitions=False)

    def _options_surface(self, p):
        ttk = self.app.ttk
        cb = ttk.Combobox(p, textvariable=self.opt["type"], state="readonly", width=38, values=[l for l, _ in TYPES_SURFACE])
        self._lg(p, 0, "Type de plan", cb, span=2)
        cb.bind("<<ComboboxSelected>>", lambda e: self._etat_options())
        self.cb_alpha = ttk.Combobox(p, textvariable=self.opt["alpha_mode"], state="readonly", width=24, values=[l for l, _ in ALPHAS])
        self._lg(p, 1, "Distance axiale α (CCC)", self.cb_alpha, span=2)
        self.cb_alpha.bind("<<ComboboxSelected>>", lambda e: self._etat_options())
        self.ent_alpha = ttk.Entry(p, textvariable=self.opt["alpha"], width=8)
        self._lg(p, 2, "Valeur de α", self.ent_alpha)
        ttk.Label(p, text="Points centraux : vide = valeur usuelle.", foreground="#666666").grid(row=3, column=0, columnspan=3, sticky="w", padx=6)
        self._commun(p, 4)
        self._etat_options()

    def _options_melange(self, p):
        ttk = self.app.ttk
        cb = ttk.Combobox(p, textvariable=self.opt["type"], state="readonly", width=42, values=[l for l, _ in TYPES_MELANGE])
        self._lg(p, 0, "Type de plan", cb, span=2)
        cb.bind("<<ComboboxSelected>>", lambda e: self._etat_options())
        self.sp_degre = ttk.Spinbox(p, from_=1, to=4, width=4, textvariable=self.opt["degre"])
        self._lg(p, 1, "Degré du réseau", self.sp_degre)
        self.ck_ar = ttk.Checkbutton(p, text="Milieux des arêtes", variable=self.opt["aretes"])
        self.ck_fa = ttk.Checkbutton(p, text="Centres des faces", variable=self.opt["faces"])
        self.ck_ar.grid(row=2, column=0, columnspan=3, sticky="w", padx=6)
        self.ck_fa.grid(row=3, column=0, columnspan=3, sticky="w", padx=6)
        ttk.Checkbutton(p, text="Points axiaux", variable=self.opt["axiaux"]).grid(row=4, column=0, columnspan=3, sticky="w", padx=6)
        self._lg(p, 5, "Points centraux", ttk.Entry(p, textvariable=self.opt["centres"], width=6))
        self._lg(p, 6, "Répétitions du plan", ttk.Entry(p, textvariable=self.opt["rep"], width=6))
        ttk.Checkbutton(p, text="Ordre d'exécution aléatoire", variable=self.opt["alea"]).grid(row=7, column=0, columnspan=3, sticky="w", padx=6)
        self._lg(p, 8, "Valeur d'initialisation (facultatif)", ttk.Entry(p, textvariable=self.opt["graine"], width=8))
        self._etat_options()

    def _etat_options(self):
        def etat(w, actif):
            try:
                w.configure(state=("readonly" if isinstance(w, self.app.ttk.Combobox) else "normal") if actif else "disabled")
            except Exception:
                pass
        if self.cle == "criblage":
            pb = self.opt["methode"].get().startswith("Plackett")
            etat(self.cb_essais, pb)
            etat(self.cb_frac, not pb)
            etat(self.ent_gens, not pb)
        elif self.cle == "factoriel":
            fr = self.opt["type"].get() == "Plan fractionnaire"
            etat(self.cb_frac, fr)
            etat(self.ent_gens, fr)
        elif self.cle == "surface":
            ccc = code_de(self.opt["type"].get(), TYPES_SURFACE) == "ccc"
            etat(self.cb_alpha, ccc)
            etat(self.ent_alpha, ccc and code_de(self.opt["alpha_mode"].get(), ALPHAS) == "personnalise")
        elif self.cle == "melange":
            t = code_de(self.opt["type"].get(), TYPES_MELANGE)
            etat(self.sp_degre, t == "reseau")
            etat(self.ck_ar, t == "sommets")
            etat(self.ck_fa, t == "sommets")

    def _maj_fractions(self):
        if self.cle not in ("criblage", "factoriel"):
            return
        try:
            k = int(self.v_nb.get())
        except ValueError:
            return
        opts = options_fractions(k, res_max=4 if self.cle == "criblage" else 99)
        self.frac_options = opts
        self.cb_frac.configure(values=[l for l, _ in opts])
        if opts:
            if self.opt["fraction"].get() not in [l for l, _ in opts]:
                self.opt["fraction"].set(opts[0][0])
            self._choisir_fraction()
        else:
            self.opt["fraction"].set("")
            self.opt["gens"].set("")

    def _choisir_fraction(self):
        for l, g in self.frac_options:
            if l == self.opt["fraction"].get():
                self.opt["gens"].set(" ; ".join(g))
                return

    # ------------------------------------------------------------------ facteurs
    def _maj_nb(self):
        try:
            n = int(float(self.v_nb.get()))
        except ValueError:
            n = self.cfg["defaut"]
        n = max(self.cfg["nmin"], min(self.cfg["nmax"], n))
        self.v_nb.set(str(n))
        tk = self.app.tk
        while len(self.lignes_fact) < n:
            i = len(self.lignes_fact)
            m = self.cle == "melange"
            self.lignes_fact.append({"nom": tk.StringVar(value=f"{self.cfg['nom']} {lettre(i)}"), "a": tk.StringVar(value="0" if m else ""),
                                     "b": tk.StringVar(value="1" if m else ""), "niv": tk.StringVar()})
        del self.lignes_fact[n:]
        self._dessiner_facteurs()
        self._maj_fractions()
        if self.cle == "taguchi":
            self.interactions = [t for t in self.interactions if max(t) < n]
            self._maj_liste_inter()

    def _dessiner_facteurs(self):
        ttk = self.app.ttk
        inner = self.defil.interieur
        for w in inner.winfo_children():
            w.destroy()
        if self.cle == "taguchi":
            ent = ["", "Nom", "Niveaux (séparés par des points-virgules)"]
        elif self.cle == "melange":
            ent = ["", "Nom", "Borne basse", "Borne haute"]
        else:
            ent = ["", "Nom", "Niveau bas", "Niveau haut"]
        for j, t in enumerate(ent):
            ttk.Label(inner, text=t, font=("TkDefaultFont", 9, "bold")).grid(row=0, column=j, padx=3, pady=(2, 4), sticky="w")
        for i, v in enumerate(self.lignes_fact):
            ttk.Label(inner, text=lettre(i), width=3).grid(row=i + 1, column=0, padx=3, pady=1)
            ttk.Entry(inner, textvariable=v["nom"], width=18).grid(row=i + 1, column=1, padx=3, pady=1)
            if self.cle == "taguchi":
                ttk.Entry(inner, textvariable=v["niv"], width=26).grid(row=i + 1, column=2, padx=3, pady=1)
            else:
                ttk.Entry(inner, textvariable=v["a"], width=10).grid(row=i + 1, column=2, padx=3, pady=1)
                ttk.Entry(inner, textvariable=v["b"], width=10).grid(row=i + 1, column=3, padx=3, pady=1)

    def _collecter(self):
        """Définition des facteurs saisis -> liste de dict pour R ; contrôles de saisie."""
        res, vus = [], set()
        for i, v in enumerate(self.lignes_fact):
            nom = v["nom"].get().strip()
            if not nom:
                raise ErreurDOE(f"Le nom du {self.cfg['nom'].lower()} {lettre(i)} est vide.")
            if nom in vus:
                raise ErreurDOE(f"Le nom « {nom} » est utilisé plusieurs fois.")
            if nom.lower() in RESERVES or any(ch in nom for ch in '",;\n\t'):
                raise ErreurDOE(f"Le nom « {nom} » n'est pas permis (nom réservé, virgule, point-virgule ou guillemet).")
            vus.add(nom)
            a, b = v["a"].get().strip(), v["b"].get().strip()
            if self.cle == "taguchi":
                brut = v["niv"].get().strip()
                niv = [x.strip() for x in (brut.split(";") if ";" in brut else brut.split())] if brut else []
                niv = [x for x in niv if x]
                if len(niv) < 2:
                    raise ErreurDOE(f"« {nom} » : saisir au moins deux niveaux, séparés par des points-virgules.")
                if len(niv) == 4 or len(niv) > 5:
                    raise ErreurDOE(f"« {nom} » : 2, 3 ou 5 niveaux sont pris en charge.")
                if len(set(niv)) != len(niv):
                    raise ErreurDOE(f"« {nom} » : les niveaux doivent être différents.")
                res.append({"nom": nom, "niveaux": niv})
            elif self.cle == "melange":
                lo, hi = lire_nombre(a or "0", f"la borne basse de « {nom} »"), lire_nombre(b or "1", f"la borne haute de « {nom} »")
                if lo < 0 or hi > 1 or lo > hi:
                    raise ErreurDOE(f"« {nom} » : les bornes doivent vérifier 0 ≤ borne basse ≤ borne haute ≤ 1.")
                res.append({"nom": nom, "bas": lo, "haut": hi})
            else:
                if not a or not b:
                    raise ErreurDOE(f"« {nom} » : saisir le niveau bas et le niveau haut.")
                if est_nombre(a) and est_nombre(b):
                    lo, hi = lire_nombre(a), lire_nombre(b)
                    if lo >= hi:
                        raise ErreurDOE(f"« {nom} » : le niveau bas doit être inférieur au niveau haut.")
                    res.append({"nom": nom, "type": "num", "bas": lo, "haut": hi})
                else:
                    if self.cle == "surface":
                        raise ErreurDOE(f"« {nom} » : un plan de surface exige des facteurs quantitatifs (niveaux numériques).")
                    if a == b:
                        raise ErreurDOE(f"« {nom} » : les deux niveaux doivent être différents.")
                    res.append({"nom": nom, "type": "cat", "niveaux": [a, b]})
        return res

    def _entier(self, var, nom, minimum=0, vide=None):
        t = var.get().strip()
        if t == "":
            return vide
        v = lire_nombre(t, nom)
        if v is None or v != int(v) or v < minimum:
            raise ErreurDOE(f"{nom} : saisir un entier supérieur ou égal à {minimum}.")
        return int(v)

    def _params_plan(self, facteurs):
        o = self.opt
        graine = lire_nombre(o["graine"].get(), "la valeur d'initialisation")
        base = {"aleatoire": bool(o["alea"].get()), "graine": int(graine) if graine is not None else None}
        if self.cle in ("criblage", "factoriel"):
            base.update({"facteurs": facteurs, "repetitions": self._entier(o["rep"], "Répétitions", 1, 1), "points_centraux": self._entier(o["centres"], "Points centraux", 0, 0)})
            if self.cle == "criblage":
                if o["methode"].get().startswith("Plackett"):
                    ne = o["nb_essais"].get()
                    base.update({"methode": "pb", "nb_essais": 0 if ne == "Automatique" else int(ne)})
                else:
                    base.update({"methode": "fraction", "generateurs": self._gens()})
                return "plan_criblage", base
            if o["type"].get() == "Plan fractionnaire":
                base.update({"type": "fractionnaire", "generateurs": self._gens()})
            else:
                base["type"] = "complet"
            return "plan_factoriel", base
        if self.cle == "taguchi":
            tab = o["tableau"].get()
            base.update({"facteurs": facteurs, "interactions": [[a + 1, b + 1] for a, b in self.interactions], "tableau": "auto" if tab == "Automatique" else tab})
            self._entier(o["nrep"], "Colonnes de réponses", 1, 1)
            return "plan_taguchi", base
        if self.cle == "surface":
            t = code_de(o["type"].get(), TYPES_SURFACE)
            base.update({"facteurs": facteurs, "type": t, "repetitions": self._entier(o["rep"], "Répétitions", 1, 1)})
            pc = self._entier(o["centres"], "Points centraux", 0, None)
            if pc is not None:
                base["points_centraux"] = pc
            if t == "ccc":
                m = code_de(o["alpha_mode"].get(), ALPHAS)
                base["alpha_mode"] = m
                if m == "personnalise":
                    a = lire_nombre(o["alpha"].get(), "α")
                    if a is None or a <= 0:
                        raise ErreurDOE("Saisir une valeur de α strictement positive.")
                    base["alpha"] = a
            return "plan_surface", base
        t = code_de(o["type"].get(), TYPES_MELANGE)
        base.update({"composants": facteurs, "type": t, "repetitions": self._entier(o["rep"], "Répétitions", 1, 1), "points_centraux": self._entier(o["centres"], "Points centraux", 0, 1),
                     "axiaux": bool(o["axiaux"].get()), "aretes": bool(o["aretes"].get()), "faces": bool(o["faces"].get())})
        if t == "reseau":
            base["degre"] = self._entier(o["degre"], "Degré", 1, 2)
        return "plan_melange", base

    def _gens(self):
        g = [x.strip() for x in re.split(r"[;,\s]+", self.opt["gens"].get()) if x.strip()]
        if not g:
            raise ErreurDOE("Aucune fraction disponible pour ce nombre de facteurs : saisir les générateurs ou modifier le nombre de facteurs.")
        return g

    # ------------------------------------------------------------------ interactions Taguchi
    def _noms_saisis(self):
        return [v["nom"].get().strip() or f"{self.cfg['nom']} {lettre(i)}" for i, v in enumerate(self.lignes_fact)]

    def _maj_noms_inter(self):
        n = self._noms_saisis()
        self.cb_ia.configure(values=n)
        self.cb_ib.configure(values=n)

    def _maj_liste_inter(self):
        n = self._noms_saisis()
        self.liste_inter.delete(0, "end")
        for a, b in self.interactions:
            self.liste_inter.insert("end", f"{n[a]} × {n[b]}")

    def ajouter_interaction(self):
        from tkinter import messagebox
        n = self._noms_saisis()
        if self.v_ia.get() not in n or self.v_ib.get() not in n:
            return
        a, b = sorted((n.index(self.v_ia.get()), n.index(self.v_ib.get())))
        if a == b:
            messagebox.showinfo(TITRE, "Choisir deux facteurs différents.")
            return
        if (a, b) not in self.interactions:
            self.interactions.append((a, b))
        self._maj_liste_inter()

    def retirer_interaction(self):
        sel = self.liste_inter.curselection()
        if sel:
            del self.interactions[sel[0]]
            self._maj_liste_inter()

    # ------------------------------------------------------------------ génération
    def generer(self):
        from tkinter import messagebox
        try:
            facteurs = self._collecter()
            commande, params = self._params_plan(facteurs)
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        if self.plan is not None and any(self.plan.nb_renseignes(r) for r in self.plan.reponses):
            if boite_dialogue(self.app.root, TITRE, "Remplacer le plan actuel ? Les réponses saisies seront perdues.", ("Remplacer", "Annuler"), defaut=1) != 0:
                return
        self.app.lancer_plan(self, commande, params, facteurs)

    def recevoir_plan(self, plan):
        self.plan = plan
        if self.cle == "taguchi":
            for i in range(max(1, int(float(self.opt["nrep"].get() or 1)))):
                plan.ajouter_reponse(f"Y{i + 1}")
        else:
            plan.ajouter_reponse("Y1")
        self.v_rep.set("")
        self.dernier_optimum = None
        self.afficher_plan()
        self.maj_formulaires()
        self.app.statut(f"Plan généré : {plan.n} essais.")

    def effacer_plan(self):
        if self.plan is None:
            return
        if boite_dialogue(self.app.root, TITRE, "Effacer le plan et les réponses saisies ?", ("Effacer", "Annuler"), defaut=1) != 0:
            return
        self.plan = None
        self.afficher_plan()
        self.maj_formulaires()

    # ------------------------------------------------------------------ tableau du plan
    def _tableau(self, parent):
        tk, ttk = self.app.tk, self.app.ttk
        from tkinter import font as tkfont
        b1 = ttk.Frame(parent)
        b1.pack(fill="x")
        ttk.Label(b1, text="Réponse :").pack(side="left")
        self.cb_edit = ttk.Combobox(b1, textvariable=self.v_edit, state="readonly", width=16)
        self.cb_edit.pack(side="left", padx=4)
        ttk.Button(b1, text="Ajouter…", command=self.ajouter_reponse).pack(side="left", padx=2)
        ttk.Button(b1, text="Renommer…", command=self.renommer_reponse).pack(side="left", padx=2)
        ttk.Button(b1, text="Supprimer", command=self.supprimer_reponse).pack(side="left", padx=2)
        b2 = ttk.Frame(parent)
        b2.pack(fill="x", pady=4)
        ttk.Button(b2, text="Saisir ou coller des réponses…", command=self.coller).pack(side="left")
        ttk.Button(b2, text="Exporter le plan (CSV)…", command=self.exporter).pack(side="left", padx=6)
        cb = ttk.Combobox(b2, textvariable=self.v_ordre, state="readonly", width=18, values=["Ordre d'exécution", "Ordre standard"])
        cb.pack(side="left", padx=6)
        cb.bind("<<ComboboxSelected>>", lambda e: self.afficher_plan())
        self.lbl_info = ttk.Label(b2, text="Aucun plan.", foreground="#555555")
        self.lbl_info.pack(side="left", padx=6)
        pan = ttk.PanedWindow(parent, orient="vertical")
        pan.pack(fill="both", expand=True)
        corps = ttk.Frame(pan)
        self.arbre = ttk.Treeview(corps, show="headings", height=12)
        ay = ttk.Scrollbar(corps, orient="vertical", command=self.arbre.yview)
        ax = ttk.Scrollbar(corps, orient="horizontal", command=self.arbre.xview)
        self.arbre.configure(yscrollcommand=ay.set, xscrollcommand=ax.set)
        ax.pack(side="bottom", fill="x")
        ay.pack(side="right", fill="y")
        self.arbre.pack(side="left", fill="both", expand=True)
        self.arbre.bind("<Double-1>", self._editer_cellule)
        pan.add(corps, weight=3)
        bas = ttk.Frame(pan)
        self.desc = tk.Text(bas, wrap="word", height=8, state="disabled", padx=8, pady=6, font=tkfont.nametofont("TkFixedFont"), background="white", relief="solid", borderwidth=1)
        sy = ttk.Scrollbar(bas, orient="vertical", command=self.desc.yview)
        self.desc.configure(yscrollcommand=sy.set)
        sy.pack(side="right", fill="y")
        self.desc.pack(side="left", fill="both", expand=True)
        pan.add(bas, weight=1)

    def _ordre_lignes(self):
        idx = list(range(self.plan.n))
        if self.v_ordre.get() == "Ordre standard" and "Ordre_standard" in self.plan.cols:
            idx.sort(key=lambda i: self.plan.cols["Ordre_standard"][i])
        return idx

    def afficher_plan(self):
        a = self.arbre
        a.delete(*a.get_children())
        self.desc.configure(state="normal")
        self.desc.delete("1.0", "end")
        p = self.plan
        if p is None:
            a.configure(columns=())
            self.lbl_info.configure(text="Aucun plan.")
            self.desc.configure(state="disabled")
            self.cb_edit.configure(values=[])
            self.v_edit.set("")
            return
        cols = p.ordre_colonnes()
        a.configure(columns=cols)
        for c in cols:
            a.heading(c, text=c + (" ✎" if c in p.reponses else ""))
            larg = max([len(c) + 3] + [len(fmt_nombre(v, 6)) for v in p.cols[c][:200]])
            a.column(c, width=max(70, min(220, 10 * larg)), anchor="center", stretch=False)
        for i in self._ordre_lignes():
            a.insert("", "end", iid=str(i), values=[fmt_nombre(p.cols[c][i], 6) if p.cols[c][i] is not None else "" for c in cols])
        self.lbl_info.configure(text=f"{p.n} essais — double-clic sur une cellule ✎ pour saisir une réponse.")
        texte = [p.titre, "=" * len(p.titre), ""] + p.lignes
        for nom, t in p.tables.items():
            cl = list(t.keys())
            rows = [list(r) for r in zip(*[_liste(t[c]) for c in cl])]
            texte += ["", libelle(nom), table_texte(cl, rows, libre=True)]
        self.desc.insert("1.0", "\n".join(texte))
        self.desc.configure(state="disabled")
        self.cb_edit.configure(values=p.reponses)
        if self.v_edit.get() not in p.reponses:
            self.v_edit.set(p.reponses[0] if p.reponses else "")

    def _editer_cellule(self, ev):
        a = self.arbre
        if self.plan is None or a.identify_region(ev.x, ev.y) != "cell":
            return
        iid, col = a.identify_row(ev.y), a.identify_column(ev.x)
        nom = self.plan.ordre_colonnes()[int(col[1:]) - 1]
        if nom not in self.plan.reponses:
            return
        self._ouvrir_editeur(iid, col, nom)

    def _ouvrir_editeur(self, iid, col, nom):
        ttk = self.app.ttk
        a = self.arbre
        bbox = a.bbox(iid, col)
        if not bbox:
            return
        x, y, w, h = bbox
        e = ttk.Entry(a, justify="center")
        e.place(x=x, y=y, width=w, height=h)
        i = int(iid)
        v = self.plan.cols[nom][i]
        e.insert(0, "" if v is None else fmt_nombre(v, 10))
        e.select_range(0, "end")
        e.focus_set()
        fait = {"ok": False}

        def valider(suivant=False):
            if fait["ok"]:
                return
            fait["ok"] = True
            try:
                val = convertir_reponse(e.get())
            except ErreurDOE as err:
                from tkinter import messagebox
                e.destroy()
                messagebox.showerror(TITRE, str(err))
                return
            e.destroy()
            self.plan.cols[nom][i] = val
            cols = self.plan.ordre_colonnes()
            a.set(iid, nom, "" if val is None else fmt_nombre(val, 6))
            if suivant:
                enfants = a.get_children()
                k = enfants.index(iid)
                if k + 1 < len(enfants):
                    nxt = enfants[k + 1]
                    a.see(nxt)
                    a.selection_set(nxt)
                    a.after(20, lambda: self._ouvrir_editeur(nxt, col, nom))

        def annuler(_e=None):
            fait["ok"] = True
            e.destroy()
        e.bind("<Return>", lambda ev: valider(True))
        e.bind("<KP_Enter>", lambda ev: valider(True))
        e.bind("<Escape>", annuler)
        e.bind("<FocusOut>", lambda ev: valider(False))

    def _reponse_choisie(self):
        from tkinter import messagebox
        if self.plan is None:
            messagebox.showinfo(TITRE, "Aucun plan.")
            return None
        r = self.v_edit.get()
        if r not in self.plan.reponses:
            messagebox.showinfo(TITRE, "Choisir une réponse dans la liste.")
            return None
        return r

    def ajouter_reponse(self):
        from tkinter import simpledialog, messagebox
        if self.plan is None:
            messagebox.showinfo(TITRE, "Générer d'abord un plan.")
            return
        nom = simpledialog.askstring(TITRE, "Nom de la nouvelle réponse :", initialvalue=self.plan.nom_reponse_libre(), parent=self.app.root)
        if nom is None:
            return
        try:
            n = self.plan.ajouter_reponse(nom)
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self.v_edit.set(n)
        self.afficher_plan()
        self.maj_formulaires()

    def renommer_reponse(self):
        from tkinter import simpledialog, messagebox
        r = self._reponse_choisie()
        if r is None:
            return
        nom = simpledialog.askstring(TITRE, f"Nouveau nom de « {r} » :", initialvalue=r, parent=self.app.root)
        if nom is None:
            return
        try:
            self.plan.renommer_reponse(r, nom.strip())
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self.v_edit.set(nom.strip())
        self.afficher_plan()
        self.maj_formulaires()

    def supprimer_reponse(self):
        r = self._reponse_choisie()
        if r is None:
            return
        if len(self.plan.reponses) == 1 and boite_dialogue(self.app.root, TITRE, "C'est la seule réponse du plan. La supprimer quand même ?", ("Supprimer", "Annuler"), defaut=1) != 0:
            return
        self.plan.supprimer_reponse(r)
        self.v_edit.set("")
        self.afficher_plan()
        self.maj_formulaires()

    def coller(self):
        from tkinter import messagebox
        if self.plan is None:
            messagebox.showinfo(TITRE, "Générer d'abord un plan.")
            return
        DialogueCollage(self)

    def integrer_reponses(self, texte, ordre_standard):
        entete, vals = lire_tableau_texte(texte)
        p = self.plan
        if len(vals) != p.n:
            raise ErreurDOE(f"Le plan compte {p.n} essais, mais {len(vals)} ligne(s) de valeurs ont été fournies.")
        ncol = max(len(r) for r in vals)
        if entete:
            noms = [(entete[j] if j < len(entete) and entete[j] else None) for j in range(ncol)]
        else:
            noms = [p.reponses[j] if j < len(p.reponses) else None for j in range(ncol)]
        # indice du plan pour chaque ligne fournie
        if ordre_standard:
            pos = {int(s): i for i, s in enumerate(p.cols["Ordre_standard"])}
            cible = [pos.get(r + 1) for r in range(p.n)]
            if None in cible:
                raise ErreurDOE("L'ordre standard du plan est incomplet.")
        else:
            cible = list(range(p.n))
        existant = set(p.reponses)
        nouveaux = []
        for j in range(ncol):
            nom = noms[j] or p.nom_reponse_libre()
            if nom not in existant and nom not in nouveaux:
                p.verifier_nom_reponse(nom)
            noms[j] = nom
            nouveaux.append(nom)
        for j, nom in enumerate(noms):
            if nom not in p.cols:
                p.ajouter_reponse(nom)
            col = p.cols[nom]
            for r, i in enumerate(cible):
                col[i] = vals[r][j] if j < len(vals[r]) else None
        self.afficher_plan()
        self.maj_formulaires()
        self.app.statut(f"{ncol} colonne(s) de réponses enregistrée(s).")

    def exporter(self):
        from tkinter import messagebox
        if self.plan is None:
            messagebox.showinfo(TITRE, "Aucun plan à exporter.")
            return
        chemin = choisir_fichier_ecriture(self.app.root, title="Exporter le plan", defaultextension=".csv", filetypes=[("CSV", "*.csv")],
                                          initialdir=str(Path.home()), initialfile=f"plan_{self.cle}.csv")
        if chemin:
            self.plan.exporter_csv(chemin)
            self.app.statut(f"Plan exporté : {Path(chemin).name}")

    # ------------------------------------------------------------------ page 2 : analyse
    def _page_analyse(self, nb):
        tk, ttk = self.app.tk, self.app.ttk
        f = ttk.Frame(nb)
        nb.add(f, text="2. Analyse")
        cad = ttk.LabelFrame(f, text="Paramètres de l'analyse")
        cad.pack(fill="x", padx=6, pady=6)
        if self.cle == "taguchi":
            ttk.Label(cad, text="Réponses utilisées\n(sélection multiple)").grid(row=0, column=0, padx=6, pady=4, sticky="ne")
            self.lb_rep = tk.Listbox(cad, selectmode="extended", height=4, width=22, exportselection=False)
            self.lb_rep.grid(row=0, column=1, rowspan=3, padx=4, pady=4, sticky="w")
            ttk.Label(cad, text="Type de rapport signal/bruit").grid(row=0, column=2, padx=(14, 4), sticky="e")
            ttk.Combobox(cad, textvariable=self.v_rsb, state="readonly", width=38, values=[l for l, _ in TYPES_RSB]).grid(row=0, column=3, sticky="w")
            ttk.Checkbutton(cad, text="Regrouper automatiquement les termes non significatifs (erreur)", variable=self.v_regroup).grid(row=1, column=2, columnspan=2, sticky="w", padx=(14, 4))
            ttk.Label(cad, text="Essais de confirmation").grid(row=2, column=2, padx=(14, 4), sticky="e")
            ttk.Entry(cad, textvariable=self.v_nconf, width=6).grid(row=2, column=3, sticky="w")
        else:
            ttk.Label(cad, text="Réponse").grid(row=0, column=0, padx=6, pady=4, sticky="e")
            self.cb_rep_a = ttk.Combobox(cad, textvariable=self.v_rep, state="readonly", width=18)
            self.cb_rep_a.grid(row=0, column=1, sticky="w")
            ttk.Label(cad, text="Modèle").grid(row=0, column=2, padx=(14, 4), sticky="e")
            ttk.Combobox(cad, textvariable=self.v_modele, state="readonly", width=40, values=[l for l, _ in MODELES[self.cle]]).grid(row=0, column=3, sticky="w")
            ttk.Checkbutton(cad, text="Réduire le modèle (élimination des termes non significatifs)", variable=self.v_reduire).grid(row=1, column=0, columnspan=3, sticky="w", padx=6)
            ttk.Label(cad, text="Seuil de sortie (α)").grid(row=1, column=3, sticky="w", padx=(14, 0))
            ttk.Entry(cad, textvariable=self.v_aout, width=7).grid(row=1, column=3, sticky="e")
            if self.cle == "melange":
                ttk.Checkbutton(cad, text="Utiliser les pseudo-composants (bornes inférieures retranchées)", variable=self.v_pseudo).grid(row=2, column=0, columnspan=4, sticky="w", padx=6)
        ttk.Button(cad, text="Analyser", command=self.analyser).grid(row=3, column=0, padx=6, pady=6, sticky="w")
        self._zone(f, "analyse")

    def _params_modele(self):
        r = self.v_rep.get()
        if self.plan is None:
            raise ErreurDOE("Générer d'abord un plan.")
        if r not in self.plan.reponses:
            raise ErreurDOE("Choisir une réponse.")
        if self.plan.nb_renseignes(r) < 2:
            raise ErreurDOE(f"La réponse « {r} » doit être renseignée pour au moins deux essais.")
        p = {"reponse": r, "modele": code_de(self.v_modele.get(), MODELES[self.cle]), "reduire": bool(self.v_reduire.get()),
             "alpha_sortie": lire_alpha(self.v_aout.get())}
        if self.cle == "melange":
            p["pseudo"] = bool(self.v_pseudo.get())
        return p

    def _params_taguchi(self):
        if self.plan is None:
            raise ErreurDOE("Générer d'abord un plan.")
        reps = [self.lb_rep.get(i) for i in self.lb_rep.curselection()]
        if not reps:
            raise ErreurDOE("Sélectionner au moins une réponse.")
        for r in reps:
            if self.plan.nb_renseignes(r) < self.plan.n:
                raise ErreurDOE(f"La réponse « {r} » n'est pas renseignée pour tous les essais.")
        nc = lire_nombre(self.v_nconf.get(), "le nombre d'essais de confirmation")
        if nc is None or nc < 1:
            raise ErreurDOE("Le nombre d'essais de confirmation doit être au moins 1.")
        return {"reponses": reps, "type_rsb": code_de(self.v_rsb.get(), TYPES_RSB), "regroupement": "auto" if self.v_regroup.get() else "aucun",
                "n_confirmation": nc}

    def analyser(self):
        from tkinter import messagebox
        try:
            p = self._params_taguchi() if self.cle == "taguchi" else self._params_modele()
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self.app.lancer(self.zones["analyse"], "analyse", p)

    # ------------------------------------------------------------------ page 3 : surfaces
    def _page_surfaces(self, nb):
        tk, ttk = self.app.tk, self.app.ttk
        f = ttk.Frame(nb)
        nb.add(f, text="3. Graphiques du modèle")
        cad = ttk.LabelFrame(f, text="Représentation du modèle ajusté (utilise le modèle choisi dans l'onglet Analyse)")
        cad.pack(fill="x", padx=6, pady=6)
        ttk.Label(cad, text="Réponse").grid(row=0, column=0, padx=6, pady=4, sticky="e")
        self.cb_rep_s = ttk.Combobox(cad, textvariable=self.v_rep, state="readonly", width=18)
        self.cb_rep_s.grid(row=0, column=1, sticky="w")
        self.frm_axes = ttk.Frame(cad)
        self.frm_axes.grid(row=1, column=0, columnspan=6, sticky="ew", padx=6, pady=4)
        ttk.Button(cad, text="Tracer", command=self.tracer).grid(row=2, column=0, padx=6, pady=6, sticky="w")
        self._zone(f, "surfaces")

    def _maj_axes(self):
        if self.cle not in ("factoriel", "surface", "melange"):
            return
        ttk = self.app.ttk
        for w in self.frm_axes.winfo_children():
            w.destroy()
        p = self.plan
        if p is None:
            return
        noms = p.noms_facteurs
        self.fixes_vars = {}
        if self.cle == "melange":
            ttk.Label(self.frm_axes, text="Composants du triangle").grid(row=0, column=0, sticky="e", padx=4)
            for j, v in enumerate((self.v_ca, self.v_cb, self.v_cc)):
                if v.get() not in noms:
                    v.set(noms[j] if j < len(noms) else "")
                ttk.Combobox(self.frm_axes, textvariable=v, state="readonly", width=14, values=noms).grid(row=0, column=1 + j, padx=3)
            if len(noms) > 3:
                ttk.Label(self.frm_axes, text="Proportions de référence des autres composants (facultatif) :").grid(row=1, column=0, columnspan=4, sticky="w", pady=(6, 2))
                for j, n in enumerate(noms):
                    self.fixes_vars[n] = self.app.tk.StringVar()
                    ttk.Label(self.frm_axes, text=n).grid(row=2 + j // 4 * 2, column=j % 4, sticky="w", padx=4)
                    ttk.Entry(self.frm_axes, textvariable=self.fixes_vars[n], width=10).grid(row=3 + j // 4 * 2, column=j % 4, sticky="w", padx=4)
            return
        quant = [f["nom"] for f in p.facteurs if f.get("type", "num") == "num"]
        if self.v_ax.get() not in quant:
            self.v_ax.set(quant[0] if quant else "")
        if self.v_ay.get() not in quant or self.v_ay.get() == self.v_ax.get():
            self.v_ay.set(quant[1] if len(quant) > 1 else "")
        ttk.Label(self.frm_axes, text="Axe horizontal").grid(row=0, column=0, sticky="e", padx=4)
        ttk.Combobox(self.frm_axes, textvariable=self.v_ax, state="readonly", width=14, values=quant).grid(row=0, column=1, padx=3)
        ttk.Label(self.frm_axes, text="Axe vertical").grid(row=0, column=2, sticky="e", padx=4)
        ttk.Combobox(self.frm_axes, textvariable=self.v_ay, state="readonly", width=14, values=quant).grid(row=0, column=3, padx=3)
        ttk.Label(self.frm_axes, text="Type").grid(row=0, column=4, sticky="e", padx=4)
        ttk.Combobox(self.frm_axes, textvariable=self.v_graphe, state="readonly", width=18, values=[l for l, _ in TYPES_GRAPHE]).grid(row=0, column=5, padx=3)
        if len(noms) > 2:
            ttk.Label(self.frm_axes, text="Valeurs des autres facteurs (vide = centre du domaine) :").grid(row=1, column=0, columnspan=6, sticky="w", pady=(6, 2))
            for j, f in enumerate(p.facteurs):
                n = f["nom"]
                self.fixes_vars[n] = self.app.tk.StringVar()
                ttk.Label(self.frm_axes, text=n).grid(row=2 + j // 4 * 2, column=j % 4 * 1 if False else j % 4, sticky="w", padx=4)
                if f.get("type") == "cat":
                    w = ttk.Combobox(self.frm_axes, textvariable=self.fixes_vars[n], state="readonly", width=12, values=[""] + list(f["niveaux"]))
                else:
                    w = ttk.Entry(self.frm_axes, textvariable=self.fixes_vars[n], width=12)
                w.grid(row=3 + j // 4 * 2, column=j % 4, sticky="w", padx=4)

    def tracer(self):
        from tkinter import messagebox
        try:
            p = self._params_modele()
            fixes = {}
            for n, v in self.fixes_vars.items():
                t = v.get().strip()
                if not t:
                    continue
                f = next(f for f in self.plan.facteurs if f["nom"] == n)
                fixes[n] = t if f.get("type") == "cat" else lire_nombre(t, f"la valeur de « {n} »")
            p["fixes"] = fixes
            if self.cle == "melange":
                if len(self.plan.facteurs) < 3:
                    raise ErreurDOE("Le diagramme ternaire exige au moins trois composants.")
                noms = self.plan.noms_facteurs
                ch = [self.v_ca.get(), self.v_cb.get(), self.v_cc.get()]
                if len(set(ch)) < 3 or any(c not in noms for c in ch):
                    raise ErreurDOE("Choisir trois composants différents pour le triangle.")
                p["composants_axes"] = [noms.index(c) + 1 for c in ch]
            else:
                if not self.v_ax.get() or not self.v_ay.get() or self.v_ax.get() == self.v_ay.get():
                    raise ErreurDOE("Choisir deux facteurs quantitatifs différents pour les axes.")
                p.update({"axe_x": self.v_ax.get(), "axe_y": self.v_ay.get(), "type": code_de(self.v_graphe.get(), TYPES_GRAPHE)})
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self.app.lancer(self.zones["surfaces"], "surfaces", p)

    # ------------------------------------------------------------------ page 4 : prévision et optimisation
    def _page_prevision(self, nb):
        tk, ttk = self.app.tk, self.app.ttk
        f = ttk.Frame(nb)
        nb.add(f, text="4. Prévision" if self.cle == "taguchi" else "4. Prévision et optimisation")
        if self.cle == "taguchi":
            self._bloc_prevision(f)
            return
        sous = ttk.Notebook(f)
        sous.pack(fill="both", expand=True, padx=2, pady=2)
        a, b = ttk.Frame(sous), ttk.Frame(sous)
        sous.add(a, text="Prévision")
        sous.add(b, text="Optimisation")
        self._bloc_prevision(a)
        self._bloc_optimisation(b)

    def _bloc_prevision(self, f):
        tk, ttk = self.app.tk, self.app.ttk
        cad = ttk.LabelFrame(f, text="Points à prévoir")
        cad.pack(fill="x", padx=6, pady=6)
        if self.cle != "taguchi":
            ttk.Label(cad, text="Réponse").grid(row=0, column=0, padx=6, pady=4, sticky="e")
            self.cb_rep_p = ttk.Combobox(cad, textvariable=self.v_rep, state="readonly", width=18)
            self.cb_rep_p.grid(row=0, column=1, sticky="w")
        self.frm_prev = ttk.Frame(cad)
        self.frm_prev.grid(row=1, column=0, columnspan=6, sticky="ew", padx=6, pady=2)
        barre = ttk.Frame(cad)
        barre.grid(row=2, column=0, columnspan=6, sticky="w", padx=6, pady=4)
        ttk.Button(barre, text="Ajouter le point", command=self.ajouter_point).pack(side="left")
        ttk.Button(barre, text="Retirer", command=self.retirer_point).pack(side="left", padx=4)
        ttk.Button(barre, text="Valeurs centrales" if self.cle != "taguchi" else "Premiers niveaux", command=self.valeurs_centrales).pack(side="left", padx=4)
        ttk.Button(barre, text="Utiliser l'optimum", command=self.utiliser_optimum).pack(side="left", padx=4)
        ttk.Button(barre, text="Prévoir", command=self.prevoir).pack(side="left", padx=12)
        self.liste_pts = tk.Listbox(cad, height=4, exportselection=False)
        self.liste_pts.grid(row=3, column=0, columnspan=6, sticky="ew", padx=6, pady=(0, 6))
        cad.columnconfigure(5, weight=1)
        self.points = []
        self._zone(f, "prevision")

    def _bloc_optimisation(self, f):
        tk, ttk = self.app.tk, self.app.ttk
        cad = ttk.LabelFrame(f, text="Objectifs (fonction de désirabilité de Derringer et Suich)")
        cad.pack(fill="x", padx=6, pady=6)
        self.frm_opt = ttk.Frame(cad)
        self.frm_opt.pack(fill="x", padx=6, pady=4)
        ttk.Label(cad, text="Valeurs de la cible et des limites : vide = déduit des valeurs observées.", foreground="#666666").pack(anchor="w", padx=6)
        ttk.Button(cad, text="Optimiser", command=self.optimiser).pack(anchor="w", padx=6, pady=6)
        self._zone(f, "optimisation")

    def _valeur_defaut(self, f):
        if self.cle == "taguchi":
            return str(f["niveaux"][0])
        if f.get("type") == "cat":
            return str(f["niveaux"][0])
        if "bas" in f and "haut" in f:
            if self.cle == "melange":
                col = [v for v in self.plan.cols[f["nom"]] if v is not None]
                return fmt_nombre(sum(col) / len(col), 6) if col else ""
            return fmt_nombre((f["bas"] + f["haut"]) / 2, 6)
        return ""

    def _maj_prevision(self):
        ttk = self.app.ttk
        for w in self.frm_prev.winfo_children():
            w.destroy()
        self.prev_vars = {}
        self.points = []
        self.liste_pts.delete(0, "end")
        p = self.plan
        if p is None:
            return
        for j, f in enumerate(p.facteurs):
            n = f["nom"]
            v = self.app.tk.StringVar(value=self._valeur_defaut(f))
            self.prev_vars[n] = v
            ttk.Label(self.frm_prev, text=n).grid(row=j // 4 * 2, column=j % 4, sticky="w", padx=4)
            if self.cle == "taguchi" or f.get("type") == "cat":
                w = ttk.Combobox(self.frm_prev, textvariable=v, state="readonly", width=12, values=[str(x) for x in f["niveaux"]])
            else:
                w = ttk.Entry(self.frm_prev, textvariable=v, width=12)
            w.grid(row=j // 4 * 2 + 1, column=j % 4, sticky="w", padx=4, pady=(0, 4))

    def valeurs_centrales(self):
        if self.plan is None:
            return
        for f in self.plan.facteurs:
            self.prev_vars[f["nom"]].set(self._valeur_defaut(f))

    def _lire_point(self):
        pt = {}
        for f in self.plan.facteurs:
            n = f["nom"]
            t = self.prev_vars[n].get().strip()
            if not t:
                raise ErreurDOE(f"Saisir une valeur pour « {n} ».")
            if self.cle == "taguchi" or f.get("type") == "cat":
                pt[n] = t
            else:
                pt[n] = lire_nombre(t, f"la valeur de « {n} »")
        return pt

    def ajouter_point(self):
        from tkinter import messagebox
        if self.plan is None:
            messagebox.showinfo(TITRE, "Générer d'abord un plan.")
            return
        try:
            pt = self._lire_point()
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self.points.append(pt)
        self.liste_pts.insert("end", " ; ".join(f"{k} = {fmt_nombre(v, 6)}" for k, v in pt.items()))

    def retirer_point(self):
        sel = self.liste_pts.curselection()
        if sel:
            del self.points[sel[0]]
            self.liste_pts.delete(sel[0])

    def utiliser_optimum(self):
        from tkinter import messagebox
        o = self.dernier_optimum
        if not o:
            messagebox.showinfo(TITRE, "Aucun optimum disponible : lancer d'abord " + ("l'analyse." if self.cle == "taguchi" else "l'optimisation."))
            return
        for n, v in o.items():
            if n in self.prev_vars:
                self.prev_vars[n].set(fmt_nombre(v, 8) if not isinstance(v, str) else v)

    def prevoir(self):
        from tkinter import messagebox
        try:
            if self.plan is None:
                raise ErreurDOE("Générer d'abord un plan.")
            pts = list(self.points) or [self._lire_point()]
            p = self._params_taguchi() if self.cle == "taguchi" else self._params_modele()
            p["points"] = pts
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self.app.lancer(self.zones["prevision"], "prevision", p)

    def _maj_optim(self):
        if self.cle == "taguchi":
            return
        tk, ttk = self.app.tk, self.app.ttk
        for w in self.frm_opt.winfo_children():
            w.destroy()
        old = self.optim_lignes
        self.optim_lignes = {}
        if self.plan is None:
            return
        for j, t in enumerate(["Utiliser", "Réponse", "But", "Limite inférieure", "Cible", "Limite supérieure", "Importance"]):
            ttk.Label(self.frm_opt, text=t, font=("TkDefaultFont", 9, "bold")).grid(row=0, column=j, padx=4, sticky="w")
        for i, r in enumerate(self.plan.reponses):
            o = old.get(r) or {"use": tk.BooleanVar(value=True), "but": tk.StringVar(value=BUTS[0][0]), "inf": tk.StringVar(), "cible": tk.StringVar(),
                               "sup": tk.StringVar(), "imp": tk.StringVar(value="1")}
            self.optim_lignes[r] = o
            ttk.Checkbutton(self.frm_opt, variable=o["use"]).grid(row=i + 1, column=0)
            ttk.Label(self.frm_opt, text=r).grid(row=i + 1, column=1, sticky="w", padx=4)
            ttk.Combobox(self.frm_opt, textvariable=o["but"], state="readonly", width=18, values=[l for l, _ in BUTS]).grid(row=i + 1, column=2, padx=4, pady=1)
            for c, k in ((3, "inf"), (4, "cible"), (5, "sup"), (6, "imp")):
                ttk.Entry(self.frm_opt, textvariable=o[k], width=11).grid(row=i + 1, column=c, padx=4)

    def optimiser(self):
        from tkinter import messagebox
        try:
            if self.plan is None:
                raise ErreurDOE("Générer d'abord un plan.")
            base = self._params_modele_sans_reponse()
            reps = []
            for r, o in self.optim_lignes.items():
                if not o["use"].get():
                    continue
                if self.plan.nb_renseignes(r) < 2:
                    raise ErreurDOE(f"La réponse « {r} » doit être renseignée pour au moins deux essais.")
                d = {"nom": r, "but": code_de(o["but"].get(), BUTS)}
                for k in ("inf", "cible", "sup"):
                    v = lire_nombre(o[k].get(), f"{k} de « {r} »")
                    if v is not None:
                        d[k] = v
                imp = lire_nombre(o["imp"].get() or "1", f"l'importance de « {r} »")
                if imp is None or imp <= 0:
                    raise ErreurDOE(f"L'importance de « {r} » doit être strictement positive.")
                d["importance"] = imp
                reps.append(d)
            if not reps:
                raise ErreurDOE("Cocher au moins une réponse à optimiser.")
            base["reponses"] = reps
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self.app.lancer(self.zones["optimisation"], "optimisation", base)

    def _params_modele_sans_reponse(self):
        p = {"modele": code_de(self.v_modele.get(), MODELES[self.cle]), "reduire": bool(self.v_reduire.get()), "alpha_sortie": lire_alpha(self.v_aout.get())}
        if self.cle == "melange":
            p["pseudo"] = bool(self.v_pseudo.get())
        return p

    # ------------------------------------------------------------------ synchronisation des formulaires
    def maj_formulaires(self):
        p = self.plan
        reps = list(p.reponses) if p else []
        for nom in ("cb_rep_a", "cb_rep_s", "cb_rep_p"):
            w = getattr(self, nom, None)
            if w is not None:
                w.configure(values=reps)
        if self.cle == "taguchi":
            sel = [self.lb_rep.get(i) for i in self.lb_rep.curselection()]
            self.lb_rep.delete(0, "end")
            for r in reps:
                self.lb_rep.insert("end", r)
            for i, r in enumerate(reps):
                if r in sel or not sel:
                    self.lb_rep.selection_set(i)
        elif self.v_rep.get() not in reps:
            self.v_rep.set(reps[0] if reps else "")
        self._maj_axes()
        self._maj_prevision()
        self._maj_optim()
        for z in self.zones.values():
            z.reinitialiser()

    # ------------------------------------------------------------------ projet
    def etat_pour_projet(self):
        return self.plan

    def restaurer(self, plan):
        self.plan = plan
        if plan is not None:
            self.v_nb.set(str(len(plan.facteurs)))
            self._maj_nb()
            for v, f in zip(self.lignes_fact, plan.facteurs):
                v["nom"].set(f["nom"])
                if self.cle == "taguchi":
                    v["niv"].set(" ; ".join(str(x) for x in f["niveaux"]))
                elif self.cle == "melange":
                    v["a"].set(fmt_nombre(f["bas"], 8))
                    v["b"].set(fmt_nombre(f["haut"], 8))
                elif f.get("type") == "cat":
                    v["a"].set(str(f["niveaux"][0]))
                    v["b"].set(str(f["niveaux"][1]))
                else:
                    v["a"].set(fmt_nombre(f["bas"], 8))
                    v["b"].set(fmt_nombre(f["haut"], 8))
        self.dernier_optimum = None
        self.afficher_plan()
        self.maj_formulaires()


# =============================================================================
#  Aide
# =============================================================================
TEXTE_AIDE = """\
DÉMARCHE GÉNÉRALE
Chaque famille de plans suit les mêmes étapes, dans les sous-onglets :
  1. Plan : définir les facteurs, choisir les options, générer la matrice d'essais, puis saisir les réponses mesurées
     (double-clic dans une cellule marquée ✎, ou « Saisir ou coller des réponses… »).
  2. Analyse : choisir la réponse et le modèle ; le programme fournit l'ANOVA, les effets, les graphiques,
     les équations et le diagnostic des résidus.
  3. Graphiques du modèle (factoriel, surface, mélange) : lignes de niveau, surface 3D ou diagramme ternaire.
  4. Prévision et optimisation : prévision avec intervalles, puis recherche des réglages optimaux.

FAMILLES
  Criblage : repérer les facteurs influents avec un minimum d'essais (Plackett-Burman ou fractions à 2 niveaux).
  Factoriel : plans complets ou fractionnaires à 2 niveaux (générateurs, résolution et structure d'alias affichées),
     analyse des effets et interactions, réduction du modèle, prévision et optimisation.
  Taguchi : tableaux orthogonaux, rapport signal/bruit, effets des niveaux, contributions, niveaux optimaux
     et intervalle de confirmation.
  Surface : plans composites centrés et Box-Behnken, modèle du second degré, point stationnaire, optimisation.
  Mélange : composants dont la somme vaut 1, modèles de Scheffé, diagramme ternaire, optimisation.

SAISIE DES NOMBRES
  La virgule et le point sont acceptés comme séparateur décimal. Dans un plan de Taguchi, les niveaux d'un facteur
  sont séparés par des points-virgules.

PLAN ALÉATOIRE
  L'ordre d'exécution est aléatoire par défaut. Une valeur d'initialisation permet de reproduire le même ordre.

RAPPORT ET PROJET
  Chaque résultat calculé est ajouté au rapport. L'onglet Rapport permet de choisir les résultats à inclure et
  d'exporter un document HTML. Le projet complet (plans, réponses et résultats) s'enregistre dans un fichier
  « .doeproj » et peut être rouvert plus tard.
"""


class EcranAccueil:
    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk
        self.win = tk.Toplevel(root)
        self.win.overrideredirect(True)
        cadre = tk.Frame(self.win, background="white", highlightbackground="#1f4e8c", highlightthickness=2)
        cadre.pack(fill="both", expand=True)
        self.img = None
        try:
            if ICONE.is_file():
                self.img = charger_image(ICONE, 140, 140)
                tk.Label(cadre, image=self.img, background="white").pack(padx=40, pady=(26, 6))
        except Exception:
            pass
        tk.Label(cadre, text="DOE", font=("TkDefaultFont", 22, "bold"), foreground="#1f4e8c", background="white").pack()
        tk.Label(cadre, text="Plans d'expériences", font=("TkDefaultFont", 11), background="white").pack(pady=(0, 6))
        self.lbl = tk.Label(cadre, text="Démarrage…", background="white", foreground="#555555")
        self.lbl.pack(pady=(4, 20))
        self.win.update_idletasks()
        w, h = self.win.winfo_reqwidth(), self.win.winfo_reqheight()
        self.win.geometry(f"{w}x{h}+{(self.win.winfo_screenwidth() - w) // 2}+{(self.win.winfo_screenheight() - h) // 3}")
        self.win.update()

    def message(self, t):
        self.lbl.configure(text=t)
        self.win.update()

    def fermer(self):
        self.win.destroy()


# =============================================================================
#  Application
# =============================================================================
class Application:
    def __init__(self):
        import tkinter as tk
        from tkinter import ttk
        self.tk, self.ttk = tk, ttk
        self.root = tk.Tk()
        self.root.withdraw()
        self.root.title(TITRE)
        self.zones = []
        self.journal = []
        self.file = queue.Queue()
        self.occupe = False
        self.chemin_projet = None
        self.v_alpha = tk.StringVar(value="0,05")
        accueil = EcranAccueil(self.root)
        ok, msg = verifier_r()
        accueil.message(msg if not ok else "Préparation de l'interface…")
        try:
            style = ttk.Style()
            if "clam" in style.theme_names():
                style.theme_use("clam")
        except Exception:
            pass
        try:
            if ICONE.is_file():
                self._icone = charger_image(ICONE, 64, 64)
                self.root.iconphoto(True, self._icone)
        except Exception:
            pass
        self._construire()
        self.root.update_idletasks()
        w = min(1280, self.root.winfo_screenwidth() - 60)
        h = min(820, self.root.winfo_screenheight() - 90)
        self.root.geometry(f"{w}x{h}+20+20")
        self.root.minsize(900, 560)
        accueil.fermer()
        self.root.deiconify()
        self.root.protocol("WM_DELETE_WINDOW", self.quitter)
        self.root.after(100, self._surveiller)
        if not ok:
            from tkinter import messagebox
            messagebox.showwarning(TITRE, msg)
        self.statut(msg if not ok else "Prêt.")

    # ------------------------------------------------------------------ construction
    def _construire(self):
        tk, ttk = self.tk, self.ttk
        haut = ttk.Frame(self.root)
        haut.pack(fill="x", padx=8, pady=(6, 0))
        ttk.Button(haut, text="Nouveau projet", command=self.nouveau).pack(side="left")
        ttk.Button(haut, text="Ouvrir un projet…", command=self.ouvrir).pack(side="left", padx=4)
        ttk.Button(haut, text="Enregistrer le projet…", command=self.enregistrer).pack(side="left")
        ttk.Label(haut, text="Risque α des tests et intervalles :").pack(side="left", padx=(24, 4))
        ttk.Entry(haut, textvariable=self.v_alpha, width=7).pack(side="left")
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill="both", expand=True, padx=6, pady=6)
        self.familles = {}
        for cle in FAMILLES:
            fam = Famille(self, cle)
            cadre = ttk.Frame(self.nb)
            self.nb.add(cadre, text=FAMILLES[cle]["titre"])
            fam.construire(cadre).pack(fill="both", expand=True)
            self.familles[cle] = fam
        self._page_rapport()
        aide = ttk.Frame(self.nb)
        self.nb.add(aide, text="Aide")
        t = tk.Text(aide, wrap="word", padx=14, pady=10, background="white", relief="flat")
        t.insert("1.0", TEXTE_AIDE)
        t.configure(state="disabled")
        t.pack(fill="both", expand=True)
        bas = ttk.Frame(self.root)
        bas.pack(fill="x", side="bottom")
        self.v_statut = tk.StringVar()
        ttk.Label(bas, textvariable=self.v_statut, anchor="w").pack(side="left", fill="x", expand=True, padx=8, pady=2)
        self.barre = ttk.Progressbar(bas, mode="indeterminate", length=120)

    def statut(self, texte):
        self.v_statut.set(texte)

    # ------------------------------------------------------------------ calculs en arrière-plan
    def _alpha(self):
        a = lire_alpha(self.v_alpha.get())
        return a

    def _occuper(self, oui, texte=""):
        self.occupe = oui
        if oui:
            self.barre.pack(side="right", padx=8)
            self.barre.start(12)
            self.statut(texte)
        else:
            self.barre.stop()
            self.barre.pack_forget()

    def lancer_plan(self, fam, commande, params, facteurs):
        if self.occupe:
            self.statut("Un calcul est déjà en cours.")
            return
        self._occuper(True, "Construction du plan…")

        def travail():
            try:
                plan = generer_plan(fam.cle, commande, params, facteurs)
                self.file.put(("plan", fam, plan, None))
            except Exception as e:
                self.file.put(("plan", fam, None, e))
        threading.Thread(target=travail, daemon=True).start()

    def lancer(self, zone, commande, params):
        from tkinter import messagebox
        if self.occupe:
            self.statut("Un calcul est déjà en cours.")
            return
        try:
            alpha = self._alpha()
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        fam = zone.famille
        if fam.plan is None:
            messagebox.showinfo(TITRE, "Générer d'abord un plan.")
            return
        ech = echelle_ecran(self.root)
        self._occuper(True, "Calcul en cours…")
        zone.vider()

        def travail():
            try:
                res = executer_r(fam.plan, commande, params, ech, alpha)
                self.file.put(("res", zone, res, None))
            except Exception as e:
                self.file.put(("res", zone, None, e))
        threading.Thread(target=travail, daemon=True).start()

    def _surveiller(self):
        from tkinter import messagebox
        try:
            while True:
                msg = self.file.get_nowait()
                if msg[0] == "statut":
                    self.statut(msg[1])
                elif msg[0] == "plan":
                    _, fam, plan, err = msg
                    self._occuper(False)
                    if err is not None:
                        self.statut("Échec de la construction du plan.")
                        messagebox.showerror(TITRE, str(err) if isinstance(err, ErreurDOE) else f"Erreur inattendue : {err}")
                    else:
                        fam.recevoir_plan(plan)
                elif msg[0] == "res":
                    _, zone, res, err = msg
                    self._occuper(False)
                    if err is not None:
                        zone.afficher_erreur(str(err) if isinstance(err, ErreurDOE) else f"Erreur inattendue : {err}")
                        self.statut("Échec du calcul.")
                    else:
                        zone.afficher(res)
                        if res.optimum:
                            zone.famille.dernier_optimum = res.optimum
                        self.ajouter_journal(zone.famille, zone, res)
                        self.statut(res.conclusion or "Calcul terminé.")
        except queue.Empty:
            pass
        self.root.after(100, self._surveiller)

    # ------------------------------------------------------------------ rapport
    def _page_rapport(self):
        tk, ttk = self.tk, self.ttk
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="Rapport")
        ttk.Label(f, wraplength=900, justify="left",
                  text="Chaque résultat calculé est ajouté ici. Cocher (✓) ceux à inclure dans le rapport HTML ; double-clic pour revoir un résultat.").pack(anchor="w", padx=8, pady=(8, 4))
        barre = ttk.Frame(f)
        barre.pack(fill="x", padx=8, pady=4)
        for t, c in (("Inclure / exclure", self.basculer), ("Revoir", self.revoir), ("Supprimer", self.supprimer_entree), ("Tout vider", self.vider_journal)):
            ttk.Button(barre, text=t, command=c).pack(side="left", padx=(0, 4))
        ttk.Button(barre, text="Exporter le rapport (HTML)…", command=self.exporter_rapport).pack(side="right")
        self.arbre_j = ttk.Treeview(f, columns=("inclus", "heure", "famille", "titre"), show="headings", height=14, selectmode="browse")
        for c, t, w in (("inclus", "Inclus", 60), ("heure", "Date et heure", 150), ("famille", "Famille", 120), ("titre", "Résultat", 600)):
            self.arbre_j.heading(c, text=t)
            self.arbre_j.column(c, width=w, anchor="center" if c == "inclus" else "w")
        self.arbre_j.pack(fill="both", expand=True, padx=8, pady=6)
        self.arbre_j.bind("<Double-1>", lambda e: self.revoir())
        self.arbre_j.bind("<space>", lambda e: self.basculer())

    def ajouter_journal(self, fam, zone, res):
        e = EntreeJournal(fam, zone, res)
        self.journal.append(e)
        self._rafraichir_journal()

    def _rafraichir_journal(self):
        a = self.arbre_j
        a.delete(*a.get_children())
        for i, e in enumerate(self.journal):
            a.insert("", "end", iid=str(i), values=("✓" if e.inclus else "", e.heure.replace("T", " "), FAMILLES[e.famille.cle]["titre"][3:], e.res.titre))

    def _entree_choisie(self):
        from tkinter import messagebox
        sel = self.arbre_j.selection()
        if not sel:
            messagebox.showinfo(TITRE, "Sélectionner un résultat dans la liste.")
            return None
        return int(sel[0])

    def basculer(self):
        i = self._entree_choisie()
        if i is not None:
            self.journal[i].inclus = not self.journal[i].inclus
            self._rafraichir_journal()
            self.arbre_j.selection_set(str(i))

    def revoir(self):
        i = self._entree_choisie()
        if i is None:
            return
        e = self.journal[i]
        win = self.tk.Toplevel(self.root)
        win.title(f"{TITRE} — {e.res.titre}")
        win.geometry("1100x680")
        z = ZoneResultat(win, self, e.famille, "revue")
        z.frame.pack(fill="both", expand=True, padx=6, pady=6)
        z.afficher(e.res)

    def supprimer_entree(self):
        i = self._entree_choisie()
        if i is not None:
            del self.journal[i]
            self._rafraichir_journal()

    def vider_journal(self):
        if self.journal and boite_dialogue(self.root, TITRE, "Supprimer tous les résultats du rapport ?", ("Supprimer", "Annuler"), defaut=1) == 0:
            self.journal.clear()
            self._rafraichir_journal()

    def exporter_rapport(self):
        from tkinter import messagebox
        inclus = [e for e in self.journal if e.inclus]
        plans = [(c, f.plan) for c, f in self.familles.items() if f.plan is not None]
        if not inclus and not plans:
            messagebox.showinfo(TITRE, "Rien à inclure dans le rapport.")
            return
        chemin = choisir_fichier_ecriture(self.root, title="Exporter le rapport", defaultextension=".html", filetypes=[("Page HTML", "*.html")],
                                          initialdir=str(Path.home()), initialfile="rapport_doe.html")
        if not chemin:
            return
        Path(chemin).write_text(self._html(inclus, plans), encoding="utf-8")
        self.statut(f"Rapport exporté : {Path(chemin).name}")

    def _html(self, inclus, plans):
        esc = html.escape
        p = ["<!DOCTYPE html><html lang='fr'><head><meta charset='utf-8'><title>Rapport — plans d'expériences</title>",
             "<style>body{font-family:Arial,sans-serif;max-width:1000px;margin:24px auto;padding:0 16px;color:#222}"
             "h1{color:#1f4e8c}h2{color:#1f4e8c;border-bottom:2px solid #1f4e8c;padding-bottom:4px;margin-top:36px}"
             "pre{background:#f6f8fa;padding:10px;overflow-x:auto;font-size:12.5px}img{max-width:100%;border:1px solid #ddd;margin:8px 0}"
             "table{border-collapse:collapse;font-size:12.5px}td,th{border:1px solid #ccc;padding:3px 8px;text-align:center}th{background:#eef3fa}"
             ".rem{color:#8a5a00}</style></head><body>",
             "<h1>Rapport — plans d'expériences</h1>",
             f"<p>Établi le {datetime.datetime.now().strftime('%d/%m/%Y à %H:%M')}.</p>"]
        for cle, plan in plans:
            p.append(f"<h2>Plan — {esc(FAMILLES[cle]['titre'][3:])}</h2><p>{esc(plan.titre)} — {plan.n} essais.</p>")
            cols = plan.ordre_colonnes()
            p.append("<table><tr>" + "".join(f"<th>{esc(c)}</th>" for c in cols) + "</tr>")
            for i in range(plan.n):
                p.append("<tr>" + "".join(f"<td>{esc(fmt_nombre(plan.cols[c][i], 6)) if plan.cols[c][i] is not None else ''}</td>" for c in cols) + "</tr>")
            p.append("</table>")
            if plan.lignes:
                p.append("<pre>" + esc("\n".join(plan.lignes)) + "</pre>")
        for e in inclus:
            r = e.res
            p.append(f"<h2>{esc(r.titre)}</h2><p><small>{esc(FAMILLES[e.famille.cle]['titre'][3:])} — {esc(e.heure.replace('T', ' '))}</small></p>")
            p.append("<pre>" + esc(r.texte(titres=False)) + "</pre>")
            for f in r.figures:
                try:
                    b = base64.b64encode(f.read_bytes()).decode("ascii")
                    p.append(f"<img src='data:image/png;base64,{b}' alt='{esc(f.name)}'>")
                except OSError:
                    pass
        p.append("</body></html>")
        return "\n".join(p)

    # ------------------------------------------------------------------ projet
    def _entrees(self):
        return [(e.famille.cle, e.heure, e.res, e.inclus) for e in self.journal]

    def enregistrer(self):
        from tkinter import messagebox
        plans = {c: f.plan for c, f in self.familles.items() if f.plan is not None}
        if not plans and not self.journal:
            messagebox.showinfo(TITRE, "Rien à enregistrer.")
            return
        chemin = choisir_fichier_ecriture(self.root, title="Enregistrer le projet", defaultextension=EXT_PROJET,
                                          filetypes=[("Projet DOE", "*" + EXT_PROJET)], initialdir=str(Path.home()),
                                          initialfile=Path(self.chemin_projet).name if self.chemin_projet else "projet" + EXT_PROJET)
        if not chemin:
            return
        try:
            enregistrer_projet(chemin, plans, self._entrees(), {"alpha": self.v_alpha.get()})
        except Exception as e:
            messagebox.showerror(TITRE, f"Enregistrement impossible : {e}")
            return
        self.chemin_projet = chemin
        self.statut(f"Projet enregistré : {Path(chemin).name}")

    def ouvrir(self):
        from tkinter import filedialog, messagebox
        chemin = filedialog.askopenfilename(parent=self.root, title="Ouvrir un projet", filetypes=[("Projet DOE", "*" + EXT_PROJET), ("Tous les fichiers", "*.*")])
        if not chemin:
            return
        try:
            meta, plans, entrees = lire_projet(chemin)
        except ErreurDOE as e:
            messagebox.showerror(TITRE, str(e))
            return
        self._vider_tout()
        if meta.get("alpha"):
            self.v_alpha.set(str(meta["alpha"]))
        for cle, plan in plans.items():
            if cle in self.familles:
                self.familles[cle].restaurer(plan)
        for cle, heure, res, inclus in entrees:
            if cle in self.familles:
                self.journal.append(EntreeJournal(self.familles[cle], None, res, heure, inclus))
        self._rafraichir_journal()
        self.chemin_projet = chemin
        self.statut(f"Projet ouvert : {Path(chemin).name}")

    def _vider_tout(self):
        for f in self.familles.values():
            f.plan = None
            f.dernier_optimum = None
            f.afficher_plan()
            f.maj_formulaires()
        self.journal.clear()
        self._rafraichir_journal()

    def nouveau(self):
        if (any(f.plan is not None for f in self.familles.values()) or self.journal) and \
                boite_dialogue(self.root, TITRE, "Abandonner le travail en cours et commencer un nouveau projet ?", ("Nouveau projet", "Annuler"), defaut=1) != 0:
            return
        self._vider_tout()
        self.chemin_projet = None
        self.statut("Nouveau projet.")

    def quitter(self):
        if any(f.plan is not None for f in self.familles.values()) and \
                boite_dialogue(self.root, TITRE, "Quitter le programme ? Le travail non enregistré sera perdu.", ("Quitter", "Annuler"), defaut=1) != 0:
            return
        self.root.destroy()


def main():
    try:
        app = Application()
    except Exception as e:
        print(f"Impossible de démarrer l'interface : {e}", file=sys.stderr)
        raise
    app.root.mainloop()


if __name__ == "__main__":
    main()

