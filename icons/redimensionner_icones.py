#!/usr/bin/env python3
"""
Redimensionne des images en icônes carrées PNG à fond transparent.

Étapes pour chaque image :
  1. le fond uniforme (couleur des coins) est rendu transparent ;
  2. l'image est rognée au plus près du motif visible ;
  3. le motif est ajusté dans un carré, centré, sans déformation, avec une marge réglable.

Les originaux ne sont jamais modifiés.

Utilisation :
    python3 redimensionner_icones.py image1.png image2.png
    python3 redimensionner_icones.py dossier_d_images
    python3 redimensionner_icones.py dossier_d_images -t 256 -m 4 -o dossier_de_sortie
    python3 redimensionner_icones.py image.png --garder-fond --pas-de-rognage

Dépendance : Pillow  (sudo apt install python3-pil  ou  pip install pillow)
"""

import argparse
import sys
from pathlib import Path

try:
    from PIL import Image, ImageChops
except ImportError:
    sys.exit("Pillow est requis : sudo apt install python3-pil  (ou pip install pillow)")

EXTENSIONS = {".png", ".jpg", ".jpeg", ".bmp", ".gif", ".webp"}


def collecter(entrees):
    """Retourne la liste des fichiers image à traiter."""
    fichiers = []
    for entree in entrees:
        chemin = Path(entree).expanduser()
        if chemin.is_dir():
            fichiers += sorted(
                f for f in chemin.iterdir()
                if f.is_file() and f.suffix.lower() in EXTENSIONS
            )
        elif chemin.is_file():
            fichiers.append(chemin)
        else:
            print(f"Ignoré (introuvable) : {chemin}")
    return fichiers


def a_du_transparent(img):
    """Indique si l'image contient déjà des pixels transparents."""
    alpha = img.getchannel("A")
    return alpha.getextrema()[0] < 250


def fond_vers_transparent(img):
    """Rend transparent le fond uniforme, d'après la couleur du coin haut-gauche."""
    fond = img.getpixel((0, 0))[:3]
    uni = Image.new("RGB", img.size, fond)
    ecart = ImageChops.difference(img.convert("RGB"), uni).convert("L")
    # Transition douce pour garder des bords lisses
    masque = ecart.point(lambda v: min(255, v * 4))
    img = img.copy()
    img.putalpha(masque)
    return img


def rogner(img):
    """Rogne l'image au plus près des pixels non transparents."""
    masque = img.getchannel("A").point(lambda a: 255 if a > 10 else 0)
    boite = masque.getbbox()
    return img.crop(boite) if boite else img


def creer_icone(source, destination, taille, marge_pct, garder_fond, pas_de_rognage):
    """Crée une icône carrée transparente à partir d'une image."""
    with Image.open(source) as original:
        img = original.convert("RGBA")

    if not garder_fond and not a_du_transparent(img):
        img = fond_vers_transparent(img)
    if not pas_de_rognage:
        img = rogner(img)

    marge = round(taille * marge_pct / 100)
    zone = max(1, taille - 2 * marge)

    # Ajustement proportionnel dans la zone utile
    ratio = min(zone / img.width, zone / img.height)
    nouvelle_taille = (
        max(1, round(img.width * ratio)),
        max(1, round(img.height * ratio)),
    )
    img = img.resize(nouvelle_taille, Image.LANCZOS)

    # Canevas transparent, image centrée
    canevas = Image.new("RGBA", (taille, taille), (0, 0, 0, 0))
    x = (taille - img.width) // 2
    y = (taille - img.height) // 2
    canevas.paste(img, (x, y), img)

    destination.parent.mkdir(parents=True, exist_ok=True)
    canevas.save(destination, "PNG")


def main():
    parser = argparse.ArgumentParser(
        description="Redimensionne des images en icônes carrées PNG transparentes."
    )
    parser.add_argument("entrees", nargs="+",
                        help="fichiers image et/ou dossiers à traiter")
    parser.add_argument("-t", "--taille", type=int, default=128,
                        help="côté de l'icône en pixels (défaut : 128)")
    parser.add_argument("-m", "--marge", type=float, default=3,
                        help="marge en %% de chaque côté (défaut : 3)")
    parser.add_argument("-o", "--sortie", default="icones",
                        help="dossier de sortie (défaut : ./icones)")
    parser.add_argument("--garder-fond", action="store_true",
                        help="ne pas rendre le fond transparent")
    parser.add_argument("--pas-de-rognage", action="store_true",
                        help="ne pas rogner autour du motif")
    args = parser.parse_args()

    if args.taille < 16:
        sys.exit("La taille doit être d'au moins 16 pixels.")
    if not 0 <= args.marge < 40:
        sys.exit("La marge doit être comprise entre 0 et 40 %.")

    fichiers = collecter(args.entrees)
    if not fichiers:
        sys.exit("Aucune image à traiter.")

    dossier_sortie = Path(args.sortie).expanduser()
    ok = 0
    for source in fichiers:
        destination = dossier_sortie / f"{source.stem}_{args.taille}.png"
        try:
            creer_icone(source, destination, args.taille, args.marge,
                        args.garder_fond, args.pas_de_rognage)
            print(f"OK  {source.name} -> {destination}")
            ok += 1
        except Exception as erreur:
            print(f"ERREUR  {source.name} : {erreur}")

    print(f"\n{ok}/{len(fichiers)} icône(s) créée(s) dans {dossier_sortie}")


if __name__ == "__main__":
    main()
