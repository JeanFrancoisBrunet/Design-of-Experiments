# -*- coding: utf-8 -*-
# =============================================================================
#  plans_doe.R — Génération des plans d'expériences
#    criblage (Plackett-Burman, fractions), factoriels complets et fractionnaires,
#    tableaux orthogonaux de Taguchi, plans de surface (composite centré, Box-Behnken),
#    plans de mélange (réseau simplexe, centroïde simplexe, sommets extrêmes).
#  Fonctions de calcul seulement : aucune lecture ni écriture de fichier ici.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# ---------------------------------------------------------------------------
#  Utilitaires communs
# ---------------------------------------------------------------------------
popcount <- function(x, k = 30) sum(as.integer(intToBits(x))[seq_len(k)])

mot_vers_texte <- function(m, k) {
  paste(LETTERS[which(as.integer(intToBits(m))[seq_len(k)] == 1)], collapse = "")
}

valeur_reelle <- function(f, code) {
  if (identical(f$type, "cat")) {
    lv <- unlist(f$niveaux)
    return(ifelse(code < 0, lv[1], lv[2]))
  }
  c0 <- (f$bas + f$haut) / 2
  h <- (f$haut - f$bas) / 2
  signif(c0 + code * h, 10)
}

verifier_facteurs <- function(facteurs, kmin = 2, kmax = 26, quantitatifs = FALSE) {
  k <- length(facteurs)
  if (k < kmin) stop(sprintf("Il faut au moins %d facteurs.", kmin), call. = FALSE)
  if (k > kmax) stop(sprintf("Le nombre de facteurs est limité à %d pour ce plan.", kmax), call. = FALSE)
  noms <- vapply(facteurs, function(f) as.character(f$nom), "")
  if (any(!nzchar(trimws(noms)))) stop("Chaque facteur doit avoir un nom.", call. = FALSE)
  if (anyDuplicated(noms)) stop("Les noms des facteurs doivent être tous différents.", call. = FALSE)
  for (f in facteurs) {
    if (identical(f$type, "cat")) {
      if (quantitatifs) stop(sprintf("Le facteur « %s » doit être quantitatif pour ce plan.", f$nom), call. = FALSE)
      lv <- unlist(f$niveaux)
      if (length(lv) != 2 || lv[1] == lv[2]) stop(sprintf("Le facteur « %s » doit avoir deux niveaux différents.", f$nom), call. = FALSE)
    } else {
      if (!is.numeric(f$bas) || !is.numeric(f$haut) || !is.finite(f$bas) || !is.finite(f$haut))
        stop(sprintf("Niveaux numériques invalides pour le facteur « %s ».", f$nom), call. = FALSE)
      if (f$bas >= f$haut) stop(sprintf("Facteur « %s » : le niveau bas doit être inférieur au niveau haut.", f$nom), call. = FALSE)
    }
  }
  noms
}

# Mise en forme finale : ordre standard, ordre d'exécution (aléatoire ou non), valeurs réelles
finaliser_plan <- function(codes, types, facteurs, aleatoire, graine) {
  codes <- as.matrix(codes)
  n <- nrow(codes)
  if (isTRUE(aleatoire)) {
    if (!is.null(graine) && is.finite(graine)) set.seed(as.integer(graine))
    ord <- sample.int(n)
  } else ord <- seq_len(n)
  df <- data.frame(Ordre_standard = seq_len(n), Type_point = types, stringsAsFactors = FALSE)
  for (j in seq_along(facteurs)) df[[facteurs[[j]]$nom]] <- valeur_reelle(facteurs[[j]], codes[, j])
  df <- df[ord, , drop = FALSE]
  df <- cbind(Ordre_essai = seq_len(n), df)
  rownames(df) <- NULL
  df
}

ajouter_repetitions_centres <- function(M, repetitions, points_centraux, facteurs, types_base = "Cube") {
  repetitions <- max(1L, as.integer(repetitions %||% 1))
  points_centraux <- max(0L, as.integer(points_centraux %||% 0))
  codes <- M[rep(seq_len(nrow(M)), repetitions), , drop = FALSE]
  types <- rep(types_base, nrow(codes))
  if (points_centraux > 0) {
    if (any(vapply(facteurs, function(f) identical(f$type, "cat"), TRUE)))
      stop("Les points centraux exigent des facteurs tous quantitatifs.", call. = FALSE)
    codes <- rbind(codes, matrix(0, points_centraux, ncol(M)))
    types <- c(types, rep("Centre", points_centraux))
  }
  list(codes = codes, types = types)
}

# ---------------------------------------------------------------------------
#  Plans à deux niveaux
# ---------------------------------------------------------------------------
base_2niveaux <- function(k) {
  m <- as.matrix(expand.grid(rep(list(c(-1, 1)), k)))
  dimnames(m) <- list(NULL, LETTERS[seq_len(k)])
  m
}

hadamard_sylvester <- function(n) {
  H <- matrix(1, 1, 1)
  while (nrow(H) < n) H <- rbind(cbind(H, H), cbind(H, -H))
  H
}

PB_LIGNES <- list(
  "8"  = c(1, 1, 1, -1, 1, -1, -1),
  "12" = c(1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1),
  "20" = c(1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, 1, -1),
  "24" = c(1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1))
TAILLES_PB <- c(8, 12, 16, 20, 24, 32)

matrice_plackett_burman <- function(n) {
  if (as.character(n) %in% names(PB_LIGNES)) {
    g <- PB_LIGNES[[as.character(n)]]
    M <- t(vapply(0:(n - 2), function(s) g[((seq_along(g) - 1 - s) %% length(g)) + 1], numeric(n - 1)))
    M <- rbind(M, rep(-1, n - 1))
  } else if (n %in% c(16, 32)) {
    M <- hadamard_sylvester(n)[, -1, drop = FALSE]
    M <- M[order(do.call(order, as.data.frame(M[, rev(seq_len(ncol(M)))]))), , drop = FALSE]
  } else stop("Taille de plan de Plackett-Burman non disponible.", call. = FALSE)
  colnames(M) <- NULL
  M
}

# --- fractions : générateurs « E=ABC », « F=-BCD » -------------------------
analyser_generateurs <- function(k, generateurs) {
  generateurs <- trimws(unlist(generateurs))
  generateurs <- generateurs[nzchar(generateurs)]
  p <- length(generateurs)
  kb <- k - p
  if (p < 1 || kb < 2) stop("Nombre de générateurs incohérent avec le nombre de facteurs.", call. = FALSE)
  bit <- function(ch) bitwShiftL(1L, utf8ToInt(ch) - 65L)
  mots <- integer(p); signes <- numeric(p); cibles <- character(p); sources <- vector("list", p)
  for (j in seq_len(p)) {
    g <- toupper(gsub("[[:space:]]", "", generateurs[j]))
    parts <- strsplit(g, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2 || nchar(parts[1]) != 1)
      stop(sprintf("Générateur invalide : « %s » (forme attendue : E=ABC).", generateurs[j]), call. = FALSE)
    cible <- parts[1]; rhs <- parts[2]; s <- 1
    if (substr(rhs, 1, 1) == "-") { s <- -1; rhs <- substring(rhs, 2) } else if (substr(rhs, 1, 1) == "+") rhs <- substring(rhs, 2)
    if (cible != LETTERS[kb + j])
      stop(sprintf("Le générateur n° %d doit définir le facteur %s (et non %s).", j, LETTERS[kb + j], cible), call. = FALSE)
    lettres <- strsplit(rhs, "")[[1]]
    if (length(lettres) < 2 || any(!(lettres %in% LETTERS[seq_len(kb + j - 1)])) || anyDuplicated(lettres))
      stop(sprintf("Générateur invalide : « %s » (il faut au moins deux facteurs distincts parmi A à %s).", generateurs[j], LETTERS[kb + j - 1]), call. = FALSE)
    m <- bit(cible)
    for (ch in lettres) m <- bitwXor(m, bit(ch))
    mots[j] <- m; signes[j] <- s; cibles[j] <- cible; sources[[j]] <- lettres
  }
  # sous-groupe de définition : toutes les combinaisons non vides de mots
  mots_tous <- integer(0); signes_tous <- numeric(0)
  for (s in seq_len(2^p - 1)) {
    sel <- which(as.integer(intToBits(s))[seq_len(p)] == 1)
    m <- 0L; sg <- 1
    for (j in sel) { m <- bitwXor(m, mots[j]); sg <- sg * signes[j] }
    mots_tous <- c(mots_tous, m); signes_tous <- c(signes_tous, sg)
  }
  longueurs <- vapply(mots_tous, popcount, 1, k = k)
  list(k = k, p = p, kb = kb, mots = mots, signes = signes, cibles = cibles, sources = sources,
       mots_tous = mots_tous, signes_tous = signes_tous, longueurs = longueurs, resolution = min(longueurs))
}

romain <- function(n) c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")[n]

fraction_factorielle <- function(k, generateurs) {
  info <- analyser_generateurs(k, generateurs)
  if (info$resolution < 3) stop("Ces générateurs donnent un plan de résolution inférieure à III : des effets principaux seraient confondus entre eux.", call. = FALSE)
  M <- base_2niveaux(info$kb)
  X <- M
  for (j in seq_len(info$p)) {
    col <- info$signes[j] * apply(X[, match(info$sources[[j]], LETTERS[seq_len(ncol(X))]), drop = FALSE], 1, prod)
    X <- cbind(X, col)
    colnames(X)[ncol(X)] <- info$cibles[j]
  }
  list(M = X, info = info)
}

# Structure des alias (effets de 1 à 2 facteurs, alias jusqu'à l'ordre 3)
structure_alias <- function(info, ordre_max_alias = 3, effets_max = 2) {
  k <- info$k
  effets <- integer(0)
  for (o in seq_len(effets_max)) {
    cmb <- utils::combn(k, o)
    effets <- c(effets, apply(cmb, 2, function(ix) sum(2^(ix - 1))))
  }
  effets <- as.integer(effets)
  vus <- integer(0); lignes <- character(0)
  for (e in effets) {
    if (e %in% vus) next
    membres <- e; signes <- 1
    for (i in seq_along(info$mots_tous)) {
      m <- bitwXor(e, info$mots_tous[i])
      if (popcount(m, k) <= ordre_max_alias && !(m %in% membres)) { membres <- c(membres, m); signes <- c(signes, info$signes_tous[i]) }
    }
    vus <- c(vus, membres)
    ord <- order(vapply(membres, popcount, 1, k = k), membres)
    membres <- membres[ord]; signes <- signes[ord]
    txt <- mot_vers_texte(membres[1], k)
    if (length(membres) > 1) for (i in 2:length(membres))
      txt <- paste0(txt, if (signes[i] > 0) " + " else " − ", mot_vers_texte(membres[i], k))
    lignes <- c(lignes, txt)
  }
  lignes
}

relation_definition <- function(info) {
  ord <- order(info$longueurs, info$mots_tous)
  paste0("I = ", paste0(ifelse(info$signes_tous[ord] > 0, "", "−"), vapply(info$mots_tous[ord], mot_vers_texte, "", k = info$k), collapse = " = "))
}

# --- criblage -------------------------------------------------------------
plan_criblage <- function(par) {
  facteurs <- par$facteurs
  noms <- verifier_facteurs(facteurs, 2, 31)
  k <- length(facteurs)
  methode <- par$methode %||% "pb"
  lignes <- character(0); tables <- list(); avert <- character(0); meta <- list()
  if (methode == "pb") {
    n <- as.integer(par$nb_essais %||% 0)
    if (n == 0) n <- TAILLES_PB[which(TAILLES_PB - 1 >= k)[1]]
    if (is.na(n) || !(n %in% TAILLES_PB)) stop("Taille de plan de Plackett-Burman non disponible (8, 12, 16, 20, 24 ou 32 essais).", call. = FALSE)
    if (k > n - 1) stop(sprintf("Un plan de %d essais accepte au plus %d facteurs.", n, n - 1), call. = FALSE)
    Mfull <- matrice_plackett_burman(n)
    M <- Mfull[, seq_len(k), drop = FALSE]
    titre <- sprintf("Plan de criblage de Plackett-Burman — %d essais", n)
    lignes <- c("Plan de résolution III : les effets principaux sont confondus avec des interactions entre deux facteurs ou plus. Il sert à repérer les facteurs influents, non à estimer les interactions.",
                sprintf("Degrés de liberté disponibles pour l'erreur : %d (avant répétitions et points centraux).", n - 1 - k))
    meta$methode <- "pb"; meta$resolution <- "III"; meta$nb_essais_base <- n
  } else {
    gens <- unlist(par$generateurs)
    if (length(gens) == 0) stop("Aucun générateur : choisir une fraction.", call. = FALSE)
    fr <- fraction_factorielle(k, gens)
    M <- fr$M; info <- fr$info
    titre <- sprintf("Plan de criblage fractionnaire 2^(%d−%d) — %d essais", k, info$p, nrow(M))
    lignes <- c(sprintf("Résolution %s.", romain(info$resolution)), paste0("Générateurs : ", paste(gens, collapse = " ; ")),
                relation_definition(info))
    tables[["Structure_des_alias"]] <- list(Alias = structure_alias(info))
    meta$methode <- "fraction"; meta$resolution <- romain(info$resolution); meta$generateurs <- as.list(gens)
  }
  ar <- ajouter_repetitions_centres(M, par$repetitions, par$points_centraux, facteurs)
  plan <- finaliser_plan(ar$codes, ar$types, facteurs, par$aleatoire, par$graine)
  lignes <- c(sprintf("%s", titre), sprintf("Facteurs : %d — essais : %d (répétitions : %d, points centraux : %d).", k, nrow(plan),
                                           as.integer(par$repetitions %||% 1), as.integer(par$points_centraux %||% 0)),
              paste0("Codes : ", paste(sprintf("%s = %s", LETTERS[seq_len(k)], noms), collapse = " ; ")), lignes)
  if (isTRUE(par$aleatoire)) lignes <- c(lignes, "Ordre d'exécution : aléatoire (l'ordre standard est conservé dans la colonne « Ordre_standard »).")
  meta$famille <- "criblage"
  list(titre = titre, lignes = lignes, tables = tables, plan = plan, meta = meta)
}

# --- factoriel ------------------------------------------------------------
plan_factoriel <- function(par) {
  facteurs <- par$facteurs
  noms <- verifier_facteurs(facteurs, 2, 11)
  k <- length(facteurs)
  type <- par$type %||% "complet"
  tables <- list(); meta <- list(famille = "factoriel")
  if (type == "complet") {
    if (k > 10) stop("Un plan factoriel complet est limité à 10 facteurs.", call. = FALSE)
    M <- base_2niveaux(k); titre <- sprintf("Plan factoriel complet 2^%d — %d essais de base", k, nrow(M))
    lignes <- "Plan de résolution complète : tous les effets principaux et toutes les interactions sont estimables."
    meta$resolution <- "complet"
  } else {
    gens <- unlist(par$generateurs)
    if (length(gens) == 0) stop("Aucun générateur : choisir une fraction.", call. = FALSE)
    fr <- fraction_factorielle(k, gens)
    M <- fr$M; info <- fr$info
    titre <- sprintf("Plan factoriel fractionnaire 2^(%d−%d), résolution %s — %d essais de base", k, info$p, romain(info$resolution), nrow(M))
    interp <- switch(as.character(info$resolution),
      "3" = "Résolution III : les effets principaux sont confondus avec des interactions entre deux facteurs.",
      "4" = "Résolution IV : les effets principaux ne sont confondus qu'avec des interactions d'ordre 3 ou plus ; les interactions entre deux facteurs sont confondues entre elles.",
      "Résolution V ou plus : les effets principaux et les interactions entre deux facteurs sont estimables, sans confusion entre eux.")
    lignes <- c(interp, paste0("Générateurs : ", paste(gens, collapse = " ; ")), relation_definition(info))
    tables[["Structure_des_alias"]] <- list(Alias = structure_alias(info))
    meta$resolution <- romain(info$resolution); meta$generateurs <- as.list(gens)
  }
  ar <- ajouter_repetitions_centres(M, par$repetitions, par$points_centraux, facteurs)
  plan <- finaliser_plan(ar$codes, ar$types, facteurs, par$aleatoire, par$graine)
  lignes <- c(titre, sprintf("Facteurs : %d — essais : %d (répétitions : %d, points centraux : %d).", k, nrow(plan),
                             as.integer(par$repetitions %||% 1), as.integer(par$points_centraux %||% 0)),
              paste0("Codes : ", paste(sprintf("%s = %s", LETTERS[seq_len(k)], noms), collapse = " ; ")), lignes)
  if (isTRUE(par$aleatoire)) lignes <- c(lignes, "Ordre d'exécution : aléatoire (l'ordre standard est conservé dans la colonne « Ordre_standard »).")
  list(titre = titre, lignes = lignes, tables = tables, plan = plan, meta = meta)
}

# ---------------------------------------------------------------------------
#  Plans de surface de réponse
# ---------------------------------------------------------------------------
plan_surface <- function(par) {
  facteurs <- par$facteurs
  noms <- verifier_facteurs(facteurs, 2, 7, quantitatifs = TRUE)
  k <- length(facteurs)
  type <- par$type %||% "ccc"
  nc_defaut <- c(5, 6, 7, 10, 15, 20)[k - 1]
  nc <- as.integer(par$points_centraux %||% nc_defaut)
  rep <- max(1L, as.integer(par$repetitions %||% 1))
  meta <- list(famille = "surface", type_plan = type)
  if (type %in% c("ccc", "ccf", "cci")) {
    cube <- if (k <= 4) base_2niveaux(k) else fraction_factorielle(k, c("E=ABCD", "F=ABCDE", "G=ABCDEF")[k - 4])$M
    nf <- nrow(cube)
    a_rot <- nf^0.25
    alpha <- switch(type, ccf = 1, cci = a_rot, {
      mode <- par$alpha_mode %||% "rotatable"
      if (mode == "personnalise") as.numeric(par$alpha) else if (mode == "orthogonal") {
        tot <- nf + 2 * k + nc
        sqrt((sqrt(nf * tot) - nf) / 2)
      } else a_rot })
    if (!is.finite(alpha) || alpha <= 0) stop("Valeur de α invalide.", call. = FALSE)
    cubec <- if (type == "cci") cube / alpha else cube
    ax <- if (type == "cci") 1 else alpha
    axiaux <- do.call(rbind, lapply(seq_len(k), function(j) { a <- matrix(0, 2, k); a[1, j] <- -ax; a[2, j] <- ax; a }))
    codes <- rbind(cubec[rep(seq_len(nf), rep), , drop = FALSE], axiaux[rep(seq_len(nrow(axiaux)), rep), , drop = FALSE], matrix(0, nc, k))
    types <- c(rep("Cube", nf * rep), rep("Axial", 2 * k * rep), rep("Centre", nc))
    nom_type <- c(ccc = "Plan composite centré circonscrit (CCC)", ccf = "Plan composite centré à faces centrées (CCF)", cci = "Plan composite centré inscrit (CCI)")[[type]]
    titre <- sprintf("%s — %d facteurs, %d essais", nom_type, k, nrow(codes))
    lignes <- c(sprintf("Partie factorielle : %d essais%s ; points axiaux : %d ; points centraux : %d.", nf, if (k >= 5) " (fraction de résolution V ou plus)" else "", 2 * k, nc),
                sprintf("Distance axiale α = %s (unités codées).", format(round(ax, 4), decimal.mark = ",")),
                if (type == "ccc") "Niveaux bas et haut saisis = niveaux de la partie factorielle ; les points axiaux sortent de ce domaine." else
                if (type == "ccf") "Niveaux bas et haut saisis = limites du domaine ; les points axiaux sont au centre des faces." else
                "Niveaux bas et haut saisis = limites du domaine (points axiaux) ; la partie factorielle est réduite d'un facteur α.",
                if (type == "ccc" && abs(alpha - a_rot) < 1e-9) "Plan isovariant par rotation (α = racine quatrième du nombre d'essais factoriels)." else NULL)
    meta$alpha <- alpha; meta$essais_cube <- nf
  } else if (type == "bb") {
    if (k < 3 || k > 5) stop("Le plan de Box-Behnken est proposé pour 3 à 5 facteurs.", call. = FALSE)
    paires <- utils::combn(k, 2)
    blocs <- lapply(seq_len(ncol(paires)), function(i) {
      m <- matrix(0, 4, k); g <- as.matrix(expand.grid(c(-1, 1), c(-1, 1)))
      m[, paires[1, i]] <- g[, 1]; m[, paires[2, i]] <- g[, 2]; m })
    base <- do.call(rbind, blocs)
    nc <- as.integer(par$points_centraux %||% c(3, 3, 6)[k - 2])
    codes <- rbind(base[rep(seq_len(nrow(base)), rep), , drop = FALSE], matrix(0, nc, k))
    types <- c(rep("Arête", nrow(base) * rep), rep("Centre", nc))
    titre <- sprintf("Plan de Box-Behnken — %d facteurs, %d essais", k, nrow(codes))
    lignes <- c(sprintf("Trois niveaux par facteur ; %d essais sur les milieux des arêtes du cube et %d points centraux.", nrow(base), nc),
                "Aucun essai aux sommets du cube : le plan convient quand les combinaisons extrêmes sont irréalisables ou dangereuses.")
  } else stop("Type de plan de surface inconnu.", call. = FALSE)
  plan <- finaliser_plan(codes, types, facteurs, par$aleatoire, par$graine)
  lignes <- c(titre, paste0("Codes : ", paste(sprintf("%s = %s", LETTERS[seq_len(k)], noms), collapse = " ; ")), lignes)
  if (isTRUE(par$aleatoire)) lignes <- c(lignes, "Ordre d'exécution : aléatoire (l'ordre standard est conservé dans la colonne « Ordre_standard »).")
  list(titre = titre, lignes = lignes, tables = list(), plan = plan, meta = meta)
}

# ---------------------------------------------------------------------------
#  Plans de mélange
# ---------------------------------------------------------------------------
sommets_extremes <- function(L, U) {
  q <- length(L); tol <- 1e-9
  fixes <- which(U - L < tol)
  res <- list()
  for (libre in seq_len(q)) {
    autres <- setdiff(seq_len(q), libre)
    g <- as.matrix(expand.grid(rep(list(c(0, 1)), q - 1)))
    for (i in seq_len(nrow(g))) {
      x <- numeric(q)
      x[autres] <- ifelse(g[i, ] == 0, L[autres], U[autres])
      x[libre] <- 1 - sum(x[autres])
      if (x[libre] >= L[libre] - tol && x[libre] <= U[libre] + tol) res[[length(res) + 1]] <- round(x, 10)
    }
  }
  if (!length(res)) return(NULL)
  m <- do.call(rbind, res)
  m <- m[!duplicated(round(m, 8)), , drop = FALSE]
  m
}

actifs_sommet <- function(x, L, U, tol = 1e-8) c(which(abs(x - L) < tol), length(L) + which(abs(x - U) < tol))

verifier_melange <- function(L, U) {
  if (any(L < -1e-12) || any(U > 1 + 1e-12) || any(L > U + 1e-12)) stop("Bornes des composants invalides (0 ≤ bas ≤ haut ≤ 1).", call. = FALSE)
  if (sum(L) > 1 + 1e-9) stop("La somme des bornes inférieures dépasse 1 : aucun mélange n'est possible.", call. = FALSE)
  if (sum(U) < 1 - 1e-9) stop("La somme des bornes supérieures est inférieure à 1 : aucun mélange n'est possible.", call. = FALSE)
}

compositions <- function(q, m) {
  res <- list()
  rec <- function(prefix, reste, pos) {
    if (pos == q) { res[[length(res) + 1]] <<- c(prefix, reste); return(invisible()) }
    for (v in 0:reste) rec(c(prefix, v), reste - v, pos + 1)
  }
  rec(numeric(0), m, 1)
  do.call(rbind, res)
}

distance_max_rayon <- function(c0, d, L, U) {
  t <- Inf
  for (j in seq_along(d)) {
    if (d[j] > 1e-12) t <- min(t, (U[j] - c0[j]) / d[j])
    if (d[j] < -1e-12) t <- min(t, (L[j] - c0[j]) / d[j])
  }
  t
}

plan_melange <- function(par) {
  comp <- par$composants
  q <- length(comp)
  if (q < 2) stop("Un mélange comporte au moins deux composants.", call. = FALSE)
  if (q > 10) stop("Le nombre de composants est limité à 10.", call. = FALSE)
  noms <- vapply(comp, function(f) as.character(f$nom), "")
  if (any(!nzchar(trimws(noms))) || anyDuplicated(noms)) stop("Les composants doivent avoir des noms non vides et différents.", call. = FALSE)
  L <- vapply(comp, function(f) as.numeric(f$bas %||% 0), 1)
  U <- vapply(comp, function(f) as.numeric(f$haut %||% 1), 1)
  verifier_melange(L, U)
  type <- par$type %||% "reseau"
  rep <- max(1L, as.integer(par$repetitions %||% 1))
  nc <- max(0L, as.integer(par$points_centraux %||% 1))
  sL <- sum(L); etendue <- 1 - sL
  Ueff <- 1 - (sL - L)                      # borne supérieure implicite quand seules les bornes inférieures sont imposées
  if (type %in% c("reseau", "centroide")) {
    if (any(U < Ueff - 1e-9)) stop("Des bornes supérieures contraignent le domaine : utiliser le plan à sommets extrêmes.", call. = FALSE)
    Ubord <- Ueff
  } else Ubord <- pmin(U, Ueff)
  pts <- NULL; types <- character(0); lignes <- character(0)
  if (type == "reseau") {
    m <- as.integer(par$degre %||% 2)
    if (m < 1 || m > 4) stop("Le degré du réseau doit être compris entre 1 et 4.", call. = FALSE)
    cp <- compositions(q, m) / m
    if (nrow(cp) > 600) stop("Ce réseau compte plus de 600 points : réduire le nombre de composants ou le degré.", call. = FALSE)
    pts <- sweep(cp * etendue, 2, L, "+")
    types <- ifelse(rowSums(cp > 1e-9) == 1, "Sommet", "Réseau")
    titre <- sprintf("Plan de mélange en réseau simplexe {%d, %d}", q, m)
    lignes <- c(sprintf("Points : proportions multiples de 1/%d (avec bornes inférieures éventuelles) : %d points.", m, nrow(pts)),
                if (m == 1) "Degré 1 : seuls les sommets, adapté au modèle linéaire." else if (m == 2) "Degré 2 : adapté au modèle quadratique." else "Degré 3 ou plus : permet le modèle cubique spécial ou complet.")
  } else if (type == "centroide") {
    if (q > 8) stop("Le centroïde simplexe est limité à 8 composants.", call. = FALSE)
    pts <- do.call(rbind, lapply(seq_len(2^q - 1), function(s) {
      b <- as.integer(intToBits(s))[seq_len(q)]; x <- b / sum(b); L + x * etendue }))
    types <- ifelse(vapply(seq_len(2^q - 1), function(s) sum(as.integer(intToBits(s))[seq_len(q)]), 1) == 1, "Sommet",
                    ifelse(seq_len(2^q - 1) == 2^q - 1, "Centre", "Centroïde"))
    titre <- sprintf("Plan de mélange en centroïde simplexe — %d composants", q)
    lignes <- sprintf("Points : sommets, milieux des arêtes, centres des faces et centre global (%d points).", nrow(pts))
  } else if (type == "sommets") {
    V <- sommets_extremes(L, U)
    if (is.null(V)) stop("Aucun sommet : les contraintes sont incompatibles.", call. = FALSE)
    if (nrow(V) > 300) stop("Plus de 300 sommets : le plan est trop grand.", call. = FALSE)
    pts <- V; types <- rep("Sommet", nrow(V))
    act <- lapply(seq_len(nrow(V)), function(i) actifs_sommet(V[i, ], L, U))
    if (isTRUE(par$aretes) && nrow(V) > 1) {
      for (i in 1:(nrow(V) - 1)) for (j in (i + 1):nrow(V))
        if (length(intersect(act[[i]], act[[j]])) >= q - 2) { pts <- rbind(pts, (V[i, ] + V[j, ]) / 2); types <- c(types, "Arête") }
    }
    if (isTRUE(par$faces)) {
      for (cst in seq_len(2 * q)) {
        ix <- which(vapply(act, function(a) cst %in% a, TRUE))
        if (length(ix) >= 3 && length(ix) < nrow(V)) { pts <- rbind(pts, colMeans(V[ix, , drop = FALSE])); types <- c(types, "Face") }
      }
    }
    titre <- sprintf("Plan de mélange à sommets extrêmes — %d composants", q)
    lignes <- c(sprintf("Sommets du domaine contraint : %d.", nrow(V)),
                if (isTRUE(par$aretes)) "Milieux des arêtes ajoutés." else NULL, if (isTRUE(par$faces)) "Centres des faces ajoutés." else NULL)
  } else stop("Type de plan de mélange inconnu.", call. = FALSE)
  # centre global et points axiaux, dans le domaine défini par les bornes effectives
  Vc <- if (type == "sommets") sommets_extremes(L, U) else sommets_extremes(L, Ubord)
  centre <- if (is.null(Vc)) (L + (1 - sL) / q) else colMeans(Vc)
  if (isTRUE(par$axiaux)) {
    for (i in seq_len(q)) {
      d <- -centre; d[i] <- d[i] + 1
      t <- distance_max_rayon(centre, d, L, if (type == "sommets") U else Ubord)
      if (is.finite(t) && t > 1e-9) { pts <- rbind(pts, centre + 0.5 * t * d); types <- c(types, "Axial") }
    }
    lignes <- c(lignes, "Points axiaux ajoutés (à mi-distance entre le centre et le sommet de chaque composant).")
  }
  pts <- pts[rep(seq_len(nrow(pts)), rep), , drop = FALSE]; types <- rep(types, rep)
  if (nc > 0) { pts <- rbind(pts, matrix(centre, nc, q, byrow = TRUE)); types <- c(types, rep("Centre", nc)) }
  pts <- round(pts, 10)
  n <- nrow(pts)
  if (isTRUE(par$aleatoire)) { if (!is.null(par$graine) && is.finite(par$graine)) set.seed(as.integer(par$graine)); ord <- sample.int(n) } else ord <- seq_len(n)
  df <- data.frame(Ordre_standard = seq_len(n), Type_point = types, stringsAsFactors = FALSE)
  for (j in seq_len(q)) df[[noms[j]]] <- pts[, j]
  df <- df[ord, , drop = FALSE]
  df <- cbind(Ordre_essai = seq_len(n), df); rownames(df) <- NULL
  lg <- c(titre, sprintf("Composants : %d — essais : %d (répétitions : %d, points centraux : %d).", q, n, rep, nc))
  if (any(L > 0) || any(U < 1)) lg <- c(lg, "Bornes : " %+% paste(sprintf("%s [%s ; %s]", noms, format(L, decimal.mark = ",", trim = TRUE), format(U, decimal.mark = ",", trim = TRUE)), collapse = " ; "))
  lg <- c(lg, lignes)
  if (isTRUE(par$aleatoire)) lg <- c(lg, "Ordre d'exécution : aléatoire (l'ordre standard est conservé dans la colonne « Ordre_standard »).")
  meta <- list(famille = "melange", type_plan = type, bornes_bas = as.list(L), bornes_haut = as.list(U))
  list(titre = titre, lignes = lg, tables = list(), plan = df, meta = meta)
}
`%+%` <- function(a, b) paste0(a, b)

# ---------------------------------------------------------------------------
#  Tableaux orthogonaux de Taguchi
# ---------------------------------------------------------------------------
# Tableaux de Galois (p premier) : colonnes dans l'ordre classique de Taguchi
vecteurs_gf <- function(p, m) {
  res <- list()
  for (j in seq_len(m)) {
    if (j == 1) { res[[length(res) + 1]] <- c(1, rep(0, m - 1)); next }
    g <- as.matrix(expand.grid(rep(list(0:(p - 1)), j - 1)))
    for (i in seq_len(nrow(g))) res[[length(res) + 1]] <- c(g[i, ], 1, rep(0, m - j))
  }
  res
}

normaliser_gf <- function(v, p) {
  nz <- which(v %% p != 0)
  if (!length(nz)) return(NULL)
  dern <- v[nz[length(nz)]] %% p
  inv <- which((dern * seq_len(p - 1)) %% p == 1)
  (v * inv) %% p
}

tableau_gf <- function(p, m) {
  n <- p^m
  base <- as.matrix(expand.grid(rep(list(0:(p - 1)), m)))[, m:1, drop = FALSE]   # première colonne = la plus lente
  vec <- vecteurs_gf(p, m)
  M <- vapply(vec, function(v) as.integer((base %*% v) %% p) + 1L, integer(n))
  list(M = M, vecteurs = vec, p = p)
}

L12 <- function() { M <- matrice_plackett_burman(12); (M + 3) / 2 }
L18 <- function() {
  txt <- "1 1 1 1 1 1 1 1
          1 1 2 2 2 2 2 2
          1 1 3 3 3 3 3 3
          1 2 1 1 2 2 3 3
          1 2 2 2 3 3 1 1
          1 2 3 3 1 1 2 2
          1 3 1 2 1 3 2 3
          1 3 2 3 2 1 3 1
          1 3 3 1 3 2 1 2
          2 1 1 3 3 2 2 1
          2 1 2 1 1 3 3 2
          2 1 3 2 2 1 1 3
          2 2 1 2 3 1 3 2
          2 2 2 3 1 2 1 3
          2 2 3 1 2 3 2 1
          2 3 1 3 2 3 1 2
          2 3 2 1 3 1 2 3
          2 3 3 2 1 2 3 1"
  matrix(as.integer(strsplit(gsub("\n", " ", txt), "[[:space:]]+")[[1]]), ncol = 8, byrow = TRUE)
}

ARRAYS <- list(
  L4  = list(runs = 4,  niv = rep(2, 3),  gf = c(2, 2)),
  L8  = list(runs = 8,  niv = rep(2, 7),  gf = c(2, 3)),
  L9  = list(runs = 9,  niv = rep(3, 4),  gf = c(3, 2)),
  L12 = list(runs = 12, niv = rep(2, 11), gf = NULL),
  L16 = list(runs = 16, niv = rep(2, 15), gf = c(2, 4)),
  L18 = list(runs = 18, niv = c(2, rep(3, 7)), gf = NULL),
  L25 = list(runs = 25, niv = rep(5, 6),  gf = c(5, 2)),
  L27 = list(runs = 27, niv = rep(3, 13), gf = c(3, 3)),
  L32 = list(runs = 32, niv = rep(2, 31), gf = c(2, 5)))

obtenir_tableau <- function(nom) {
  a <- ARRAYS[[nom]]
  if (is.null(a)) stop("Tableau orthogonal inconnu.", call. = FALSE)
  if (!is.null(a$gf)) { t <- tableau_gf(a$gf[1], a$gf[2]); list(M = t$M, niv = a$niv, vecteurs = t$vecteurs, p = a$gf[1]) }
  else if (nom == "L12") list(M = matrix(as.integer(L12()), 12), niv = a$niv, vecteurs = NULL, p = 2)
  else list(M = L18(), niv = a$niv, vecteurs = NULL, p = NA)
}

colonnes_interaction <- function(arr, ci, cj) {
  p <- arr$p
  res <- integer(0)
  for (s in seq_len(p - 1)) {
    v <- normaliser_gf(arr$vecteurs[[ci]] + s * arr$vecteurs[[cj]], p)
    idx <- which(vapply(arr$vecteurs, function(w) all(w == v), TRUE))
    res <- c(res, idx)
  }
  sort(unique(res))
}

# Recherche d'une affectation des facteurs aux colonnes (retour arrière)
affecter_colonnes <- function(arr, niv_f, inter, tolere_fictif) {
  nf <- length(niv_f); nc <- length(arr$niv)
  noeuds <- 0L
  choix <- integer(nf)
  utilisees_ic <- integer(0)
  fictif <- logical(nf)
  rec <- function(i) {
    noeuds <<- noeuds + 1L
    if (noeuds > 300000L) return(FALSE)
    if (i > nf) return(TRUE)
    cand <- which(arr$niv == niv_f[i] & !(seq_len(nc) %in% c(choix[seq_len(i - 1)], utilisees_ic)))
    cand_f <- integer(0)
    if (tolere_fictif && niv_f[i] == 2) cand_f <- which(arr$niv == 3 & !(seq_len(nc) %in% c(choix[seq_len(i - 1)], utilisees_ic)))
    for (c in c(cand, cand_f)) {
      choix[i] <<- c; fictif[i] <<- !(arr$niv[c] == niv_f[i])
      # interactions dont les deux facteurs sont maintenant affectés
      ok <- TRUE; ajout <- integer(0)
      for (it in inter) {
        if (max(it) == i) {
          ci <- choix[it[1]]; cj <- choix[it[2]]
          cols <- colonnes_interaction(arr, ci, cj)
          if (length(cols) == 0 || any(cols %in% c(choix[seq_len(i)], utilisees_ic, ajout))) { ok <- FALSE; break }
          ajout <- c(ajout, cols)
        }
      }
      if (ok) {
        sauv <- utilisees_ic
        utilisees_ic <<- c(utilisees_ic, ajout)
        if (rec(i + 1)) return(TRUE)
        utilisees_ic <<- sauv
      }
    }
    choix[i] <<- 0L
    FALSE
  }
  # interactions : seulement dans les tableaux de Galois, et entre facteurs de même nombre de niveaux
  if (length(inter) && (is.null(arr$vecteurs))) return(NULL)
  if (rec(1)) list(colonnes = choix, fictif = fictif) else NULL
}

plan_taguchi <- function(par) {
  fac <- par$facteurs
  k <- length(fac)
  if (k < 2) stop("Il faut au moins deux facteurs.", call. = FALSE)
  noms <- vapply(fac, function(f) as.character(f$nom), "")
  if (any(!nzchar(trimws(noms))) || anyDuplicated(noms)) stop("Les facteurs doivent avoir des noms non vides et différents.", call. = FALSE)
  niveaux <- lapply(fac, function(f) as.character(unlist(f$niveaux)))
  nv <- vapply(niveaux, length, 1L)
  if (any(nv < 2 | nv > 5 | nv == 4)) stop("Chaque facteur doit avoir 2, 3 ou 5 niveaux (4 niveaux : non pris en charge).", call. = FALSE)
  for (i in seq_len(k)) if (anyDuplicated(niveaux[[i]])) stop(sprintf("Les niveaux du facteur « %s » doivent être différents.", noms[i]), call. = FALSE)
  inter <- lapply(par$interactions %||% list(), function(it) as.integer(unlist(it)))
  inter <- lapply(inter, function(it) sort(it))
  for (it in inter) {
    if (length(it) != 2 || any(it < 1 | it > k) || it[1] == it[2]) stop("Interaction invalide.", call. = FALSE)
    if (nv[it[1]] != nv[it[2]]) stop("Une interaction ne peut être étudiée qu'entre facteurs ayant le même nombre de niveaux.", call. = FALSE)
  }
  choix_t <- par$tableau %||% "auto"
  candidats <- if (choix_t == "auto") names(ARRAYS)[order(vapply(ARRAYS, function(a) a$runs, 1))] else choix_t
  trouve <- NULL
  for (nm in candidats) {
    arr <- obtenir_tableau(nm); arr$nom <- nm
    if (length(arr$niv) < k) next
    for (fict in c(FALSE, TRUE)) {
      aff <- affecter_colonnes(arr, nv, inter, fict)
      if (!is.null(aff)) { trouve <- list(arr = arr, aff = aff); break }
    }
    if (!is.null(trouve)) break
  }
  if (is.null(trouve)) stop(if (choix_t == "auto") "Aucun tableau orthogonal disponible ne convient à ces facteurs (et interactions)." else
    sprintf("Le tableau %s ne convient pas à ces facteurs et interactions.", choix_t), call. = FALSE)
  arr <- trouve$arr; aff <- trouve$aff
  Mcol <- arr$M[, aff$colonnes, drop = FALSE]
  n <- nrow(Mcol)
  codes_niv <- Mcol
  # niveau fictif : le niveau 3 reprend le niveau 1
  for (i in seq_len(k)) if (aff$fictif[i]) codes_niv[codes_niv[, i] == 3L, i] <- 1L
  df <- data.frame(Ordre_standard = seq_len(n), stringsAsFactors = FALSE)
  for (i in seq_len(k)) df[[noms[i]]] <- niveaux[[i]][codes_niv[, i]]
  if (isTRUE(par$aleatoire)) { if (!is.null(par$graine) && is.finite(par$graine)) set.seed(as.integer(par$graine)); ord <- sample.int(n) } else ord <- seq_len(n)
  df <- df[ord, , drop = FALSE]
  df <- cbind(Ordre_essai = seq_len(n), df); rownames(df) <- NULL
  # description des colonnes
  desc_col <- rep("", length(arr$niv)); ncol_total <- length(arr$niv)
  for (i in seq_len(k)) desc_col[aff$colonnes[i]] <- noms[i]
  intinfo <- list()
  for (it in inter) {
    cols <- colonnes_interaction(arr, aff$colonnes[it[1]], aff$colonnes[it[2]])
    for (cc in cols) desc_col[cc] <- paste0(noms[it[1]], " × ", noms[it[2]])
    intinfo[[length(intinfo) + 1]] <- list(facteurs = as.list(it), colonnes = as.list(cols))
  }
  desc_col[desc_col == ""] <- "(libre)"
  tab_cols <- list(Colonne = seq_len(ncol_total), Niveaux = arr$niv, Affectation = desc_col)
  titre <- sprintf("Plan de Taguchi — tableau %s (%d essais)", arr$nom, n)
  lg <- c(titre, sprintf("Facteurs : %d — essais : %d.", k, n),
          paste0("Niveaux : ", paste(sprintf("%s (%s)", noms, vapply(niveaux, paste, "", collapse = " / ")), collapse = " ; ")))
  if (any(aff$fictif)) lg <- c(lg, paste0("Niveau fictif : ", paste(noms[aff$fictif], collapse = ", "),
                                          " (facteur à deux niveaux placé dans une colonne à trois niveaux ; le troisième niveau reprend le premier, ce qui rompt légèrement l'orthogonalité)."))
  if (length(inter)) lg <- c(lg, "Interactions étudiées : " %+% paste(vapply(inter, function(it) paste0(noms[it[1]], " × ", noms[it[2]]), ""), collapse = " ; ") %+% ".")
  else if (is.null(arr$vecteurs) || arr$nom %in% c("L12", "L18")) lg <- c(lg, "Ce tableau répartit les interactions sur l'ensemble des colonnes : il n'est pas adapté à l'étude d'une interaction précise.")
  if (isTRUE(par$aleatoire)) lg <- c(lg, "Ordre d'exécution : aléatoire (l'ordre standard est conservé dans la colonne « Ordre_standard »).")
  meta <- list(famille = "taguchi", tableau = arr$nom, colonnes = as.list(aff$colonnes), fictif = as.list(aff$fictif),
               interactions = intinfo, niveaux = lapply(niveaux, as.list))
  list(titre = titre, lignes = lg, tables = list(Affectation_des_colonnes = tab_cols), plan = df, meta = meta)
}

# Tableaux disponibles pour l'interface
liste_tableaux <- function() {
  list(noms = names(ARRAYS), essais = vapply(ARRAYS, function(a) a$runs, 1), colonnes = vapply(ARRAYS, function(a) length(a$niv), 1))
}
