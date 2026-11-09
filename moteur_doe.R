# -*- coding: utf-8 -*-
# =============================================================================
#  moteur_doe.R — Analyse des plans d'expériences
#    modèle linéaire (facteurs codés), effets, ANOVA, réduction du modèle,
#    prévision, optimisation (désirabilité), plans de surface et de mélange.
#  Les fonctions renvoient des listes (tables, lignes de texte, figures) ; la
#  conversion en JSON est faite par cli_doe.R.
# =============================================================================

# ---------------------------------------------------------------------------
#  Mise en forme
# ---------------------------------------------------------------------------
fr <- function(x, d = 4) {
  vapply(x, function(v) {
    if (is.na(v)) return("—")
    if (v == 0) return("0")
    s <- format(signif(v, d), trim = TRUE, digits = d, scientific = FALSE)
    if (nchar(s) > 12) s <- format(signif(v, d), trim = TRUE, digits = d)
    gsub(".", ",", s, fixed = TRUE)
  }, "")
}
txt_p <- function(p) ifelse(is.na(p), "p-value non calculable", ifelse(p < 1e-4, "p-value < 0,0001", paste0("p-value = ", gsub(".", ",", formatC(p, format = "f", digits = 4), fixed = TRUE))))
fp <- function(p) ifelse(is.na(p), "—", ifelse(p < 1e-4, "< 0,0001", gsub(".", ",", formatC(p, format = "f", digits = 4), fixed = TRUE)))
sg <- function(v, d = 4) {   # coefficient signé pour une équation
  paste0(ifelse(v < 0, "− ", "+ "), fr(abs(v), d))
}
liste_valeurs <- function(v) as.list(v)    # colonne à types mélangés

# ---------------------------------------------------------------------------
#  Contexte : plan + réponses + codage des facteurs
# ---------------------------------------------------------------------------
charger_contexte <- function(dossier, par) {
  meta <- par$meta
  if (is.null(meta) || is.null(meta$facteurs)) stop("Aucun plan chargé.", call. = FALSE)
  plan <- utils::read.csv(file.path(dossier, "plan.csv"), colClasses = "character", check.names = FALSE,
                          na.strings = c("", "NA"), encoding = "UTF-8", stringsAsFactors = FALSE)
  fac <- meta$facteurs
  noms <- vapply(fac, function(f) as.character(f$nom), "")
  famille <- meta$famille
  reps <- as.character(unlist(meta$reponses))
  for (j in seq_along(fac)) {
    if (famille == "taguchi" || identical(fac[[j]]$type, "cat")) next
    plan[[noms[j]]] <- suppressWarnings(as.numeric(plan[[noms[j]]]))
  }
  for (r in reps) if (r %in% names(plan)) plan[[r]] <- suppressWarnings(as.numeric(plan[[r]]))
  ctx <- list(plan = plan, fac = fac, noms = noms, k = length(noms), famille = famille, reps = reps, meta = meta,
              ech = as.numeric(par$echelle %||% 1), alpha = as.numeric(par$alpha_test %||% 0.05), dossier = dossier)
  if (famille == "melange") {
    ctx$L <- vapply(fac, function(f) as.numeric(f$bas %||% 0), 1)
    ctx$U <- vapply(fac, function(f) as.numeric(f$haut %||% 1), 1)
    ctx$sL <- sum(ctx$L)
    ctx$pseudo <- isTRUE(par$pseudo %||% TRUE) && ctx$sL > 1e-12
  }
  ctx
}

coder_matrice <- function(ctx, d) {
  X <- matrix(0, nrow(d), ctx$k, dimnames = list(NULL, ctx$noms))
  for (j in seq_len(ctx$k)) {
    f <- ctx$fac[[j]]; v <- d[[ctx$noms[j]]]
    if (ctx$famille == "melange") {
      v <- as.numeric(v)
      X[, j] <- if (ctx$pseudo) (v - ctx$L[j]) / (1 - ctx$sL) else v
    } else if (identical(f$type, "cat")) {
      lv <- as.character(unlist(f$niveaux))
      X[, j] <- ifelse(as.character(v) == lv[1], -1, 1)
    } else X[, j] <- (as.numeric(v) - (f$bas + f$haut) / 2) / ((f$haut - f$bas) / 2)
  }
  X
}

decoder_reel <- function(ctx, X) {
  X <- as.matrix(X)
  d <- list()
  for (j in seq_len(ctx$k)) {
    f <- ctx$fac[[j]]
    if (ctx$famille == "melange") d[[ctx$noms[j]]] <- if (ctx$pseudo) ctx$L[j] + X[, j] * (1 - ctx$sL) else X[, j]
    else if (identical(f$type, "cat")) d[[ctx$noms[j]]] <- ifelse(X[, j] < 0, unlist(f$niveaux)[1], unlist(f$niveaux)[2])
    else d[[ctx$noms[j]]] <- (f$bas + f$haut) / 2 + X[, j] * (f$haut - f$bas) / 2
  }
  as.data.frame(d, check.names = FALSE, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
#  Termes du modèle
# ---------------------------------------------------------------------------
nom_groupe_ordre <- function(o) if (o == 1) "Effets principaux" else sprintf("Interactions à %d facteurs", o)

generer_termes <- function(ctx, modele, courbure = FALSE) {
  k <- ctx$k; noms <- ctx$noms
  termes <- list()
  mk <- function(expo, label, groupe, fn = NULL, special = NULL) list(expo = expo, label = label, groupe = groupe, fn = fn, special = special)
  modele <- modele %||% (if (ctx$famille == "melange") "lineaire" else "ordre2")
  if (ctx$famille == "melange") {
    for (i in seq_len(k)) { e <- integer(k); e[i] <- 1L; termes[[length(termes) + 1]] <- mk(e, noms[i], "Mélange linéaire") }
    if (modele %in% c("quadratique", "cubique_special", "cubique") && k >= 2) {
      for (c in utils::combn(k, 2, simplify = FALSE)) { e <- integer(k); e[c] <- 1L
        termes[[length(termes) + 1]] <- mk(e, paste(noms[c], collapse = "*"), "Quadratique") }
    }
    if (modele %in% c("cubique_special", "cubique") && k >= 3) {
      for (c in utils::combn(k, 3, simplify = FALSE)) { e <- integer(k); e[c] <- 1L
        termes[[length(termes) + 1]] <- mk(e, paste(noms[c], collapse = "*"), "Cubique spécial") }
    }
    if (modele == "cubique" && k >= 2) {
      for (c in utils::combn(k, 2, simplify = FALSE)) {
        local({ a <- c[1]; b <- c[2]; e <- integer(k); e[c] <- 1L
          termes[[length(termes) + 1]] <<- mk(e, sprintf("%s*%s*(%s−%s)", noms[a], noms[b], noms[a], noms[b]), "Cubique complet",
                                              fn = function(X) X[, a] * X[, b] * (X[, a] - X[, b])) })
      }
    }
    return(termes)
  }
  ordre <- switch(modele, lineaire = 1, interactions = 2, quadratique = 2, ordre1 = 1, ordre2 = 2, ordre3 = 3, ordre4 = 4, ordre5 = 5, 2)
  ordre <- min(ordre, k)
  for (o in seq_len(ordre)) for (c in utils::combn(k, o, simplify = FALSE)) {
    e <- integer(k); e[c] <- 1L
    termes[[length(termes) + 1]] <- mk(e, paste(noms[c], collapse = "*"), nom_groupe_ordre(o))
  }
  if (modele == "quadratique") for (i in seq_len(k)) { e <- integer(k); e[i] <- 2L
    termes[[length(termes) + 1]] <- mk(e, paste0(noms[i], "*", noms[i]), "Termes quadratiques") }
  if (courbure) termes[[length(termes) + 1]] <- mk(integer(k), "Courbure (points centraux)", "Courbure", special = "courbure")
  termes
}

matrice_termes <- function(termes, X, intercept = TRUE) {
  n <- nrow(X)
  M <- matrix(0, n, length(termes))
  for (i in seq_along(termes)) {
    t <- termes[[i]]
    M[, i] <- if (identical(t$special, "courbure")) as.numeric(rowSums(abs(X)) < 1e-9)
              else if (!is.null(t$fn)) t$fn(X)
              else { idx <- which(t$expo > 0); if (!length(idx)) rep(1, n) else apply(X[, idx, drop = FALSE]^matrix(t$expo[idx], n, length(idx), byrow = TRUE), 1, prod) }
  }
  if (intercept) M <- cbind(1, M)
  M
}

# t est « exigé » par un autre terme retenu : divisibilité des monômes (hiérarchie)
est_requis <- function(termes, i, retenus, famille) {
  t <- termes[[i]]
  if (famille == "melange") return(identical(t$groupe, "Mélange linéaire"))
  if (identical(t$special, "courbure")) return(FALSE)
  for (j in setdiff(retenus, i)) {
    u <- termes[[j]]
    if (identical(u$special, "courbure")) next
    if (all(u$expo >= t$expo) && any(u$expo > t$expo)) return(TRUE)
  }
  FALSE
}

# ---------------------------------------------------------------------------
#  Moindres carrés
# ---------------------------------------------------------------------------
ols <- function(M, y) {
  n <- nrow(M); p <- ncol(M)
  q <- qr(M)
  coef <- as.vector(qr.coef(q, y))
  fitted <- as.vector(M %*% coef)
  res <- y - fitted
  df <- n - p; rss <- sum(res^2)
  s2 <- if (df > 0) rss / df else NA_real_
  XtXi <- chol2inv(qr.R(q))
  h <- rowSums(qr.Q(q)^2)
  list(coef = coef, fitted = fitted, res = res, df = df, rss = rss, s2 = s2, XtXi = XtXi, h = h, n = n, p = p,
       se = if (df > 0) sqrt(s2 * diag(XtXi)) else rep(NA_real_, p))
}

rss_sans <- function(M, y, colonnes) {
  keep <- setdiff(seq_len(ncol(M)), colonnes)
  if (!length(keep)) return(sum(y^2))
  ols(M[, keep, drop = FALSE], y)$rss
}

# Ajustement d'un jeu de termes (colonnes dépendantes écartées, par ordre de priorité)
ajuster <- function(termes, X, y, intercept = TRUE) {
  ok <- !is.na(y)
  M <- matrice_termes(termes, X, intercept)
  Mo <- M[ok, , drop = FALSE]
  q <- qr(Mo, tol = 1e-7)
  gardes <- sort(q$pivot[seq_len(q$rank)])
  off <- if (intercept) 1L else 0L
  idx_termes <- gardes[gardes > off] - off
  perdus <- setdiff(seq_along(termes), idx_termes)
  Mk <- Mo[, gardes, drop = FALSE]
  mod <- ols(Mk, y[ok])
  mod$termes <- termes[idx_termes]
  mod$perdus <- termes[perdus]
  mod$intercept <- intercept
  mod$ok <- ok
  mod$M <- Mk
  mod$y <- y[ok]
  mod$tval <- mod$coef / mod$se
  mod$pval <- if (mod$df > 0) 2 * stats::pt(-abs(mod$tval), mod$df) else rep(NA_real_, mod$p)
  mod
}

# Prévision avec intervalles
predire <- function(mod, Xnew, niveau = 0.95) {
  Mn <- matrice_termes(mod$termes, Xnew, mod$intercept)
  fit <- as.vector(Mn %*% mod$coef)
  v <- rowSums((Mn %*% mod$XtXi) * Mn)
  if (mod$df > 0) {
    s <- sqrt(mod$s2); tq <- stats::qt(1 - (1 - niveau) / 2, mod$df)
    se <- s * sqrt(v); sep <- s * sqrt(1 + v)
    list(fit = fit, se = se, ic_inf = fit - tq * se, ic_sup = fit + tq * se, ip_inf = fit - tq * sep, ip_sup = fit + tq * sep)
  } else list(fit = fit, se = rep(NA_real_, length(fit)), ic_inf = NA * fit, ic_sup = NA * fit, ip_inf = NA * fit, ip_sup = NA * fit)
}

pred_simple <- function(mod, Xnew) as.vector(matrice_termes(mod$termes, Xnew, mod$intercept) %*% mod$coef)

# Méthode de Lenth (plans sans degré de liberté pour l'erreur)
lenth <- function(effets, alpha = 0.05) {
  e <- effets[is.finite(effets)]
  m <- length(e)
  s0 <- 1.5 * stats::median(abs(e))
  pse <- 1.5 * stats::median(abs(e)[abs(e) < 2.5 * s0])
  if (!is.finite(pse) || pse <= 0) pse <- s0
  d <- m / 3
  me <- stats::qt(1 - alpha / 2, d) * pse
  gam <- (1 + (1 - alpha)^(1 / m)) / 2
  sme <- stats::qt(gam, d) * pse
  list(pse = pse, me = me, sme = sme, d = d,
       t = effets / pse, p = 2 * stats::pt(-abs(effets / pse), d))
}

# ---------------------------------------------------------------------------
#  Réduction du modèle (retrait pas à pas, hiérarchie respectée)
# ---------------------------------------------------------------------------
reduire_modele <- function(ctx, termes, X, y, alpha_out, intercept = TRUE) {
  etapes <- list(); note <- character(0)
  courant <- termes
  repeat {
    mod <- ajuster(courant, X, y, intercept)
    courant <- mod$termes                               # termes aliasés déjà écartés
    if (mod$df <= 0) {
      # pas d'erreur estimable : sélection par la méthode de Lenth (effets à deux niveaux seulement)
      if (ctx$famille %in% c("criblage", "factoriel")) {
        off <- if (intercept) 1L else 0L
        ef <- 2 * mod$coef[(off + 1):mod$p]
        est_court <- vapply(courant, function(t) identical(t$special, "courbure"), TRUE)
        idx_nc <- which(!est_court)
        L <- lenth(ef[idx_nc], alpha_out)
        sig <- rep(FALSE, length(courant)); sig[idx_nc] <- abs(ef[idx_nc]) > L$me
        pv <- rep(NA_real_, length(courant)); pv[idx_nc] <- L$p
        garder <- which(sig)
        repeat {
          ajout <- Filter(function(j) est_requis(courant, j, garder, ctx$famille), setdiff(seq_along(courant), garder))
          if (!length(ajout)) break
          garder <- sort(c(garder, ajout))
        }
        for (j in setdiff(seq_along(courant), garder)) etapes[[length(etapes) + 1]] <- list(terme = courant[[j]]$label, p = pv[j])
        courant <- courant[garder]
        note <- "Aucun degré de liberté pour l'erreur : les termes ont été retenus d'après la marge d'erreur de Lenth."
      } else note <- "Aucun degré de liberté pour l'erreur : le modèle ne peut pas être réduit."
      break
    }
    off <- if (intercept) 1L else 0L
    p <- mod$pval[(off + 1):mod$p]
    retenus <- seq_along(courant)
    cand <- which(vapply(retenus, function(i) !est_requis(courant, i, retenus, ctx$famille), TRUE) & p > alpha_out)
    if (!length(cand)) break
    pire <- cand[which.max(p[cand])]
    etapes[[length(etapes) + 1]] <- list(terme = courant[[pire]]$label, p = p[pire])
    courant <- courant[-pire]
    if (!length(courant)) break
  }
  list(termes = courant, etapes = etapes, note = note)
}

# Modèle prêt à l'emploi pour une réponse
preparer_modele <- function(ctx, reponse, modele, reduire, alpha_out, courbure_auto = TRUE) {
  if (!(reponse %in% names(ctx$plan))) stop(sprintf("Réponse « %s » introuvable.", reponse), call. = FALSE)
  y <- ctx$plan[[reponse]]
  if (sum(!is.na(y)) < 3) stop(sprintf("La réponse « %s » compte moins de 3 valeurs renseignées.", reponse), call. = FALSE)
  X <- coder_matrice(ctx, ctx$plan)
  ok <- !is.na(y)
  if (any(is.na(X[ok, ]))) stop("Valeurs de facteurs manquantes ou invalides dans le plan.", call. = FALSE)
  courbure <- courbure_auto && ctx$famille %in% c("criblage", "factoriel") && any(rowSums(abs(X)) < 1e-9 & ok) && any(rowSums(abs(X)) > 1e-9 & ok) &&
    !any(vapply(ctx$fac, function(f) identical(f$type, "cat"), TRUE))
  termes <- generer_termes(ctx, modele, courbure)
  intercept <- ctx$famille != "melange"
  red <- NULL
  complet <- ajuster(termes, X, y, intercept)
  if (isTRUE(reduire)) {
    red <- reduire_modele(ctx, termes, X, y, alpha_out, intercept)
    mod <- ajuster(red$termes, X, y, intercept)
  } else mod <- complet
  list(mod = mod, complet = complet, red = red, X = X, y = y, intercept = intercept, termes_initiaux = termes, reponse = reponse)
}

# ---------------------------------------------------------------------------
#  Équations
# ---------------------------------------------------------------------------
equation_texte <- function(reponse, coefs, labels, intercept = TRUE, d = 4) {
  coefs[abs(coefs) < 1e-9 * max(abs(coefs), 1e-300)] <- 0
  parts <- character(0)
  if (intercept) parts <- sub("^-", "−", fr(coefs[1], d))
  off <- if (intercept) 1L else 0L
  for (i in seq_along(labels)) {
    c <- coefs[off + i]
    t <- sprintf("%s %s", sg(c, d), labels[i])
    parts <- c(parts, if (length(parts) == 0) sub("^\\+ ", "", t) else t)
  }
  paste0(reponse, " = ", paste(parts, collapse = " "))
}

# Équation en unités réelles : on ré-ajuste les mêmes monômes sur les valeurs réelles (modèle hiérarchique)
equation_reelle <- function(ctx, B) {
  mod <- B$mod
  if (any(vapply(ctx$fac, function(f) identical(f$type, "cat"), TRUE))) return(NULL)
  if (any(vapply(mod$termes, function(t) !is.null(t$fn) || identical(t$special, "courbure"), TRUE))) return(NULL)
  d <- decoder_reel(ctx, B$X[B$mod$ok, , drop = FALSE])
  Xr <- as.matrix(d); colnames(Xr) <- ctx$noms
  Mr <- matrice_termes(mod$termes, Xr, mod$intercept)
  fit <- ols(Mr, mod$fitted)
  if (any(!is.finite(fit$coef)) || sqrt(fit$rss / max(1, length(mod$fitted))) > 1e-6 * (1 + sd(mod$fitted))) return(NULL)
  fit$coef
}

# ---------------------------------------------------------------------------
#  Sorties : tableaux
# ---------------------------------------------------------------------------
tableau_resume <- function(mod, y, reponse) {
  n <- mod$n; tss <- sum((y[mod$ok] - mean(y[mod$ok]))^2)
  r2 <- 1 - mod$rss / tss
  p <- mod$p
  r2aj <- if (mod$df > 0) 1 - (mod$rss / mod$df) / (tss / (n - 1)) else NA
  press <- if (mod$df > 0 && all(mod$h < 1 - 1e-8)) sum((mod$res / (1 - mod$h))^2) else NA
  r2pr <- if (is.na(press)) NA else 1 - press / tss
  list(r2 = r2, r2aj = r2aj, r2pr = r2pr, s = if (mod$df > 0) sqrt(mod$s2) else NA, tss = tss)
}

nom_ctl <- function(x) x

tableau_anova <- function(ctx, mod, y, groupes_ordre = NULL) {
  n <- mod$n; tss <- sum((mod$y - mean(mod$y))^2)
  off <- if (mod$intercept) 1L else 0L
  nt <- length(mod$termes)
  sources <- character(0); ddl <- numeric(0); ss <- numeric(0)
  add <- function(nom, d, s) { sources <<- c(sources, nom); ddl <<- c(ddl, d); ss <<- c(ss, s) }
  add("Modèle", mod$p - off, tss - mod$rss)
  if (nt > 0) {
    grp <- vapply(mod$termes, function(t) t$groupe, "")
    ordre_g <- unique(grp)
    for (g in ordre_g) {
      ix <- which(grp == g)
      if (identical(g, "Courbure")) next
      ssg <- rss_sans(mod$M, mod$y, off + ix) - mod$rss
      add(g, length(ix), ssg)
      for (i in ix) add(paste0("   ", mod$termes[[i]]$label), 1, rss_sans(mod$M, mod$y, off + i) - mod$rss)
    }
    if ("Courbure" %in% grp) { i <- which(grp == "Courbure"); add("Courbure", 1, rss_sans(mod$M, mod$y, off + i) - mod$rss) }
  }
  if (mod$df > 0) {
    add("Erreur résiduelle", mod$df, mod$rss)
    Xo <- ctx$.Xutil
    if (!is.null(Xo)) {
      cle <- apply(round(Xo, 6), 1, paste, collapse = "|")
      grp <- split(seq_along(cle), cle)
      dfpe <- sum(lengths(grp) - 1)
      if (dfpe > 0 && mod$df - dfpe > 0) {
        sspe <- sum(vapply(grp, function(ix) sum((mod$y[ix] - mean(mod$y[ix]))^2), 0))
        add("   Manque d'ajustement", mod$df - dfpe, mod$rss - sspe)
        add("   Erreur pure", dfpe, sspe)
      }
    }
  }
  add("Total", n - 1, tss)
  ms <- ss / ddl
  F <- rep(NA_real_, length(ss)); p <- rep(NA_real_, length(ss))
  if (mod$df > 0) {
    mse <- mod$rss / mod$df
    for (i in seq_along(sources)) {
      if (sources[i] %in% c("Erreur résiduelle", "Total", "   Erreur pure", "   Manque d'ajustement")) next
      F[i] <- ms[i] / mse; p[i] <- stats::pf(F[i], ddl[i], mod$df, lower.tail = FALSE)
    }
    i_lof <- match("   Manque d'ajustement", sources); i_pe <- match("   Erreur pure", sources)
    if (!is.na(i_lof)) { F[i_lof] <- ms[i_lof] / ms[i_pe]; p[i_lof] <- stats::pf(F[i_lof], ddl[i_lof], ddl[i_pe], lower.tail = FALSE) }
  }
  list(Source = sources, ddl = ddl, Somme_carrés = ss, Carré_moyen = ms, F = F, p_value = p)
}

tableau_coefficients <- function(ctx, mod, deux_niveaux) {
  off <- if (mod$intercept) 1L else 0L
  labels <- vapply(mod$termes, function(t) t$label, "")
  idx <- seq_len(mod$p)
  noms <- c(if (mod$intercept) "Constante", labels)
  coef <- mod$coef
  effet <- if (deux_niveaux) c(if (mod$intercept) NA, ifelse(vapply(mod$termes, function(t) identical(t$special, "courbure"), TRUE), NA, 2 * coef[(off + 1):mod$p])) else NULL
  vif <- rep(NA_real_, mod$p)
  nt <- mod$p - off
  if (nt >= 2 && ctx$famille != "melange") {
    Mt <- mod$M[, (off + 1):mod$p, drop = FALSE]
    cm <- suppressWarnings(stats::cor(Mt)); cm[is.na(cm)] <- 0; diag(cm) <- 1
    iv <- tryCatch(diag(solve(cm)), error = function(e) rep(NA_real_, nt))
    vif[(off + 1):mod$p] <- iv
  } else if (nt == 1 && ctx$famille != "melange") vif[mod$p] <- 1
  tab <- list(Terme = noms)
  if (!is.null(effet)) tab$Effet <- effet
  tab$Coefficient <- coef; tab$Erreur_type <- mod$se; tab$t <- mod$tval; tab$p_value <- mod$pval
  if (ctx$famille != "melange") tab$VIF <- vif
  tab
}

# ---------------------------------------------------------------------------
#  Analyse d'un plan factoriel / criblage / surface / mélange
# ---------------------------------------------------------------------------
analyser_modele <- function(ctx, par) {
  reponse <- par$reponse
  alpha <- ctx$alpha
  modele <- par$modele %||% (if (ctx$famille == "melange") "quadratique" else if (ctx$famille == "criblage") "ordre1" else if (ctx$famille == "surface") "quadratique" else "ordre2")
  reduire <- isTRUE(par$reduire); alpha_out <- as.numeric(par$alpha_sortie %||% 0.05)
  B <- preparer_modele(ctx, reponse, modele, reduire, alpha_out)
  mod <- B$mod; y <- B$y
  ctx$.Xutil <- B$X[mod$ok, , drop = FALSE]
  deux <- ctx$famille %in% c("criblage", "factoriel")
  off <- if (mod$intercept) 1L else 0L
  res <- tableau_resume(mod, y, reponse)
  tables <- list(); lignes <- character(0); avert <- character(0); figures <- character(0)
  fam_nom <- c(criblage = "Plan de criblage", factoriel = "Plan factoriel", surface = "Plan de surface de réponse", melange = "Plan de mélange")[[ctx$famille]]
  titre <- sprintf("%s — analyse de « %s »", fam_nom, reponse)
  labels <- vapply(mod$termes, function(t) t$label, "")
  n <- mod$n

  # -- résumé
  lignes <- c(lignes, sprintf("Réponse : %s — %d essais renseignés sur %d.", reponse, n, nrow(ctx$plan)))
  mdl_txt <- c(ordre1 = "effets principaux", ordre2 = "effets principaux et interactions à 2 facteurs", ordre3 = "effets et interactions jusqu'à l'ordre 3",
               ordre4 = "effets et interactions jusqu'à l'ordre 4", ordre5 = "effets et interactions jusqu'à l'ordre 5", lineaire = "modèle linéaire",
               interactions = "linéaire et interactions à 2 facteurs", quadratique = if (ctx$famille == "melange") "quadratique de Scheffé" else "quadratique complet",
               cubique_special = "cubique spécial de Scheffé", cubique = "cubique complet de Scheffé")
  lignes <- c(lignes, sprintf("Modèle demandé : %s%s.", mdl_txt[[modele]], if (reduire) sprintf(" — réduit pas à pas (retrait des termes dont la p-value dépasse %s)", fr(alpha_out, 3)) else ""))
  if (ctx$famille == "melange") lignes <- c(lignes, if (ctx$pseudo) "Proportions codées en pseudo-composants (bornes inférieures retranchées) ; modèle de Scheffé sans constante." else "Proportions réelles ; modèle de Scheffé sans constante.")
  if (length(mod$perdus)) avert <- c(avert, sprintf("Termes non estimables (confondus avec d'autres termes du modèle) : %s.", paste(vapply(mod$perdus, function(t) t$label, ""), collapse = ", ")))
  if (ctx$famille == "factoriel" && !is.null(ctx$meta$resolution) && !(ctx$meta$resolution %in% c("complet"))) avert <- c(avert, "Plan fractionnaire : les effets estimés sont confondus avec leurs alias (voir la structure des alias du plan).")
  if (ctx$famille == "criblage") avert <- c(avert, "Plan de criblage : à utiliser pour retenir les facteurs influents ; les interactions éventuelles sont confondues avec les effets principaux.")

  rm <- list(Indicateur = c("Écart-type résiduel S", "R²", "R² ajusté", "R² prévu", "Degrés de liberté de l'erreur", "Nombre de termes dans le modèle"),
             Valeur = liste_valeurs(c(res$s, res$r2, res$r2aj, res$r2pr, mod$df, length(mod$termes))))
  tables[["Résumé_du_modèle"]] <- rm

  # -- Lenth si pas de degré de liberté
  lth <- NULL
  if (deux && mod$df <= 0) {
    est_c <- vapply(mod$termes, function(t) identical(t$special, "courbure"), TRUE)
    ef <- 2 * mod$coef[(off + 1):mod$p]
    lth <- lenth(ef[!est_c], alpha)
    sig <- abs(ef[!est_c]) > lth$me
    lab_l <- labels[!est_c]
    tables[["Effets_de_Lenth"]] <- list(Terme = lab_l, Effet = ef[!est_c], t_Lenth = lth$t, p_value = lth$p, Significatif = ifelse(sig, "oui", "non"))
    lignes <- c(lignes, sprintf("Aucun degré de liberté pour l'erreur : significativité jugée par la méthode de Lenth (erreur type pseudo-estimée PSE = %s, marge d'erreur = %s).", fr(lth$pse), fr(lth$me)))
    signif_labels <- lab_l[sig]
  } else if (mod$df > 0) {
    p <- mod$pval[(off + 1):mod$p]
    signif_labels <- labels[!is.na(p) & p < alpha]
  } else signif_labels <- character(0)

  tables[["Coefficients"]] <- tableau_coefficients(ctx, mod, deux)
  if (mod$df >= 0) tables[["Analyse_de_la_variance"]] <- tableau_anova(ctx, mod, y)

  if (!is.null(B$red) && length(B$red$etapes)) {
    tables[["Étapes_de_la_réduction"]] <- list(Étape = seq_along(B$red$etapes), Terme_retiré = vapply(B$red$etapes, function(e) e$terme, ""), p_value = vapply(B$red$etapes, function(e) e$p, 0))
    lignes <- c(lignes, sprintf("Réduction du modèle : %d terme(s) retiré(s) ; R² ajusté du modèle complet %s, du modèle réduit %s.", length(B$red$etapes),
                                fr(tableau_resume(B$complet, y, reponse)$r2aj, 3), fr(res$r2aj, 3)))
  } else if (isTRUE(reduire)) lignes <- c(lignes, "Réduction du modèle : tous les termes sont déjà significatifs (ou exigés par la hiérarchie).")
  if (!is.null(B$red) && nzchar(paste(B$red$note, collapse = ""))) avert <- c(avert, B$red$note)

  # -- conclusion
  if (length(signif_labels)) {
    lignes <- c(lignes, "", sprintf("Termes significatifs au risque α = %s : %s.", fr(alpha, 3), if (length(signif_labels)) paste(signif_labels, collapse = ", ") else "aucun"))
    concl <- sprintf("Significatifs (α = %s) : %s ; R² ajusté = %s.", fr(alpha, 3), paste(signif_labels, collapse = ", "), fr(res$r2aj, 3))
  } else { lignes <- c(lignes, "", sprintf("Aucun terme n'est significatif au risque α = %s.", fr(alpha, 3))); concl <- sprintf("Aucun terme significatif (α = %s).", fr(alpha, 3)) }
  if ("Courbure (points centraux)" %in% labels && mod$df > 0) {
    pc <- mod$pval[off + which(labels == "Courbure (points centraux)")]
    lignes <- c(lignes, sprintf("Courbure (points centraux) : %s → %s.", txt_p(pc), if (pc < alpha) "courbure significative : un modèle du premier degré ne suffit pas, envisager un plan de surface" else "pas de courbure significative"))
  }
  if (!is.na(res$r2pr) && !is.na(res$r2aj) && res$r2aj - res$r2pr > 0.2) avert <- c(avert, "Le R² prévu est nettement inférieur au R² ajusté : le modèle risque d'être sur-ajusté ou de prévoir mal.")
  if (mod$df > 0 && mod$df < 3) avert <- c(avert, "Très peu de degrés de liberté pour l'erreur : les tests ont une puissance faible ; ajouter des répétitions ou des points centraux.")

  # -- équations
  eq <- equation_texte(reponse, mod$coef, labels, mod$intercept)
  lignes <- c(lignes, "", "Équation de régression en unités codées :", paste0("  ", eq))
  er <- equation_reelle(ctx, B)
  if (!is.null(er)) { lignes <- c(lignes, "Équation de régression en unités réelles :", paste0("  ", equation_texte(reponse, er, labels, mod$intercept, 5))) }
  if (ctx$famille != "melange") lignes <- c(lignes, "(unités codées : −1 au niveau bas, +1 au niveau haut, 0 au centre du domaine)")

  # -- observations inhabituelles
  if (mod$df > 0) {
    s <- sqrt(mod$s2); rs <- mod$res / (s * sqrt(pmax(1 - mod$h, 1e-12)))
    levier <- mod$h > 3 * mod$p / n
    mauv <- which(abs(rs) > 2 | levier)
    if (length(mauv)) {
      ess <- ctx$plan$Ordre_essai[mod$ok][mauv]
      tables[["Observations_inhabituelles"]] <- list(Essai = ess, Réponse = mod$y[mauv], Valeur_ajustée = mod$fitted[mauv], Résidu_standardisé = rs[mauv], Levier = mod$h[mauv],
                                                      Remarque = ifelse(abs(rs[mauv]) > 2, "résidu élevé", "levier élevé"))
    }
    if (n >= 3 && n <= 5000 && sd(mod$res) > 0) {
      sw <- tryCatch(stats::shapiro.test(mod$res)$p.value, error = function(e) NA)
      if (!is.na(sw)) lignes <- c(lignes, sprintf("Normalité des résidus (Shapiro-Wilk) : %s → %s.", txt_p(sw), if (sw < alpha) "hypothèse de normalité rejetée, interpréter les p-values avec prudence" else "normalité non rejetée"))
    }
  }

  # -- point stationnaire (surface)
  if (ctx$famille == "surface" && modele == "quadratique") {
    st <- point_stationnaire(ctx, B)
    if (!is.null(st)) {
      tables[["Point_stationnaire"]] <- st$table
      lignes <- c(lignes, "", st$texte)
    }
  }

  # -- figures
  ech <- ctx$ech
  if (ctx$famille %in% c("criblage", "factoriel")) {
    if (!is.null(lth)) {
      f <- fig_pareto(ctx$dossier, "01_pareto.png", lab_l, abs(ef[!est_c]), lth$me, "Diagramme de Pareto des effets (méthode de Lenth)", "Effet absolu", ech, ligne_lib = sprintf("Marge d'erreur = %s", fr(lth$me)))
      figures <- c(figures, f)
      figures <- c(figures, fig_daniel(ctx$dossier, "02_daniel.png", lab_l, ef[!est_c], abs(ef[!est_c]) > lth$me, ech))
    } else if (mod$df > 0) {
      tcrit <- stats::qt(1 - alpha / 2, mod$df)
      ef_t <- abs(mod$tval[(off + 1):mod$p])
      figures <- c(figures, fig_pareto(ctx$dossier, "01_pareto.png", labels, ef_t, tcrit, "Diagramme de Pareto des effets standardisés", "|t| de l'effet", ech,
                                        ligne_lib = sprintf("t critique = %s (α = %s)", fr(tcrit), fr(alpha, 3))))
      est_c <- vapply(mod$termes, function(t) identical(t$special, "courbure"), TRUE)
      ef <- 2 * mod$coef[(off + 1):mod$p]
      figures <- c(figures, fig_daniel(ctx$dossier, "02_daniel.png", labels[!est_c], ef[!est_c], abs(mod$tval[(off + 1):mod$p])[!est_c] > tcrit, ech))
    }
    figures <- c(figures, fig_effets_principaux(ctx, y, mod, B, "03_effets_principaux.png"))
    if (ctx$famille == "factoriel" || (ctx$famille == "criblage" && any(vapply(mod$termes, function(t) sum(t$expo) == 2 && is.null(t$special), TRUE)))) {
      fi <- fig_interactions(ctx, y, mod, "04_interactions.png", signif_labels)
      if (!is.null(fi)) figures <- c(figures, fi)
    }
    if (ctx$k == 3 && ctx$famille == "factoriel") figures <- c(figures, fig_cube(ctx, y, "05_cube.png"))
  } else if (ctx$famille == "surface") {
    if (mod$df > 0) {
      tcrit <- stats::qt(1 - alpha / 2, mod$df)
      figures <- c(figures, fig_pareto(ctx$dossier, "01_pareto.png", labels, abs(mod$tval[(off + 1):mod$p]), tcrit, "Diagramme de Pareto des effets standardisés", "|t| du coefficient", ech,
                                        ligne_lib = sprintf("t critique = %s (α = %s)", fr(tcrit), fr(alpha, 3))))
    }
    figures <- c(figures, fig_effets_modele(ctx, B, "03_effets_principaux.png"))
  } else if (ctx$famille == "melange") {
    figures <- c(figures, fig_trace_melange(ctx, B, "01_trace.png"))
  }
  if (mod$df > 0) figures <- c(figures, fig_residus(ctx$dossier, "09_residus.png", mod, ctx$plan$Ordre_essai[mod$ok], reponse, ech))

  list(titre = titre, lignes = lignes, tables = tables, figures = figures, avertissements = avert, conclusion = concl, tables_libres = list("Coefficients"))
}

# Point stationnaire et analyse canonique du modèle quadratique
point_stationnaire <- function(ctx, B) {
  mod <- B$mod; k <- ctx$k
  labels <- vapply(mod$termes, function(t) t$label, "")
  b <- numeric(k); Bm <- matrix(0, k, k)
  off <- 1L
  for (i in seq_along(mod$termes)) {
    t <- mod$termes[[i]]; c <- mod$coef[off + i]
    ix <- which(t$expo > 0)
    if (sum(t$expo) == 1) b[ix] <- c
    else if (length(ix) == 1 && t$expo[ix] == 2) Bm[ix, ix] <- c
    else if (length(ix) == 2 && all(t$expo[ix] == 1)) { Bm[ix[1], ix[2]] <- c / 2; Bm[ix[2], ix[1]] <- c / 2 }
  }
  if (all(Bm == 0)) return(NULL)
  sol <- tryCatch(solve(2 * Bm, -b), error = function(e) NULL)
  if (is.null(sol)) return(NULL)
  ev <- eigen(Bm, symmetric = TRUE)$values
  nature <- if (all(ev < 0)) "maximum" else if (all(ev > 0)) "minimum" else "point col (ni maximum ni minimum)"
  Xs <- matrix(sol, 1); colnames(Xs) <- ctx$noms
  ys <- pred_simple(mod, Xs)
  reel <- decoder_reel(ctx, Xs)
  rng <- apply(B$X, 2, range)
  dans <- all(sol >= rng[1, ] - 1e-9 & sol <= rng[2, ] + 1e-9)
  tab <- list(Facteur = ctx$noms, Valeur_codée = as.numeric(sol), Valeur_réelle = as.numeric(unlist(reel)))
  txt <- c(sprintf("Analyse canonique : le point stationnaire est un %s (valeurs propres : %s).", nature, paste(fr(ev, 3), collapse = " ; ")),
           sprintf("Réponse prévue au point stationnaire : %s%s.", fr(ys, 5), if (dans) "" else " — attention : ce point est hors du domaine expérimenté, extrapolation."))
  list(table = tab, texte = txt)
}

# ---------------------------------------------------------------------------
#  Domaine de recherche, désirabilité et optimisation
# ---------------------------------------------------------------------------
domaine_facteurs <- function(ctx, X) {
  # bornes codées observées pour chaque facteur
  rng <- apply(X, 2, range)
  rng
}

desirabilite <- function(y, but, inf, cible, sup, poids) {
  d <- rep(0, length(y)); w <- poids %||% 1
  if (but == "max") { d <- ifelse(y >= cible, 1, ifelse(y <= inf, 0, ((y - inf) / (cible - inf))^w)) }
  else if (but == "min") { d <- ifelse(y <= cible, 1, ifelse(y >= sup, 0, ((sup - y) / (sup - cible))^w)) }
  else { d <- ifelse(y < inf | y > sup, 0, ifelse(y <= cible, ((y - inf) / (cible - inf))^w, ((sup - y) / (sup - cible))^w)) }
  d[is.na(d)] <- 0
  d
}

echantillon_melange <- function(ctx, n) {
  V <- sommets_extremes(ctx$L, ctx$U)
  if (is.null(V)) stop("Domaine du mélange vide.", call. = FALSE)
  nv <- nrow(V)
  W <- matrix(stats::rexp(n * nv, rate = 1), n, nv)
  W <- W / rowSums(W)
  P <- W %*% V
  rbind(V, colMeans(V), P)
}

optimiser_modeles <- function(ctx, modeles, buts, par) {
  set.seed(12345)
  k <- ctx$k
  nr <- length(modeles)
  evalue <- function(X) {
    Y <- vapply(modeles, function(B) pred_simple(B$mod, X), numeric(nrow(X)))
    if (is.null(dim(Y))) Y <- matrix(Y, nrow = nrow(X))
    D <- matrix(0, nrow(X), nr)
    for (r in seq_len(nr)) D[, r] <- desirabilite(Y[, r], buts[[r]]$but, buts[[r]]$inf, buts[[r]]$cible, buts[[r]]$sup, buts[[r]]$poids)
    imp <- vapply(buts, function(b) b$importance, 1)
    Dg <- apply(D, 1, function(d) if (any(d <= 0)) 0 else exp(sum(imp * log(d)) / sum(imp)))
    list(Y = Y, D = D, Dg = Dg)
  }
  X0 <- modeles[[1]]$X
  if (ctx$famille == "melange") {
    cand <- echantillon_melange(ctx, 60000)
    if (ctx$pseudo) cand <- sweep(sweep(cand, 2, ctx$L, "-"), 1, 1 - ctx$sL, "/")
    colnames(cand) <- ctx$noms
  } else {
    rng <- domaine_facteurs(ctx, X0)
    cat_f <- vapply(ctx$fac, function(f) identical(f$type, "cat"), TRUE)
    n <- 40000
    cand <- matrix(0, n, k)
    for (j in seq_len(k)) cand[, j] <- if (cat_f[j]) sample(c(-1, 1), n, TRUE) else stats::runif(n, rng[1, j], rng[2, j])
    # points de la grille de base : sommets, centre
    cand <- rbind(cand, as.matrix(expand.grid(lapply(seq_len(k), function(j) if (cat_f[j]) c(-1, 1) else c(rng[1, j], mean(rng[, j]), rng[2, j])))))
    colnames(cand) <- ctx$noms
  }
  ev <- evalue(cand)
  crit <- ev$Dg
  repli <- FALSE
  if (max(crit) <= 0) {
    crit <- rowMeans(ev$D); repli <- TRUE
    if (max(crit) <= 0) stop("Aucune solution n'atteint les limites acceptables fixées pour les réponses : élargir les limites ou revoir les objectifs.", call. = FALSE)
  }
  ordre <- order(crit, decreasing = TRUE)
  meilleur <- cand[ordre[1], ]
  # affinage local
  if (ctx$famille == "melange") {
    x <- meilleur; fx <- crit[ordre[1]]
    Lc <- if (ctx$pseudo) rep(0, k) else ctx$L
    Uc <- if (ctx$pseudo) (ctx$U - ctx$L) / (1 - ctx$sL) else ctx$U
    pas <- 0.05
    for (it in 1:4000) {
      i <- sample(k, 1); j <- sample(setdiff(seq_len(k), i), 1); t <- stats::runif(1, 0, pas)
      xn <- x; xn[i] <- xn[i] - t; xn[j] <- xn[j] + t
      if (xn[i] < Lc[i] - 1e-12 || xn[j] > Uc[j] + 1e-12) next
      en <- evalue(matrix(xn, 1, dimnames = list(NULL, ctx$noms)))
      fn <- if (repli) mean(en$D) else en$Dg
      if (fn > fx) { x <- xn; fx <- fn }
      if (it %% 1000 == 0) pas <- pas / 2
    }
    meilleur <- x
  } else {
    cat_f <- vapply(ctx$fac, function(f) identical(f$type, "cat"), TRUE)
    obj <- function(z) { zz <- meilleur; zz[!cat_f] <- z; e <- evalue(matrix(zz, 1, dimnames = list(NULL, ctx$noms))); -(if (repli) mean(e$D) else e$Dg) }
    if (any(!cat_f)) {
      op <- tryCatch(stats::optim(meilleur[!cat_f], obj, method = "L-BFGS-B", lower = rng[1, !cat_f], upper = rng[2, !cat_f]), error = function(e) NULL)
      if (!is.null(op) && -op$value >= crit[ordre[1]]) meilleur[!cat_f] <- op$par
    }
  }
  Xo <- matrix(meilleur, 1, dimnames = list(NULL, ctx$noms))
  eo <- evalue(Xo)
  list(Xo = Xo, ev = eo, repli = repli, evalue = evalue, X0 = X0)
}

optimiser_doe <- function(ctx, par) {
  reps <- par$reponses
  if (!length(reps)) stop("Aucune réponse à optimiser.", call. = FALSE)
  modele <- par$modele; reduire <- isTRUE(par$reduire); alpha_out <- as.numeric(par$alpha_sortie %||% 0.05)
  modeles <- list(); buts <- list(); avert <- character(0)
  for (r in reps) {
    B <- preparer_modele(ctx, r$nom, modele, reduire, alpha_out)
    y <- B$y[!is.na(B$y)]
    but <- r$but %||% "max"
    inf <- r$inf; cible <- r$cible; sup <- r$sup
    mn <- min(y); mx <- max(y)
    if (but == "max") { inf <- inf %||% mn; cible <- cible %||% mx; sup <- NA }
    else if (but == "min") { cible <- cible %||% mn; sup <- sup %||% mx; inf <- NA }
    else { cible <- cible %||% mean(c(mn, mx)); inf <- inf %||% mn; sup <- sup %||% mx }
    if (but == "max" && inf >= cible) stop(sprintf("Réponse « %s » : la valeur minimale acceptable doit être inférieure à la cible.", r$nom), call. = FALSE)
    if (but == "min" && sup <= cible) stop(sprintf("Réponse « %s » : la valeur maximale acceptable doit être supérieure à la cible.", r$nom), call. = FALSE)
    if (but == "cible" && !(inf < cible && cible < sup)) stop(sprintf("Réponse « %s » : il faut inf. < cible < sup.", r$nom), call. = FALSE)
    if (B$mod$df <= 0) avert <- c(avert, sprintf("Réponse « %s » : aucun degré de liberté pour l'erreur, le modèle est exact sur les essais mais non validé.", r$nom))
    modeles[[length(modeles) + 1]] <- B
    buts[[length(buts) + 1]] <- list(nom = r$nom, but = but, inf = inf, cible = cible, sup = sup, poids = as.numeric(r$poids %||% 1), importance = as.numeric(r$importance %||% 1))
  }
  op <- optimiser_modeles(ctx, modeles, buts, par)
  reel <- decoder_reel(ctx, op$Xo)
  nr <- length(buts)
  lignes <- c(sprintf("Optimisation de %d réponse(s) par la fonction de désirabilité de Derringer et Suich.", nr),
              sprintf("Désirabilité composite D = %s.", fr(op$ev$Dg, 4)))
  if (op$repli) avert <- c(avert, "Aucune solution ne satisfait simultanément toutes les réponses : la solution proposée est un compromis (moyenne des désirabilités).")
  lignes <- c(lignes, "", "Réglages optimaux :")
  valeurs <- lapply(seq_len(ctx$k), function(j) { v <- reel[[ctx$noms[j]]]; if (is.character(v)) v else signif(v, 6) })
  tab_f <- list(Facteur = ctx$noms, Valeur = valeurs)
  nom_but <- c(max = "Maximiser", min = "Minimiser", cible = "Atteindre la cible")
  tab_r <- list(Réponse = vapply(buts, function(b) b$nom, ""), But = nom_but[vapply(buts, function(b) b$but, "")],
                Prévision = as.numeric(op$ev$Y[1, ]), Désirabilité = as.numeric(op$ev$D[1, ]))
  pr <- lapply(seq_len(nr), function(r) predire(modeles[[r]]$mod, op$Xo))
  tab_r$IC_inf <- vapply(pr, function(p) p$ic_inf, 0); tab_r$IC_sup <- vapply(pr, function(p) p$ic_sup, 0)
  limites <- vapply(buts, function(b) switch(b$but, max = sprintf("inf. %s → cible %s", fr(b$inf), fr(b$cible)), min = sprintf("cible %s → sup. %s", fr(b$cible), fr(b$sup)),
                                              sprintf("%s ← %s → %s", fr(b$inf), fr(b$cible), fr(b$sup))), "")
  tab_r$Limites <- limites
  figs <- fig_optimisation(ctx, modeles, buts, op, "01_optimisation.png")
  optimum <- as.list(setNames(lapply(seq_len(ctx$k), function(j) { v <- reel[[ctx$noms[j]]]; if (is.character(v)) v else signif(v, 8) }), ctx$noms))
  list(titre = "Optimisation des réponses", lignes = lignes, tables = list(Réglages_optimaux = tab_f, Réponses_prévues = tab_r),
       figures = figs, avertissements = avert, conclusion = sprintf("Optimum : D = %s.", fr(op$ev$Dg, 3)), optimum = optimum, tables_libres = list())
}

# ---------------------------------------------------------------------------
#  Prévision
# ---------------------------------------------------------------------------
prevoir_doe <- function(ctx, par) {
  B <- preparer_modele(ctx, par$reponse, par$modele, isTRUE(par$reduire), as.numeric(par$alpha_sortie %||% 0.05))
  pts <- par$points
  if (!length(pts)) stop("Aucun point à prévoir.", call. = FALSE)
  d <- list()
  for (j in seq_len(ctx$k)) {
    nm <- ctx$noms[j]; f <- ctx$fac[[j]]
    v <- vapply(pts, function(p) { x <- p[[nm]]; if (is.null(x)) NA_character_ else as.character(x) }, "")
    if (any(is.na(v))) stop(sprintf("Valeur manquante pour « %s ».", nm), call. = FALSE)
    d[[nm]] <- if (identical(f$type, "cat")) v else as.numeric(v)
  }
  d <- as.data.frame(d, check.names = FALSE, stringsAsFactors = FALSE)
  if (ctx$famille == "melange") {
    s <- rowSums(as.matrix(d))
    if (any(abs(s - 1) > 1e-6)) stop(sprintf("Les proportions doivent totaliser 1 (somme saisie : %s).", paste(fr(s, 6), collapse = " ; ")), call. = FALSE)
    for (j in seq_len(ctx$k)) if (any(d[[j]] < ctx$L[j] - 1e-9 | d[[j]] > ctx$U[j] + 1e-9)) stop(sprintf("« %s » est hors de ses bornes [%s ; %s].", ctx$noms[j], fr(ctx$L[j]), fr(ctx$U[j])), call. = FALSE)
  }
  Xn <- coder_matrice(ctx, d)
  if (any(is.na(Xn))) stop("Valeurs de facteurs invalides.", call. = FALSE)
  pr <- predire(B$mod, Xn, 1 - ctx$alpha)
  rng <- apply(B$X[B$mod$ok, , drop = FALSE], 2, range)
  hors <- vapply(seq_len(nrow(Xn)), function(i) any(Xn[i, ] < rng[1, ] - 1e-9 | Xn[i, ] > rng[2, ] + 1e-9), TRUE)
  niv <- 100 * (1 - ctx$alpha)
  tab <- list(Point = seq_len(nrow(Xn)))
  for (j in seq_len(ctx$k)) tab[[ctx$noms[j]]] <- d[[ctx$noms[j]]]
  tab$Prévision <- pr$fit; tab$Erreur_type <- pr$se
  tab[[sprintf("IC_%g%%_inf", niv)]] <- pr$ic_inf; tab[[sprintf("IC_%g%%_sup", niv)]] <- pr$ic_sup
  tab[[sprintf("IP_%g%%_inf", niv)]] <- pr$ip_inf; tab[[sprintf("IP_%g%%_sup", niv)]] <- pr$ip_sup
  avert <- character(0)
  if (any(hors)) avert <- c(avert, sprintf("Point(s) %s hors du domaine expérimenté : extrapolation peu fiable.", paste(which(hors), collapse = ", ")))
  if (B$mod$df <= 0) avert <- c(avert, "Aucun degré de liberté pour l'erreur : intervalles de confiance et de prévision non calculables.")
  lignes <- c(sprintf("Réponse : %s.", par$reponse),
              sprintf("Modèle : %d terme(s)%s.", length(B$mod$termes), if (isTRUE(par$reduire)) " (modèle réduit)" else ""),
              sprintf("IC : intervalle de confiance de la valeur moyenne à %g %% ; IP : intervalle de prévision d'une valeur individuelle à %g %%.", niv, niv),
              "", "Équation (unités codées) :", paste0("  ", equation_texte(par$reponse, B$mod$coef, vapply(B$mod$termes, function(t) t$label, ""), B$mod$intercept)))
  list(titre = sprintf("Prévision de « %s »", par$reponse), lignes = lignes, tables = list(Prévisions = tab), figures = character(0), avertissements = avert,
       conclusion = sprintf("Prévision : %s.", paste(fr(pr$fit, 5), collapse = " ; ")), tables_libres = list("Prévisions"))
}

# ---------------------------------------------------------------------------
#  Surfaces (contours, 3D, ternaire, trace)
# ---------------------------------------------------------------------------
surfaces_doe <- function(ctx, par) {
  B <- preparer_modele(ctx, par$reponse, par$modele, isTRUE(par$reduire), as.numeric(par$alpha_sortie %||% 0.05))
  fs <- list(); figures <- character(0); lignes <- character(0); avert <- character(0)
  fixes <- par$fixes %||% list()
  type <- par$type %||% "contour"
  if (ctx$famille == "melange") {
    comp <- as.integer(unlist(par$composants_axes))
    if (length(comp) < 3 && ctx$k >= 3) comp <- 1:3
    if (ctx$k == 2) stop("Avec deux composants, utiliser le tracé des effets (un seul degré de liberté).", call. = FALSE)
    ref <- rep(NA_real_, ctx$k)
    for (j in seq_len(ctx$k)) { v <- fixes[[ctx$noms[j]]]; if (!is.null(v) && nzchar(as.character(v))) ref[j] <- as.numeric(v) }
    figures <- c(figures, fig_ternaire(ctx, B, comp[1:3], ref, "01_ternaire.png", par$reponse))
    figures <- c(figures, fig_trace_melange(ctx, B, "02_trace.png"))
    lignes <- c("Diagramme ternaire du modèle ajusté : les lignes de niveau sont les valeurs prévues de la réponse.",
                if (ctx$k > 3) "Les composants qui ne figurent pas sur le triangle sont maintenus à leur proportion de référence (centre du domaine par défaut) ; le triangle représente la part restante." else NULL)
  } else {
    xs <- par$axe_x; ys <- par$axe_y
    ix <- match(xs, ctx$noms); iy <- match(ys, ctx$noms)
    if (is.na(ix) || is.na(iy) || ix == iy) stop("Choisir deux facteurs différents pour les axes.", call. = FALSE)
    if (any(vapply(ctx$fac[c(ix, iy)], function(f) identical(f$type, "cat"), TRUE))) stop("Les axes doivent être des facteurs quantitatifs.", call. = FALSE)
    base <- rep(0, ctx$k)
    for (j in seq_len(ctx$k)) {
      v <- fixes[[ctx$noms[j]]]
      if (!is.null(v) && nzchar(as.character(v)) && j != ix && j != iy) {
        if (identical(ctx$fac[[j]]$type, "cat")) base[j] <- if (as.character(v) == unlist(ctx$fac[[j]]$niveaux)[1]) -1 else 1
        else { f <- ctx$fac[[j]]; base[j] <- (as.numeric(v) - (f$bas + f$haut) / 2) / ((f$haut - f$bas) / 2) }
      }
    }
    if (type %in% c("contour", "les_deux")) figures <- c(figures, fig_contour(ctx, B, ix, iy, base, "01_contour.png", par$reponse))
    if (type %in% c("surface", "les_deux")) figures <- c(figures, fig_surface3d(ctx, B, ix, iy, base, "02_surface.png", par$reponse))
    if (ctx$k > 2) {
      autres <- setdiff(seq_len(ctx$k), c(ix, iy))
      lignes <- c(lignes, sprintf("Facteurs maintenus constants : %s.", paste(sprintf("%s = %s", ctx$noms[autres], fr(unlist(decoder_reel(ctx, matrix(base, 1, dimnames = list(NULL, ctx$noms)))[autres]), 5)), collapse = " ; ")))
    }
  }
  list(titre = sprintf("Représentation du modèle — %s", par$reponse), lignes = c(sprintf("Réponse : %s.", par$reponse), lignes), tables = list(), figures = figures, avertissements = avert, conclusion = "", tables_libres = list())
}
