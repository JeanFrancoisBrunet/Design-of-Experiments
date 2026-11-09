# -*- coding: utf-8 -*-
# =============================================================================
#  graphiques_doe.R — Graphiques des plans d'expériences (R de base, fichiers PNG)
# =============================================================================
COL_BLEU <- "#1f4e8c"; COL_ROUGE <- "#c0392b"; COL_GRIS <- "#808080"; COL_VERT <- "#2e7d32"; COL_ORANGE <- "#e08a00"

ouvrir_png <- function(dossier, nom, ech = 1, larg = 900, haut = 600) {
  f <- file.path(dossier, nom)
  args <- list(filename = f, width = round(larg * ech), height = round(haut * ech), res = round(96 * ech))
  ok <- tryCatch({ do.call(grDevices::png, c(args, list(type = "cairo"))); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) do.call(grDevices::png, args)
  invisible(nom)
}
fermer_png <- function() invisible(grDevices::dev.off())

palette_surface <- function(n = 40) grDevices::hcl.colors(n, "YlGnBu", rev = TRUE)

valeur_axe <- function(ctx, j, code) {
  f <- ctx$fac[[j]]
  if (ctx$famille == "melange") return(if (ctx$pseudo) ctx$L[j] + code * (1 - ctx$sL) else code)
  if (identical(f$type, "cat")) return(ifelse(code < 0, unlist(f$niveaux)[1], unlist(f$niveaux)[2]))
  (f$bas + f$haut) / 2 + code * (f$haut - f$bas) / 2
}
etiquette_niveau <- function(ctx, j, code) { v <- valeur_axe(ctx, j, code); if (is.character(v)) v else fr(v, 4) }

grille_panneaux <- function(n, max_col = 4) { nc <- min(n, max_col); c(ceiling(n / nc), nc) }

# ---------------------------------------------------------------------------
fig_pareto <- function(dossier, nom, labels, valeurs, seuil, titre, xlab, ech, ligne_lib = NULL) {
  o <- order(valeurs, decreasing = TRUE)
  o <- o[seq_len(min(length(o), 25))]
  v <- valeurs[o]; lab <- labels[o]
  ouvrir_png(dossier, nom, ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mai = c(0.9, 0.3, 0.8, 0.3))
  wl <- max(graphics::strwidth(lab, units = "inches", cex = 0.9)) + 0.25
  graphics::par(mai = c(0.9, min(wl, 4), 0.8, 0.4))
  xmax <- max(c(v, seuil), na.rm = TRUE) * 1.12
  cols <- ifelse(v > seuil, COL_ROUGE, COL_BLEU)
  graphics::barplot(rev(v), horiz = TRUE, names.arg = rev(lab), las = 1, col = rev(cols), border = NA, xlim = c(0, xmax), cex.names = 0.9, xlab = xlab)
  graphics::abline(v = seuil, col = COL_ROUGE, lty = 2, lwd = 2)
  graphics::title(main = titre, cex.main = 1.1)
  if (!is.null(ligne_lib)) graphics::mtext(ligne_lib, side = 3, line = 0, cex = 0.85, col = COL_ROUGE)
  graphics::box(col = COL_GRIS)
  nom
}

fig_daniel <- function(dossier, nom, labels, effets, signif, ech) {
  m <- length(effets)
  o <- order(effets); e <- effets[o]; lab <- labels[o]; sg <- signif[o]
  q <- stats::qnorm(stats::ppoints(m))
  ouvrir_png(dossier, nom, ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mar = c(4.5, 4.5, 3.5, 1.5))
  graphics::plot(e, q, pch = 19, col = ifelse(sg, COL_ROUGE, COL_BLEU), cex = 1.2, xlab = "Effet", ylab = "Score normal", main = "Diagramme de Daniel (probabilité normale des effets)",
                 xlim = range(e) + c(-0.15, 0.15) * diff(range(e)))
  nsg <- which(!sg)
  if (length(nsg) >= 2 && stats::sd(e[nsg]) > 0) { fit <- stats::lm(q[nsg] ~ e[nsg]); graphics::abline(a = stats::coef(fit)[1], b = stats::coef(fit)[2], col = COL_GRIS, lty = 2, lwd = 2) }
  graphics::abline(v = 0, col = "#cccccc")
  vis <- if (m <= 30) seq_len(m) else which(sg)
  graphics::text(e[vis], q[vis], lab[vis], pos = ifelse(e[vis] >= 0, 4, 2), cex = 0.85, col = ifelse(sg[vis], COL_ROUGE, "#333333"), xpd = TRUE)
  graphics::legend("topleft", legend = c("Effet significatif", "Effet non significatif"), pch = 19, col = c(COL_ROUGE, COL_BLEU), bty = "n", cex = 0.9)
  graphics::box(col = COL_GRIS)
  nom
}

# Effets principaux d'après les données (moyennes aux niveaux)
fig_effets_principaux <- function(ctx, y, mod, B, nom) {
  X <- B$X; ok <- !is.na(y); k <- ctx$k
  moy <- mean(y[ok])
  dims <- grille_panneaux(k)
  ouvrir_png(ctx$dossier, nom, ctx$ech, 900, if (dims[1] > 2) 700 else 600)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mfrow = dims, mar = c(4, 3.5, 2.5, 0.8), oma = c(0, 1.5, 2, 0))
  mi <- vapply(seq_len(k), function(j) { r <- c(mean(y[ok & X[, j] < -0.5]), mean(y[ok & X[, j] > 0.5])); r }, numeric(2))
  yc <- if (any(ok & rowSums(abs(X)) < 1e-9)) mean(y[ok & rowSums(abs(X)) < 1e-9]) else NA
  yl <- range(c(mi, yc), na.rm = TRUE); yl <- yl + c(-0.08, 0.08) * max(diff(yl), 1e-9)
  for (j in seq_len(k)) {
    graphics::plot(c(-1, 1), mi[, j], type = "b", pch = 19, lwd = 2, col = COL_BLEU, xlim = c(-1.3, 1.3), ylim = yl, xaxt = "n", xlab = "", ylab = if ((j - 1) %% dims[2] == 0) "Moyenne" else "", main = ctx$noms[j], cex.main = 1)
    graphics::axis(1, at = c(-1, 1), labels = c(etiquette_niveau(ctx, j, -1), etiquette_niveau(ctx, j, 1)))
    if (!is.na(yc) && !identical(ctx$fac[[j]]$type, "cat")) { graphics::points(0, yc, pch = 18, cex = 1.6, col = COL_ORANGE); graphics::axis(1, at = 0, labels = etiquette_niveau(ctx, j, 0)) }
    graphics::abline(h = moy, lty = 2, col = COL_GRIS)
    graphics::box(col = COL_GRIS)
  }
  graphics::mtext(sprintf("Effets principaux sur %s (moyennes des essais)", B$reponse), outer = TRUE, font = 2, cex = 1)
  if (!is.na(yc)) graphics::mtext("◆ moyenne des points centraux", outer = TRUE, side = 1, line = -1, cex = 0.7, col = COL_ORANGE, adj = 1)
  nom
}

fig_interactions <- function(ctx, y, mod, nom, signif_labels = character(0)) {
  X <- coder_matrice(ctx, ctx$plan); ok <- !is.na(y)
  paires <- list(); cf <- numeric(0)
  off <- if (mod$intercept) 1L else 0L
  for (i in seq_along(mod$termes)) {
    t <- mod$termes[[i]]
    if (is.null(t$special) && is.null(t$fn) && sum(t$expo) == 2 && sum(t$expo > 0) == 2) { paires[[length(paires) + 1]] <- which(t$expo > 0); cf <- c(cf, abs(mod$coef[off + i])) }
  }
  if (!length(paires)) return(NULL)
  o <- order(cf, decreasing = TRUE); paires <- paires[o[seq_len(min(12, length(o)))]]
  dims <- grille_panneaux(length(paires))
  ouvrir_png(ctx$dossier, nom, ctx$ech, 900, if (dims[1] > 2) 700 else 600)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mfrow = dims, mar = c(4, 3.5, 2.5, 0.8), oma = c(0, 1.5, 2, 0))
  cell <- function(a, b, sa, sb) { s <- ok & X[, a] * sa > 0.5 & X[, b] * sb > 0.5; if (any(s)) mean(y[s]) else NA }
  allv <- unlist(lapply(paires, function(p) c(cell(p[1], p[2], -1, -1), cell(p[1], p[2], 1, -1), cell(p[1], p[2], -1, 1), cell(p[1], p[2], 1, 1))))
  yl <- range(allv, na.rm = TRUE); yl <- yl + c(-0.08, 0.08) * max(diff(yl), 1e-9)
  for (p in paires) {
    a <- p[1]; b <- p[2]
    graphics::plot(c(-1, 1), c(cell(a, b, -1, -1), cell(a, b, 1, -1)), type = "b", pch = 19, lwd = 2, col = COL_BLEU, xlim = c(-1.3, 1.3), ylim = yl, xaxt = "n",
                   xlab = ctx$noms[a], ylab = "", main = paste(ctx$noms[a], "×", ctx$noms[b]), cex.main = 0.95)
    graphics::lines(c(-1, 1), c(cell(a, b, -1, 1), cell(a, b, 1, 1)), type = "b", pch = 17, lwd = 2, col = COL_ROUGE)
    graphics::axis(1, at = c(-1, 1), labels = c(etiquette_niveau(ctx, a, -1), etiquette_niveau(ctx, a, 1)))
    graphics::legend("top", legend = paste(ctx$noms[b], "=", c(etiquette_niveau(ctx, b, -1), etiquette_niveau(ctx, b, 1))), col = c(COL_BLEU, COL_ROUGE), pch = c(19, 17), lty = 1, bty = "n", cex = 0.75, horiz = FALSE)
    graphics::box(col = COL_GRIS)
  }
  graphics::mtext(sprintf("Interactions sur %s (moyennes des essais)", y_nom(ctx, y)), outer = TRUE, font = 2, cex = 1)
  nom
}
y_nom <- function(ctx, y) { for (r in ctx$reps) if (identical(ctx$plan[[r]], y)) return(r); "la réponse" }

fig_cube <- function(ctx, y, nom) {
  X <- coder_matrice(ctx, ctx$plan); ok <- !is.na(y)
  ouvrir_png(ctx$dossier, nom, ctx$ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mar = c(1, 1, 3, 1))
  graphics::plot(NA, xlim = c(-0.25, 1.7), ylim = c(-0.25, 1.55), axes = FALSE, xlab = "", ylab = "", asp = 1, main = sprintf("Diagramme du cube — moyennes de %s", y_nom(ctx, y)))
  dx <- 0.45; dy <- 0.35
  pos <- function(a, b, c) c(a + c * dx, b + c * dy)
  for (a in 0:1) for (b in 0:1) for (c in 0:1) {
    s <- ok & X[, 1] * (2 * a - 1) > 0.5 & X[, 2] * (2 * b - 1) > 0.5 & X[, 3] * (2 * c - 1) > 0.5
    m <- if (any(s)) mean(y[s]) else NA
    p <- pos(a, b, c)
    graphics::points(p[1], p[2], pch = 19, cex = 1.6, col = COL_BLEU)
    graphics::text(p[1] + (if (a == 1) 0.04 else -0.04), p[2] + (if (b == 1) 0.05 else -0.05), if (is.na(m)) "—" else fr(m, 4), cex = 1, font = 2, adj = c(if (a == 1) 0 else 1, 0.5))
  }
  for (a in 0:1) for (b in 0:1) graphics::segments(pos(a, b, 0)[1], pos(a, b, 0)[2], pos(a, b, 1)[1], pos(a, b, 1)[2], col = COL_GRIS)
  for (b in 0:1) for (c in 0:1) graphics::segments(pos(0, b, c)[1], pos(0, b, c)[2], pos(1, b, c)[1], pos(1, b, c)[2], col = COL_GRIS)
  for (a in 0:1) for (c in 0:1) graphics::segments(pos(a, 0, c)[1], pos(a, 0, c)[2], pos(a, 1, c)[1], pos(a, 1, c)[2], col = COL_GRIS)
  graphics::text(0.5, -0.18, sprintf("%s : %s → %s", ctx$noms[1], etiquette_niveau(ctx, 1, -1), etiquette_niveau(ctx, 1, 1)), cex = 0.95)
  graphics::text(-0.2, 0.5, sprintf("%s\n%s → %s", ctx$noms[2], etiquette_niveau(ctx, 2, -1), etiquette_niveau(ctx, 2, 1)), cex = 0.95, srt = 90)
  graphics::text(1.5, 1.45, sprintf("%s : %s → %s (profondeur)", ctx$noms[3], etiquette_niveau(ctx, 3, -1), etiquette_niveau(ctx, 3, 1)), cex = 0.95, adj = 1)
  nom
}

fig_residus <- function(dossier, nom, mod, ordre, reponse, ech) {
  s <- sqrt(mod$s2); rs <- mod$res / (s * sqrt(pmax(1 - mod$h, 1e-12)))
  ouvrir_png(dossier, nom, ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  qq <- stats::qqnorm(rs, main = "Probabilité normale des résidus", xlab = "Score normal", ylab = "Résidu standardisé", pch = 19, col = COL_BLEU); stats::qqline(rs, col = COL_ROUGE, lty = 2)
  graphics::plot(mod$fitted, rs, pch = 19, col = COL_BLEU, xlab = "Valeur ajustée", ylab = "Résidu standardisé", main = "Résidus selon les valeurs ajustées"); graphics::abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c(COL_ROUGE, COL_GRIS, COL_ROUGE))
  graphics::hist(rs, col = "#9db7d9", border = "white", main = "Histogramme des résidus", xlab = "Résidu standardisé", ylab = "Fréquence")
  graphics::plot(ordre, rs, type = "b", pch = 19, col = COL_BLEU, xlab = "Ordre d'exécution", ylab = "Résidu standardisé", main = "Résidus selon l'ordre d'exécution"); graphics::abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c(COL_ROUGE, COL_GRIS, COL_ROUGE))
  graphics::mtext(sprintf("Analyse des résidus — %s", reponse), outer = TRUE, font = 2)
  nom
}

# Effets principaux d'après le modèle (surface de réponse)
fig_effets_modele <- function(ctx, B, nom) {
  k <- ctx$k; rng <- apply(B$X, 2, range)
  dims <- grille_panneaux(k)
  ouvrir_png(ctx$dossier, nom, ctx$ech, 900, if (dims[1] > 2) 700 else 600)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mfrow = dims, mar = c(4, 3.5, 2.5, 0.8), oma = c(0, 1.5, 2, 0))
  grids <- lapply(seq_len(k), function(j) {
    g <- seq(rng[1, j], rng[2, j], length.out = 50); Xg <- matrix(0, 50, k, dimnames = list(NULL, ctx$noms)); Xg[, j] <- g; list(g = g, y = pred_simple(B$mod, Xg)) })
  yl <- range(unlist(lapply(grids, `[[`, "y"))); yl <- yl + c(-0.05, 0.05) * max(diff(yl), 1e-9)
  for (j in seq_len(k)) {
    graphics::plot(valeur_axe(ctx, j, grids[[j]]$g), grids[[j]]$y, type = "l", lwd = 2, col = COL_BLEU, ylim = yl, xlab = ctx$noms[j], ylab = if ((j - 1) %% dims[2] == 0) "Réponse prévue" else "", main = ctx$noms[j])
    graphics::box(col = COL_GRIS)
  }
  graphics::mtext(sprintf("Effet de chaque facteur sur %s (autres facteurs au centre)", B$reponse), outer = TRUE, font = 2, cex = 1)
  nom
}

# ---------------------------------------------------------------------------
#  Contour, surface 3D
# ---------------------------------------------------------------------------
barre_couleurs <- function(zlim, pal) {
  graphics::par(mar = c(4, 0.5, 3, 3.2))
  graphics::image(1, seq(zlim[1], zlim[2], length.out = length(pal)), matrix(seq_along(pal), 1), col = pal, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, las = 1); graphics::box(col = COL_GRIS)
}

fig_contour <- function(ctx, B, ix, iy, base, nom, reponse) {
  rng <- apply(B$X, 2, range)
  gx <- seq(rng[1, ix], rng[2, ix], length.out = 70); gy <- seq(rng[1, iy], rng[2, iy], length.out = 70)
  G <- expand.grid(gx, gy)
  Xg <- matrix(rep(base, each = nrow(G)), nrow(G), ctx$k, dimnames = list(NULL, ctx$noms))
  Xg[, ix] <- G[, 1]; Xg[, iy] <- G[, 2]
  z <- matrix(pred_simple(B$mod, Xg), 70, 70)
  rx <- valeur_axe(ctx, ix, gx); ry <- valeur_axe(ctx, iy, gy)
  pal <- palette_surface(); zl <- range(z)
  ouvrir_png(ctx$dossier, nom, ctx$ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::layout(matrix(1:2, 1), widths = c(6, 1))
  graphics::par(mar = c(4.5, 4.5, 3.5, 0.5))
  graphics::image(rx, ry, z, col = pal, zlim = zl, xlab = ctx$noms[ix], ylab = ctx$noms[iy], main = sprintf("Lignes de niveau de %s", reponse))
  graphics::contour(rx, ry, z, add = TRUE, labcex = 0.9, col = "#222222", nlevels = 10)
  Xp <- B$X[B$mod$ok, , drop = FALSE]
  graphics::points(valeur_axe(ctx, ix, Xp[, ix]), valeur_axe(ctx, iy, Xp[, iy]), pch = 21, bg = "white", col = "#222222", cex = 1.1)
  graphics::box(col = COL_GRIS)
  barre_couleurs(zl, pal)
  nom
}

fig_surface3d <- function(ctx, B, ix, iy, base, nom, reponse) {
  rng <- apply(B$X, 2, range)
  n <- 35
  gx <- seq(rng[1, ix], rng[2, ix], length.out = n); gy <- seq(rng[1, iy], rng[2, iy], length.out = n)
  G <- expand.grid(gx, gy)
  Xg <- matrix(rep(base, each = nrow(G)), nrow(G), ctx$k, dimnames = list(NULL, ctx$noms))
  Xg[, ix] <- G[, 1]; Xg[, iy] <- G[, 2]
  z <- matrix(pred_simple(B$mod, Xg), n, n)
  rx <- valeur_axe(ctx, ix, gx); ry <- valeur_axe(ctx, iy, gy)
  pal <- palette_surface(60)
  zf <- (z[-1, -1] + z[-n, -1] + z[-1, -n] + z[-n, -n]) / 4
  idx <- pmax(1, pmin(60, 1 + floor(59 * (zf - min(z)) / max(diff(range(z)), 1e-12))))
  ouvrir_png(ctx$dossier, paste0(nom), ctx$ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mar = c(1, 1, 3, 1))
  graphics::persp(rx, ry, z, theta = 35, phi = 25, expand = 0.7, col = pal[idx], border = "#ffffff66", ticktype = "detailed", nticks = 5, shade = 0.25,
                  xlab = ctx$noms[ix], ylab = ctx$noms[iy], zlab = reponse, main = sprintf("Surface de réponse — %s", reponse), cex.axis = 0.75)
  nom
}

# ---------------------------------------------------------------------------
#  Mélange : trace de Cox et diagramme ternaire
# ---------------------------------------------------------------------------
centre_melange <- function(ctx) {
  V <- sommets_extremes(ctx$L, ctx$U)
  if (is.null(V)) rep(1 / ctx$k, ctx$k) else colMeans(V)
}

pred_melange_reel <- function(ctx, B, M) {
  d <- as.data.frame(M, check.names = FALSE); names(d) <- ctx$noms
  pred_simple(B$mod, coder_matrice(ctx, d))
}

fig_trace_melange <- function(ctx, B, nom) {
  q <- ctx$k; c0 <- centre_melange(ctx)
  cols <- grDevices::hcl.colors(max(q, 3), "Dark 3")
  tr <- lapply(seq_len(q), function(i) {
    xi <- seq(ctx$L[i], ctx$U[i], length.out = 60)
    M <- t(vapply(xi, function(v) { x <- c0 * (1 - v) / (1 - c0[i]); x[i] <- v; x }, numeric(q)))
    ok <- apply(M, 1, function(r) all(r >= ctx$L - 1e-9 & r <= ctx$U + 1e-9))
    list(x = xi[ok] - c0[i], y = pred_melange_reel(ctx, B, M[ok, , drop = FALSE]))
  })
  yr <- range(unlist(lapply(tr, `[[`, "y")), na.rm = TRUE); xr <- range(unlist(lapply(tr, `[[`, "x")))
  ouvrir_png(ctx$dossier, nom, ctx$ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mar = c(4.5, 4.5, 3.5, 1.5))
  graphics::plot(NA, xlim = xr, ylim = yr, xlab = "Écart de proportion par rapport au mélange de référence", ylab = "Réponse prévue", main = sprintf("Trace de la réponse (direction de Cox) — %s", B$reponse))
  graphics::abline(v = 0, lty = 2, col = COL_GRIS)
  for (i in seq_len(q)) graphics::lines(tr[[i]]$x, tr[[i]]$y, col = cols[i], lwd = 2.5)
  graphics::legend("topleft", legend = ctx$noms, col = cols[seq_len(q)], lwd = 2.5, bty = "n", cex = 0.9)
  graphics::mtext(sprintf("Mélange de référence : %s", paste(sprintf("%s = %s", ctx$noms, fr(c0, 3)), collapse = " ; ")), side = 1, line = 3.3, cex = 0.7, col = COL_GRIS)
  graphics::box(col = COL_GRIS)
  nom
}

fig_ternaire <- function(ctx, B, comp, ref, nom, reponse) {
  q <- ctx$k; a <- comp[1]; b <- comp[2]; c <- comp[3]
  c0 <- centre_melange(ctx)
  ref[is.na(ref)] <- c0[is.na(ref)]
  autres <- setdiff(seq_len(q), comp)
  s <- 1 - sum(ref[autres])
  if (s <= 1e-6) stop("La somme des autres composants laisse une part nulle aux trois composants du triangle.", call. = FALSE)
  H <- sqrt(3) / 2
  nx <- 150; ny <- 130
  gx <- seq(0, 1, length.out = nx); gy <- seq(0, H, length.out = ny)
  G <- expand.grid(gx, gy)
  ua <- G[, 2] / H; uc <- G[, 1] - 0.5 * ua; ub <- 1 - ua - uc
  dedans <- ua >= 0 & ub >= -1e-12 & uc >= -1e-12
  M <- matrix(0, nrow(G), q)
  M[, a] <- ua * s; M[, b] <- ub * s; M[, c] <- uc * s
  for (j in autres) M[, j] <- ref[j]
  faisable <- dedans & apply(M, 1, function(r) all(r >= ctx$L - 1e-9 & r <= ctx$U + 1e-9))
  z <- rep(NA_real_, nrow(G))
  if (any(faisable)) z[faisable] <- pred_melange_reel(ctx, B, M[faisable, , drop = FALSE])
  if (!any(faisable)) stop("Aucun mélange faisable sur ce triangle avec ces proportions de référence.", call. = FALSE)
  zm <- matrix(z, nx, ny)
  pal <- palette_surface(); zl <- range(z, na.rm = TRUE)
  if (diff(zl) < 1e-12) zl <- zl + c(-0.5, 0.5)
  ouvrir_png(ctx$dossier, nom, ctx$ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::layout(matrix(1:2, 1), widths = c(6, 1))
  graphics::par(mar = c(3, 2, 3.5, 2))
  graphics::plot(NA, xlim = c(-0.12, 1.12), ylim = c(-0.1, H + 0.09), asp = 1, axes = FALSE, xlab = "", ylab = "", main = sprintf("Diagramme ternaire — %s", reponse))
  graphics::image(gx, gy, zm, col = pal, zlim = zl, add = TRUE)
  if (sum(!is.na(zm)) > 20 && diff(range(zm, na.rm = TRUE)) > 1e-9) try(graphics::contour(gx, gy, zm, add = TRUE, labcex = 0.8, col = "#222222", nlevels = 8), silent = TRUE)
  graphics::polygon(c(0, 1, 0.5), c(0, 0, H), border = "#222222", lwd = 2)
  # graduations : proportions de chaque composant (a : lignes horizontales, b et c : lignes obliques)
  for (u in c(0.2, 0.4, 0.6, 0.8)) {
    graphics::segments(0.5 * u, H * u, 1 - 0.5 * u, H * u, col = "#ffffffaa", lty = 3)
    graphics::segments(u, 0, 0.5 + 0.5 * u, H * (1 - u), col = "#ffffffaa", lty = 3)
    graphics::segments(1 - u, 0, 0.5 * (1 - u), H * (1 - u), col = "#ffffffaa", lty = 3)
    graphics::text(1 - 0.5 * u + 0.012, H * u, fr(u * s, 3), adj = c(0, 0.5), cex = 0.7, xpd = TRUE, col = "#333333")
    graphics::text(u, -0.022, fr(u * s, 3), adj = c(0.5, 1), cex = 0.7, xpd = TRUE, col = "#333333")
    graphics::text(0.5 * (1 - u) - 0.012, H * (1 - u), fr(u * s, 3), adj = c(1, 0.5), cex = 0.7, xpd = TRUE, col = "#333333")
  }
  graphics::text(0.5, H + 0.06, ctx$noms[a], font = 2, cex = 1); graphics::text(0, -0.07, ctx$noms[b], font = 2, cex = 1); graphics::text(1, -0.07, ctx$noms[c], font = 2, cex = 1)
  graphics::mtext(sprintf("Sommets : proportion maximale de chaque composant dans le triangle = %s%s", fr(s, 3),
                          if (length(autres)) sprintf(" ; autres composants fixés : %s", paste(sprintf("%s = %s", ctx$noms[autres], fr(ref[autres], 3)), collapse = ", ")) else ""), side = 1, line = 1.5, cex = 0.75, col = COL_GRIS)
  if (q == 3) {
    Xp <- as.matrix(ctx$plan[, ctx$noms]); mode(Xp) <- "numeric"
    ua <- Xp[, a] / s; uc <- Xp[, c] / s
    graphics::points(0.5 * ua + uc, H * ua, pch = 21, bg = "white", col = "#222222", cex = 1)
  }
  barre_couleurs(zl, pal)
  nom
}

# ---------------------------------------------------------------------------
#  Diagramme d'optimisation
# ---------------------------------------------------------------------------
fig_optimisation <- function(ctx, modeles, buts, op, nom) {
  k <- ctx$k; nr <- length(modeles)
  xo <- as.numeric(op$Xo)
  larg <- max(700, min(1500, 170 * k + 150)); haut <- max(520, min(1100, 150 * (nr + 1) + 120))
  ouvrir_png(ctx$dossier, nom, ctx$ech, larg, haut)
  oldp <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(oldp); fermer_png() })
  graphics::par(mfrow = c(nr + 1, k), mar = c(3, 3.2, 1.6, 0.6), oma = c(0, 4.5, 2.5, 0))
  rng <- if (ctx$famille == "melange") rbind(ctx$L, ctx$U) else apply(op$X0, 2, range)
  for (r in seq_len(nr + 1)) for (j in seq_len(k)) {
    cat_j <- identical(ctx$fac[[j]]$type, "cat")
    if (ctx$famille == "melange") {
      c0 <- xo_reel <- as.numeric(unlist(decoder_reel(ctx, op$Xo)))
      xi <- sort(unique(c(seq(ctx$L[j], ctx$U[j], length.out = 40), xo_reel[j])))
      Mr <- t(vapply(xi, function(v) { x <- xo_reel * (1 - v) / max(1 - xo_reel[j], 1e-9); x[j] <- v; x }, numeric(k)))
      okm <- apply(Mr, 1, function(rw) all(rw >= ctx$L - 1e-7 & rw <= ctx$U + 1e-7))
      Mr <- Mr[okm, , drop = FALSE]; xi <- xi[okm]
      d <- as.data.frame(Mr, check.names = FALSE); names(d) <- ctx$noms
      Xc <- coder_matrice(ctx, d); xplot <- xi; xopt <- xo_reel[j]
    } else {
      g <- if (cat_j) c(-1, 1) else seq(rng[1, j], rng[2, j], length.out = 40)
      Xc <- matrix(rep(xo, each = length(g)), length(g), k, dimnames = list(NULL, ctx$noms)); Xc[, j] <- g
      xplot <- if (cat_j) c(-1, 1) else valeur_axe(ctx, j, g); xopt <- if (cat_j) xo[j] else valeur_axe(ctx, j, xo[j])
    }
    if (r <= nr) { yv <- pred_simple(modeles[[r]]$mod, Xc); yopt <- op$ev$Y[1, r] } else {
      ev <- op$evalue(Xc); yv <- ev$Dg; yopt <- op$ev$Dg }
    xl <- if (length(xplot) < 2 || diff(range(xplot)) < 1e-9) xopt + c(-0.02, 0.02) else range(xplot)
    graphics::plot(xplot, yv, xlim = xl, type = if (cat_j) "b" else "l", pch = 19, lwd = 2, col = if (r <= nr) COL_BLEU else COL_VERT, xaxt = if (cat_j) "n" else "s", xlab = "", ylab = "",
                   main = if (r == 1) ctx$noms[j] else "", cex.main = 1, ylim = if (r > nr) c(0, 1) else range(yv, yopt), cex.axis = 0.8)
    if (cat_j) graphics::axis(1, at = c(-1, 1), labels = unlist(ctx$fac[[j]]$niveaux), cex.axis = 0.8)
    graphics::abline(v = xopt, col = COL_ROUGE, lwd = 1.5); graphics::abline(h = yopt, col = COL_ROUGE, lty = 2)
    graphics::mtext(if (cat_j) as.character(valeur_axe(ctx, j, xo[j])) else fr(xopt, 4), side = 1, line = 1.8, cex = 0.7, col = COL_ROUGE)
    if (j == 1) graphics::mtext(if (r <= nr) sprintf("%s\n= %s", buts[[r]]$nom, fr(yopt, 4)) else sprintf("D composite\n= %s", fr(yopt, 3)), side = 2, line = 3, cex = 0.75, font = 2)
    graphics::box(col = COL_GRIS)
  }
  graphics::mtext("Diagramme d'optimisation (ligne rouge : réglage optimal)", outer = TRUE, font = 2)
  c(nom, fig_desirabilites(ctx, buts, op, "02_desirabilites.png"))
}

fig_desirabilites <- function(ctx, buts, op, nom) {
  nr <- length(buts)
  d <- c(as.numeric(op$ev$D[1, ]), op$ev$Dg)
  lab <- c(vapply(buts, function(b) b$nom, ""), "Composite D")
  ouvrir_png(ctx$dossier, nom, ctx$ech)
  oldp <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(oldp); fermer_png() })
  graphics::par(mar = c(4.5, max(6, 0.6 * max(nchar(lab))), 3.5, 2))
  b <- graphics::barplot(rev(d), horiz = TRUE, names.arg = rev(lab), las = 1, col = rev(c(rep(COL_BLEU, nr), COL_VERT)), border = NA, xlim = c(0, 1.1), xlab = "Désirabilité (0 = inacceptable, 1 = idéal)",
                         main = "Désirabilités individuelle et composite")
  graphics::text(rev(d) + 0.02, b, fr(rev(d), 3), adj = 0, cex = 0.9)
  graphics::abline(v = 1, lty = 2, col = COL_GRIS)
  nom
}
