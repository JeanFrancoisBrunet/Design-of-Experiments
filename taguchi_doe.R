# -*- coding: utf-8 -*-
# =============================================================================
#  taguchi_doe.R — Analyse des plans de Taguchi
#    rapport signal/bruit (RSB), tableaux des réponses, ANOVA avec regroupement,
#    niveaux optimaux, prévision et intervalle de confirmation.
# =============================================================================
TYPES_RSB <- c(plus_grand = "Plus grand est meilleur", plus_petit = "Plus petit est meilleur",
               nominal_1 = "Cible nominale (moyenne et dispersion)", nominal_2 = "Cible nominale (dispersion seule)")

rsb_essais <- function(Y, type) {
  moy <- rowMeans(Y, na.rm = TRUE); n <- rowSums(!is.na(Y))
  s2 <- apply(Y, 1, stats::var, na.rm = TRUE)
  if (type == "plus_grand") {
    if (any(Y <= 0, na.rm = TRUE)) stop("« Plus grand est meilleur » exige des réponses strictement positives.", call. = FALSE)
    return(list(moy = moy, s = sqrt(s2), rsb = -10 * log10(rowMeans(1 / Y^2, na.rm = TRUE))))
  }
  if (type == "plus_petit") return(list(moy = moy, s = sqrt(s2), rsb = -10 * log10(rowMeans(Y^2, na.rm = TRUE))))
  if (any(n < 2)) stop("Ce type de rapport signal/bruit exige au moins deux réponses par essai (plusieurs colonnes de réponses ou répétitions).", call. = FALSE)
  if (any(!is.finite(s2) | s2 <= 0)) stop("Un essai a une dispersion nulle : le rapport signal/bruit n'est pas défini.", call. = FALSE)
  if (type == "nominal_1") { if (any(moy == 0)) stop("Une moyenne est nulle : rapport signal/bruit non défini.", call. = FALSE); return(list(moy = moy, s = sqrt(s2), rsb = 10 * log10(moy^2 / s2))) }
  list(moy = moy, s = sqrt(s2), rsb = -10 * log10(s2))
}

# Décomposition de la variance par terme (somme des carrés ajustée)
ss_taguchi <- function(mm, assign, y, termes_idx) {
  n <- length(y)
  full <- ols_rang(mm, y)
  res <- vector("list", length(termes_idx));
  for (t in seq_along(termes_idx)) {
    cols <- which(assign %in% termes_idx[[t]])
    sans <- ols_rang(mm[, -cols, drop = FALSE], y)
    res[[t]] <- list(ss = sans$rss - full$rss, df = full$rang - sans$rang)
  }
  list(full = full, termes = res)
}
ols_rang <- function(M, y) {
  q <- qr(M, tol = 1e-8)
  fit <- stats::lm.fit(M, y)
  list(rss = sum(fit$residuals^2, na.rm = TRUE), rang = q$rank)
}

# Regroupement (mise en erreur) des termes les plus faibles
regrouper <- function(ss, df, ss_res, df_res, auto = TRUE, seuil_F = 2) {
  nt <- length(ss); poole <- rep(FALSE, nt)
  if (!auto) return(poole)
  repeat {
    e_ss <- ss_res + sum(ss[poole]); e_df <- df_res + sum(df[poole])
    libres <- which(!poole & df > 0)
    if (!length(libres)) break
    ms <- ss[libres] / df[libres]
    j <- libres[which.min(ms)]
    if (e_df < 1) { poole[j] <- TRUE; next }
    ve <- e_ss / e_df
    if (ve <= 0 || (ss[j] / df[j]) / ve < seuil_F) { poole[j] <- TRUE; next }
    break
  }
  poole
}

preparer_taguchi <- function(ctx, par) {
  reps <- as.character(unlist(par$reponses))
  if (!length(reps)) stop("Choisir au moins une colonne de réponses.", call. = FALSE)
  for (r in reps) if (!(r %in% names(ctx$plan))) stop(sprintf("Réponse « %s » introuvable.", r), call. = FALSE)
  type <- par$type_rsb %||% "plus_grand"
  Y <- as.matrix(ctx$plan[, reps, drop = FALSE]); mode(Y) <- "numeric"
  if (any(rowSums(!is.na(Y)) == 0)) stop("Chaque essai doit avoir au moins une réponse renseignée.", call. = FALSE)
  st <- rsb_essais(Y, type)
  k <- ctx$k; noms <- ctx$noms
  D <- data.frame(row.names = seq_len(nrow(Y)))
  for (j in seq_len(k)) D[[paste0("f", j)]] <- factor(as.character(ctx$plan[[noms[j]]]), levels = as.character(unlist(ctx$fac[[j]]$niveaux)))
  if (any(is.na(D))) stop("Niveaux de facteurs invalides dans le plan.", call. = FALSE)
  # termes : facteurs principaux puis interactions
  termes <- lapply(seq_len(k), function(j) list(label = noms[j], vars = j, formule = paste0("f", j)))
  for (it in ctx$meta$interactions %||% list()) {
    v <- as.integer(unlist(it$facteurs))
    termes[[length(termes) + 1]] <- list(label = paste0(noms[v[1]], " × ", noms[v[2]]), vars = v, formule = paste0("f", v[1], ":f", v[2]))
  }
  list(Y = Y, type = type, st = st, D = D, termes = termes, N = nrow(Y), reps = reps)
}

analyse_unitaire <- function(P, y, auto, alpha) {
  # ANOVA d'une variable (RSB ou moyenne) : retourne tableau, termes regroupés et modèle final
  termes <- P$termes; nt <- length(termes)
  form_full <- stats::as.formula(paste("y ~", paste(vapply(termes, `[[`, "", "formule"), collapse = " + ")))
  d <- cbind(P$D, y = y)
  mm <- stats::model.matrix(form_full, d)
  asg <- attr(mm, "assign")
  # une entrée par terme (indice dans la formule)
  idx <- as.list(seq_len(nt))
  sst <- sum((y - mean(y))^2)
  sc <- ss_taguchi(mm, asg, y, idx)
  ss <- vapply(sc$termes, `[[`, 0, "ss"); df <- vapply(sc$termes, `[[`, 0, "df")
  df_res <- P$N - sc$full$rang; ss_res <- sc$full$rss
  poole <- regrouper(ss, df, ss_res, df_res, auto)
  e_ss <- ss_res + sum(ss[poole]); e_df <- df_res + sum(df[poole])
  ve <- if (e_df >= 1) e_ss / e_df else NA_real_
  retenus <- which(!poole & df > 0)
  # modèle final
  rhs <- if (length(retenus)) paste(vapply(termes[retenus], `[[`, "", "formule"), collapse = " + ") else "1"
  fit <- stats::lm(stats::as.formula(paste("y ~", rhs)), d)
  F <- ifelse(poole | df == 0 | is.na(ve), NA_real_, (ss / df) / ve)
  p <- ifelse(is.na(F), NA_real_, stats::pf(F, df, e_df, lower.tail = FALSE))
  contrib <- 100 * ss / sst
  pure <- ifelse(poole | is.na(ve), NA_real_, 100 * (ss - df * ve) / sst)
  list(termes = termes, ss = ss, df = df, ss_res = ss_res, df_res = df_res, poole = poole, e_ss = e_ss, e_df = e_df, ve = ve, F = F, p = p,
       contrib = contrib, pure = pure, sst = sst, retenus = retenus, fit = fit, n_ret_df = sum(df[retenus]))
}

tableau_anova_t <- function(A, N) {
  nt <- length(A$termes)
  src <- character(0); ddl <- numeric(0); ss <- numeric(0); ms <- numeric(0); F <- numeric(0); p <- numeric(0); co <- numeric(0); cp <- numeric(0); reg <- character(0)
  for (i in order(-A$ss)) {
    if (A$df[i] == 0) next
    src <- c(src, A$termes[[i]]$label); ddl <- c(ddl, A$df[i]); ss <- c(ss, A$ss[i]); ms <- c(ms, A$ss[i] / A$df[i]); F <- c(F, A$F[i]); p <- c(p, A$p[i])
    co <- c(co, A$contrib[i]); cp <- c(cp, A$pure[i]); reg <- c(reg, if (A$poole[i]) "regroupé dans l'erreur" else "")
  }
  src <- c(src, "Erreur"); ddl <- c(ddl, A$e_df); ss <- c(ss, A$e_ss); ms <- c(ms, if (A$e_df >= 1) A$e_ss / A$e_df else NA); F <- c(F, NA); p <- c(p, NA)
  co <- c(co, 100 * A$e_ss / A$sst); cp <- c(cp, NA); reg <- c(reg, "")
  src <- c(src, "Total"); ddl <- c(ddl, N - 1); ss <- c(ss, A$sst); ms <- c(ms, NA); F <- c(F, NA); p <- c(p, NA); co <- c(co, 100); cp <- c(cp, NA); reg <- c(reg, "")
  list(Source = src, ddl = ddl, Somme_carrés = ss, Carré_moyen = ms, F = F, p_value = p, Contribution_pct = co, Contribution_pure_pct = cp, Remarque = reg)
}

tableau_reponses <- function(ctx, P, y) {
  k <- ctx$k; nmax <- max(vapply(seq_len(k), function(j) nlevels(P$D[[j]]), 1L))
  M <- matrix(NA_real_, k, nmax)
  for (j in seq_len(k)) { m <- tapply(y, P$D[[j]], mean); M[j, seq_along(m)] <- as.numeric(m) }
  delta <- apply(M, 1, function(r) max(r, na.rm = TRUE) - min(r, na.rm = TRUE))
  rang <- rank(-delta, ties.method = "min")
  tab <- list(Facteur = ctx$noms)
  for (l in seq_len(nmax)) tab[[paste0("Niveau_", l)]] <- M[, l]
  tab$Delta <- delta; tab$Rang <- rang
  list(tab = tab, M = M)
}

niveaux_optimaux <- function(ctx, P, rsb, A_rsb) {
  k <- ctx$k
  moy_niv <- lapply(seq_len(k), function(j) tapply(rsb, P$D[[j]], mean))
  choix <- vapply(seq_len(k), function(j) names(moy_niv[[j]])[which.max(moy_niv[[j]])], "")
  fixe <- rep(FALSE, k)
  # interactions significatives : meilleure combinaison
  its <- which(vapply(P$termes, function(t) length(t$vars) == 2, TRUE) & !A_rsb$poole)
  for (i in its[order(-A_rsb$ss[its])]) {
    v <- P$termes[[i]]$vars
    if (any(fixe[v])) next
    cell <- tapply(rsb, list(P$D[[v[1]]], P$D[[v[2]]]), mean)
    w <- which(cell == max(cell, na.rm = TRUE), arr.ind = TRUE)[1, ]
    choix[v[1]] <- rownames(cell)[w[1]]; choix[v[2]] <- colnames(cell)[w[2]]; fixe[v] <- TRUE
  }
  list(choix = choix, moy_niv = moy_niv)
}

prediction_cfg <- function(ctx, P, A, cfg) {
  nd <- data.frame(row.names = 1)
  for (j in seq_len(ctx$k)) nd[[paste0("f", j)]] <- factor(cfg[j], levels = levels(P$D[[j]]))
  as.numeric(stats::predict(A$fit, nd))
}

intervalle_confirmation <- function(A, N, alpha, r) {
  if (is.na(A$ve) || A$e_df < 1) return(c(NA_real_, NA_real_))
  neff <- N / (1 + A$n_ret_df)
  Fq <- stats::qf(1 - alpha, 1, A$e_df)
  c(conf = sqrt(Fq * A$ve * (1 / neff + 1 / r)), pop = sqrt(Fq * A$ve / neff))
}

analyser_taguchi <- function(ctx, par) {
  alpha <- ctx$alpha
  P <- preparer_taguchi(ctx, par)
  auto <- (par$regroupement %||% "auto") == "auto"
  r_conf <- as.numeric(par$n_confirmation %||% 3)
  k <- ctx$k; noms <- ctx$noms
  rsb <- P$st$rsb; moy <- P$st$moy
  A_r <- analyse_unitaire(P, rsb, auto, alpha)
  A_m <- analyse_unitaire(P, moy, auto, alpha)
  avert <- character(0)
  lignes <- c(sprintf("Type de rapport signal/bruit : %s.", TYPES_RSB[[P$type]]),
              sprintf("Réponses utilisées : %s (%d essais).", paste(P$reps, collapse = ", "), P$N),
              "Le rapport signal/bruit (RSB) est toujours à maximiser : le niveau optimal d'un facteur est celui dont le RSB moyen est le plus élevé.")
  if (length(P$reps) == 1 && P$type %in% c("plus_grand", "plus_petit")) avert <- c(avert, "Une seule réponse par essai : le RSB ne mesure que le niveau de la réponse, pas sa dispersion ; des répétitions ou des conditions de bruit sont recommandées.")
  if (any(unlist(ctx$meta$fictif %||% list()))) avert <- c(avert, "Un facteur à deux niveaux occupe une colonne à trois niveaux (niveau fictif) : l'orthogonalité est légèrement rompue, les effets sont estimés par régression.")
  tables <- list()
  tables[["Résultats_par_essai"]] <- c(list(Essai = ctx$plan$Ordre_essai), setNames(lapply(noms, function(n) ctx$plan[[n]]), noms), list(Moyenne = moy, Écart_type = P$st$s, RSB = rsb))
  tr <- tableau_reponses(ctx, P, rsb); tm <- tableau_reponses(ctx, P, moy)
  tables[["Tableau_des_réponses_—_RSB_moyen_par_niveau"]] <- tr$tab
  tables[["Tableau_des_réponses_—_moyenne_par_niveau"]] <- tm$tab
  tables[["ANOVA_du_RSB"]] <- tableau_anova_t(A_r, P$N)
  tables[["ANOVA_de_la_moyenne"]] <- tableau_anova_t(A_m, P$N)
  if (any(A_r$poole)) lignes <- c(lignes, sprintf("RSB : termes regroupés dans l'erreur (F < 2 ou degrés de liberté insuffisants) : %s.", paste(vapply(P$termes[A_r$poole], `[[`, "", "label"), collapse = ", ")))
  if (any(A_m$poole)) lignes <- c(lignes, sprintf("Moyenne : termes regroupés dans l'erreur : %s.", paste(vapply(P$termes[A_m$poole], `[[`, "", "label"), collapse = ", ")))
  sig_r <- vapply(P$termes[A_r$retenus], `[[`, "", "label")
  sig_r <- sig_r[!is.na(A_r$p[A_r$retenus]) & A_r$p[A_r$retenus] < alpha]
  sig_m <- vapply(P$termes[A_m$retenus], `[[`, "", "label")
  sig_m <- sig_m[!is.na(A_m$p[A_m$retenus]) & A_m$p[A_m$retenus] < alpha]
  lignes <- c(lignes, "", sprintf("Termes significatifs sur le RSB (α = %s) : %s.", fr(alpha, 3), if (length(sig_r)) paste(sig_r, collapse = ", ") else "aucun"),
              sprintf("Termes significatifs sur la moyenne (α = %s) : %s.", fr(alpha, 3), if (length(sig_m)) paste(sig_m, collapse = ", ") else "aucun"))
  if (P$type == "nominal_1") {
    aj <- setdiff(sig_m, sig_r)
    if (length(aj)) lignes <- c(lignes, sprintf("Facteurs d'ajustement de la moyenne (agissent sur la moyenne sans dégrader le RSB) : %s.", paste(aj, collapse = ", ")))
  }
  # niveaux optimaux
  op <- niveaux_optimaux(ctx, P, rsb, A_r)
  sig_fact <- vapply(seq_len(k), function(j) {
    t <- which(vapply(P$termes, function(t) identical(t$vars, j), TRUE))
    t <- t[1]; (!A_r$poole[t] && !is.na(A_r$p[t]) && A_r$p[t] < alpha) ||
      any(vapply(seq_along(P$termes), function(i) length(P$termes[[i]]$vars) == 2 && j %in% P$termes[[i]]$vars && !A_r$poole[i] && !is.na(A_r$p[i]) && A_r$p[i] < alpha, TRUE)) }, TRUE)
  tables[["Niveaux_optimaux"]] <- list(Facteur = noms, Niveau_optimal = op$choix,
                                       RSB_moyen_du_niveau = vapply(seq_len(k), function(j) as.numeric(op$moy_niv[[j]][op$choix[j]]), 0),
                                       Influence_sur_le_RSB = ifelse(sig_fact, "significative", "non significative"))
  lignes <- c(lignes, "", "Configuration optimale (RSB maximal) : " %+% paste(sprintf("%s = %s", noms, op$choix), collapse = " ; ") %+% ".")
  if (any(!sig_fact)) lignes <- c(lignes, sprintf("Pour les facteurs sans influence significative sur le RSB (%s), le niveau retenu est indicatif : choisir le moins coûteux.", paste(noms[!sig_fact], collapse = ", ")))
  # prévision au niveau optimal et à la configuration de référence
  ref <- vapply(seq_len(k), function(j) { v <- par$reference[[noms[j]]]; if (is.null(v) || !(as.character(v) %in% levels(P$D[[j]]))) levels(P$D[[j]])[1] else as.character(v) }, "")
  pr_o <- c(prediction_cfg(ctx, P, A_r, op$choix), prediction_cfg(ctx, P, A_m, op$choix))
  pr_r <- c(prediction_cfg(ctx, P, A_r, ref), prediction_cfg(ctx, P, A_m, ref))
  ic_r <- intervalle_confirmation(A_r, P$N, alpha, r_conf); ic_m <- intervalle_confirmation(A_m, P$N, alpha, r_conf)
  tables[["Prévision_au_niveau_optimal"]] <- list(
    Configuration = c("Optimale", "De référence", "Gain"),
    RSB_prévu = c(pr_o[1], pr_r[1], pr_o[1] - pr_r[1]),
    Moyenne_prévue = c(pr_o[2], pr_r[2], pr_o[2] - pr_r[2]))
  conf_txt <- c(sprintf("Intervalle de confirmation du RSB optimal (%d essai(s) de confirmation, risque %s) : %s ± %s.", as.integer(r_conf), fr(alpha, 3), fr(pr_o[1], 4), fr(ic_r[1], 3)),
                sprintf("Intervalle de confirmation de la moyenne optimale : %s ± %s.", fr(pr_o[2], 4), fr(ic_m[1], 3)))
  if (any(is.na(ic_r))) { conf_txt <- "Intervalle de confirmation non calculable : aucun degré de liberté pour l'erreur, même après regroupement."; }
  lignes <- c(lignes, "", sprintf("RSB prévu à l'optimum : %s (référence : %s, gain : %s dB).", fr(pr_o[1], 4), fr(pr_r[1], 4), fr(pr_o[1] - pr_r[1], 3)),
              sprintf("Moyenne prévue à l'optimum : %s (référence : %s).", fr(pr_o[2], 4), fr(pr_r[2], 4)), conf_txt,
              "Le modèle de prévision est additif sur les facteurs (et interactions) significatifs ; il doit être confirmé par des essais à la configuration optimale.")
  # figures
  figs <- c(fig_effets_taguchi(ctx, P, tr$M, rsb, "RSB moyen (dB)", "Effets principaux sur le rapport signal/bruit", "01_effets_rsb.png"),
            fig_effets_taguchi(ctx, P, tm$M, moy, "Moyenne", "Effets principaux sur la moyenne", "02_effets_moyenne.png"),
            fig_contributions(ctx, A_r, "03_contributions.png"))
  fi <- fig_interactions_taguchi(ctx, P, rsb, "04_interactions.png"); if (!is.null(fi)) figs <- c(figs, fi)
  concl <- sprintf("Optimum : %s ; RSB prévu %s dB.", paste(sprintf("%s=%s", noms, op$choix), collapse = ", "), fr(pr_o[1], 3))
  list(titre = sprintf("Plan de Taguchi — analyse (%s)", TYPES_RSB[[P$type]]), lignes = lignes, tables = tables, figures = figs, avertissements = avert, conclusion = concl,
       optimum = as.list(setNames(op$choix, noms)), tables_libres = list("Résultats_par_essai"))
}

prevoir_taguchi <- function(ctx, par) {
  P <- preparer_taguchi(ctx, par)
  auto <- (par$regroupement %||% "auto") == "auto"; alpha <- ctx$alpha
  A_r <- analyse_unitaire(P, P$st$rsb, auto, alpha); A_m <- analyse_unitaire(P, P$st$moy, auto, alpha)
  pts <- par$points; if (!length(pts)) stop("Aucune configuration à prévoir.", call. = FALSE)
  r_conf <- as.numeric(par$n_confirmation %||% 3)
  ic_r <- intervalle_confirmation(A_r, P$N, alpha, r_conf); ic_m <- intervalle_confirmation(A_m, P$N, alpha, r_conf)
  tab <- list(Configuration = seq_along(pts))
  cfgs <- lapply(pts, function(p) vapply(seq_len(ctx$k), function(j) {
    v <- p[[ctx$noms[j]]]; if (is.null(v) || !(as.character(v) %in% levels(P$D[[j]]))) stop(sprintf("Niveau invalide ou manquant pour « %s ».", ctx$noms[j]), call. = FALSE); as.character(v) }, ""))
  for (j in seq_len(ctx$k)) tab[[ctx$noms[j]]] <- vapply(cfgs, function(c) c[j], "")
  rs <- vapply(cfgs, function(c) prediction_cfg(ctx, P, A_r, c), 0); mo <- vapply(cfgs, function(c) prediction_cfg(ctx, P, A_m, c), 0)
  tab$RSB_prévu <- rs; tab$Marge_RSB <- rep(ic_r[1], length(rs)); tab$Moyenne_prévue <- mo; tab$Marge_moyenne <- rep(ic_m[1], length(rs))
  avert <- if (is.na(ic_r[1])) "Marges non calculables : aucun degré de liberté pour l'erreur." else character(0)
  list(titre = "Prévision — plan de Taguchi", tables = list(Prévisions = tab), figures = character(0), avertissements = avert,
       lignes = c(sprintf("Type de RSB : %s.", TYPES_RSB[[P$type]]),
                  sprintf("Marges : intervalle de confirmation (%d essai(s) de confirmation, risque %s).", as.integer(r_conf), fr(alpha, 3)),
                  "Modèle additif sur les termes significatifs du RSB et de la moyenne."),
       conclusion = sprintf("RSB prévu : %s.", paste(fr(rs, 4), collapse = " ; ")), tables_libres = list("Prévisions"))
}

# ---------------------------------------------------------------------------
#  Graphiques Taguchi
# ---------------------------------------------------------------------------
fig_effets_taguchi <- function(ctx, P, M, valeurs, ylab, titre, nom) {
  k <- ctx$k; dims <- grille_panneaux(k)
  ouvrir_png(ctx$dossier, nom, ctx$ech, 900, if (dims[1] > 2) 700 else 600)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mfrow = dims, mar = c(4, 3.8, 2.5, 0.8), oma = c(0, 1, 2, 0))
  yl <- range(M, na.rm = TRUE); yl <- yl + c(-0.08, 0.08) * max(diff(yl), 1e-9); moy <- mean(valeurs)
  for (j in seq_len(k)) {
    lv <- levels(P$D[[j]]); y <- M[j, seq_along(lv)]
    graphics::plot(seq_along(lv), y, type = "b", pch = 19, lwd = 2, col = COL_BLEU, ylim = yl, xlim = c(0.6, length(lv) + 0.4), xaxt = "n", xlab = "", ylab = if ((j - 1) %% dims[2] == 0) ylab else "", main = ctx$noms[j], cex.main = 1)
    graphics::axis(1, at = seq_along(lv), labels = lv); graphics::abline(h = moy, lty = 2, col = COL_GRIS); graphics::box(col = COL_GRIS)
  }
  graphics::mtext(titre, outer = TRUE, font = 2)
  nom
}

fig_contributions <- function(ctx, A, nom) {
  ok <- A$df > 0
  lab <- c(vapply(A$termes[ok], `[[`, "", "label")[!A$poole[ok]], "Erreur")
  v <- c(A$contrib[ok][!A$poole[ok]], 100 * A$e_ss / A$sst)
  o <- order(v, decreasing = TRUE); lab <- lab[o]; v <- v[o]
  ouvrir_png(ctx$dossier, nom, ctx$ech)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  wl <- max(graphics::strwidth(lab, units = "inches", cex = 0.9)) + 0.25
  graphics::par(mai = c(0.9, min(wl, 4), 0.8, 0.5))
  b <- graphics::barplot(rev(v), horiz = TRUE, names.arg = rev(lab), las = 1, col = rev(ifelse(lab == "Erreur", COL_GRIS, COL_BLEU)), border = NA, xlim = c(0, max(v) * 1.15), xlab = "Contribution à la variation du RSB (%)", main = "Contribution de chaque terme (ANOVA du RSB)")
  graphics::text(rev(v), b, paste0(gsub(".", ",", formatC(rev(v), format = "f", digits = 1), fixed = TRUE), " %"), pos = 4, cex = 0.9)
  nom
}

fig_interactions_taguchi <- function(ctx, P, rsb, nom) {
  its <- Filter(function(t) length(t$vars) == 2, P$termes)
  if (!length(its)) return(NULL)
  dims <- grille_panneaux(length(its))
  ouvrir_png(ctx$dossier, nom, ctx$ech, 900, if (dims[1] > 2) 700 else 600)
  op <- graphics::par(no.readonly = TRUE); on.exit({ graphics::par(op); fermer_png() })
  graphics::par(mfrow = dims, mar = c(4, 3.8, 2.5, 0.8), oma = c(0, 1, 2, 0))
  cols <- grDevices::hcl.colors(5, "Dark 3")
  cells <- lapply(its, function(t) tapply(rsb, list(P$D[[t$vars[1]]], P$D[[t$vars[2]]]), mean))
  yl <- range(unlist(cells), na.rm = TRUE); yl <- yl + c(-0.08, 0.08) * max(diff(yl), 1e-9)
  for (i in seq_along(its)) {
    t <- its[[i]]; cell <- cells[[i]]; a <- t$vars[1]; b <- t$vars[2]
    graphics::plot(NA, xlim = c(0.7, nrow(cell) + 0.3), ylim = yl, xaxt = "n", xlab = ctx$noms[a], ylab = "RSB moyen (dB)", main = t$label, cex.main = 1)
    graphics::axis(1, at = seq_len(nrow(cell)), labels = rownames(cell))
    for (l in seq_len(ncol(cell))) graphics::lines(seq_len(nrow(cell)), cell[, l], type = "b", pch = 19 + (l %% 3), lwd = 2, col = cols[l])
    graphics::legend("top", legend = paste(ctx$noms[b], "=", colnames(cell)), col = cols[seq_len(ncol(cell))], lty = 1, pch = 19, bty = "n", cex = 0.75)
    graphics::box(col = COL_GRIS)
  }
  graphics::mtext("Interactions étudiées — RSB moyen", outer = TRUE, font = 2)
  nom
}
