# -*- coding: utf-8 -*-
# =============================================================================
#  cli_doe.R — Passerelle en ligne de commande vers le moteur R des plans d'expériences
#  Usage : Rscript cli_doe.R <dossier de travail>
#    entrée : <dossier>/params.json (commande + paramètres), <dossier>/plan.csv (plan et réponses, si nécessaire)
#    sortie : <dossier>/resultat.json et les graphiques PNG
# =============================================================================
suppressMessages(library(jsonlite))
args_tous <- commandArgs(trailingOnly = FALSE)
f_arg <- sub("^--file=", "", args_tous[grep("^--file=", args_tous)][1])
rep_script <- dirname(normalizePath(f_arg))
for (f in c("plans_doe.R", "moteur_doe.R", "graphiques_doe.R", "taguchi_doe.R"))
  source(file.path(rep_script, f), encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)
dossier <- args[1]
if (is.na(dossier) || !dir.exists(dossier)) stop("Dossier de travail manquant.")
par <- fromJSON(file.path(dossier, "params.json"), simplifyVector = FALSE)

resultat_plan <- function(r) {
  list(ok = TRUE, titre = r$titre, lignes = r$lignes, tables = r$tables, plan = as.list(r$plan), meta = r$meta, figures = character(0),
       avertissements = character(0), conclusion = "", tables_libres = list())
}

executer <- function(par, dossier) {
  cmd <- par$commande
  if (is.null(cmd)) stop("Commande manquante.", call. = FALSE)
  if (cmd == "plan_criblage") return(resultat_plan(plan_criblage(par)))
  if (cmd == "plan_factoriel") return(resultat_plan(plan_factoriel(par)))
  if (cmd == "plan_surface") return(resultat_plan(plan_surface(par)))
  if (cmd == "plan_melange") return(resultat_plan(plan_melange(par)))
  if (cmd == "plan_taguchi") return(resultat_plan(plan_taguchi(par)))
  if (cmd == "verifier_r") return(list(ok = TRUE, titre = "R prêt", lignes = list(R.version.string)))
  ctx <- charger_contexte(dossier, par)
  r <- switch(cmd,
    analyse = if (ctx$famille == "taguchi") analyser_taguchi(ctx, par) else analyser_modele(ctx, par),
    prevision = if (ctx$famille == "taguchi") prevoir_taguchi(ctx, par) else prevoir_doe(ctx, par),
    optimisation = optimiser_doe(ctx, par),
    surfaces = surfaces_doe(ctx, par),
    stop(sprintf("Commande inconnue : %s", cmd), call. = FALSE))
  r$ok <- TRUE; r$commande <- cmd
  r
}

res <- tryCatch(executer(par, dossier), error = function(e) list(ok = FALSE, erreur = conditionMessage(e)))
if (is.null(res$commande)) res$commande <- par$commande
if (!is.null(par$plan_id)) res$plan_id <- par$plan_id
txt <- toJSON(res, auto_unbox = TRUE, na = "null", null = "null", digits = NA, force = TRUE)
writeLines(txt, file.path(dossier, "resultat.json"), useBytes = TRUE)
