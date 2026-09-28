#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else
  normalizePath(".", mustWork = TRUE)
installer <- file.path(root, "scripts", "R", "task004d_install_deps.R")
if (!file.exists(installer)) stop("Task 004d dependency installer not found: ", installer)
source(installer, local = globalenv())
