# Export the ei.Datasets New Zealand tables (same construction as eiPartialID::eiDatasetsTables)
# to two CSV files for the exploratory diagnostics in nz_diagnostics.py.
# Usage: Rscript export_nz_units.R <outdir>
args <- commandArgs(TRUE); outdir <- if (length(args)) args[1] else "."
src <- parse("package/eiPartialID/R/rxc_api.R")
i <- which(sapply(src, function(e) is.call(e) && identical(as.character(e[[2]]), "eiDatasetsTables")))
eval(src[[i]])
tabs <- eiDatasetsTables()
units <- list(); truth <- list()
for (k in names(tabs)) {
  tb <- tabs[[k]]
  if (length(tb$rows) != 3 || length(tb$cols) != 3) next          # keep the 3x3 tables
  d <- tb$data[, c(tb$rows, tb$cols)]
  names(d) <- c("x_LAB", "x_NAT", "x_OTH", "t_LAB", "t_NAT", "t_OTH")
  units[[k]] <- cbind(table = k, election = tb$meta$election, district = tb$meta$district, d)
  tr <- tb$truth
  truth[[k]] <- data.frame(table = k, row = rep(c("LAB", "NAT", "OTH"), 3),
                           col = rep(c("LAB", "NAT", "OTH"), each = 3), count = as.vector(tr))
}
write.csv(do.call(rbind, units), file.path(outdir, "nz_units.csv"), row.names = FALSE)
write.csv(do.call(rbind, truth), file.path(outdir, "nz_truth.csv"), row.names = FALSE)
cat("tables:", length(units), "\n")
