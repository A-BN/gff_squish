cli::boxx(label = crayon::bold("Wilkommen bienvenue!"))
args <- commandArgs(trailingOnly = TRUE)

arga <- list()
for (i in seq(1, length(args), 2)) {
  arga[args[i]] <- args[i + 1]
}

print(arga)
file_in <- arga[["--gff_in"]]

profvis::profvis({
# file_in <- file("to_squish.gff", open = "r")
file_in <- file("11111_a_trim.gff", open = "r")
file_out <- file("gff_out_R.gff", open = "w")
# n_line <- 0
contig_list <- list()
to_add <- 0

while (TRUE) {
  linus <- readLines(con = file_in, n = 1)
  # n_line <- n_line + 1
  if (linus == "##FASTA") break
}

while (TRUE) {
  linus <- readLines(con = file_in, n = 1)
  if (length(linus) == 0) break
  if (startsWith(x = linus, prefix = ">")) {
    contig_name <- linus
    contig_list[[linus]] <- to_add
  }
  else {
    to_add <- to_add + nchar(linus)
  }
}

seek(con = file_in, where = 0)

while (TRUE) {
  linus <- readLines(con = file_in, n = 1)
  if (linus == "##FASTA") {
    writeLines(text = linus, con = file_out)
    break
  }
  if (startsWith(x = linus, prefix = "##")) next
  annot <- unlist(strsplit(x = linus, split = "\t", fixed = TRUE))
  annot[4] <- as.numeric(annot[4]) + contig_list[[paste0(">", annot[1])]]
  annot[5] <- as.numeric(annot[4]) + contig_list[[paste0(">", annot[1])]]
  annot[1] <- "squouich"
  writeLines(text = paste(annot, collapse = "\t"), con = file_out)
}



while (TRUE) {
  linus <- readLines(con = file_in, n = 1)
  if (length(linus) == 0) break
  if (!startsWith(x = linus, prefix = ">")) {
    writeLines(text = linus, con = file_out)
  }
}

close(file_in)
close(file_out)
})

rm(list = ls())