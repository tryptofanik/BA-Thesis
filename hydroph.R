#! /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

library(Peptides)

con = file(args[1], "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    seq = unlist(strsplit(line, "\t"))[1]
      nonidp_z = unlist(strsplit(line, "\t"))[2]
      #idp_z = unlist(strsplit(line, "\t"))[3]
      cat(seq, "\t")
      cat(nonidp_z, "\t")
      #cat(idp_z, "\t")
      cat(hydrophobicity(seq, scale = "Guy"), "\n")
  }