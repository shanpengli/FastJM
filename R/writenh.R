writenh<-function(newdata,newfile) {
  rowss=apply(newdata, 1, paste, collapse=" ")
  for (i in 1:nrow(newdata)-1)
    rowss[i]=paste0(rowss[i],'\n')
  writeLines(rowss, sep="", newfile)
}
