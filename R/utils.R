### =========================================================================
### Helper functions not exported
### =========================================================================

#supposed to merge paths and file irrespective OS and presence of trailing slash
.mergeDirAndFilename <- function(dir, files){
	#check for presence of / in filename in that case remove
	files <- sub("/","",files)
	paste(normalizePath(dir),"/",files, sep="")
}

