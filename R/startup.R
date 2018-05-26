# .onAttach is executed upon attaching the TargetedQCMS packing and prints the package startup message to the console

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("TargetedMSQC provides a semi-automated workflow for quality control (QC) of chromatographic peaks in targeted proteomics experiments using supervised learning.")
}

# .onLoad is executed upon loading the TargetedQCMS packing

.onLoad <- function(libname, pkgname) {

}
