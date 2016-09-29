
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    'This is BioGeoBias version ', utils::packageVersion("BioGeoBias"),
    '\nsee <url> for detailed description of the available methods.',
    '\ncontact laurensgeffert@gmail.com for bug reports, feedback, and suggestions')
}
