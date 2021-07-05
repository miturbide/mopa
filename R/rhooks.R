
#'@references http://thecoatlessprofessor.com/programming/r-data-packages-in-external-data-repositories-using-the-additional_repositories-field/

.onLoad <- function(...) {
      repos = getOption("repos")
      repos["datarepo"] = "http://miturbide.github.io/datarepo"
      options(repos = repos)
      invisible(repos)
}



