# Resolve the dataset independently from the Git repository location.
#
# Set MARINE_P_DATA_DIR in a project-local .Renviron when large data files live
# on another disk. Without that variable, analyses use the repository's data/
# directory, which keeps the workflow portable for other installations.
marine_project_dir <- function() {
  working_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- unique(c(working_dir, dirname(working_dir)))
  marker <- "marine-p-deficiency-metat.rproj"
  matches <- candidates[file.exists(file.path(candidates, marker))]

  if (!length(matches)) {
    stop("Cannot locate the Marine P project root from: ", working_dir)
  }

  matches[[1]]
}

PROJECT_DIR <- marine_project_dir()

# R only reads .Renviron automatically from its startup directory. Explicitly
# read the project-local file as well so command-line runs from analysis/ work.
local_renviron <- file.path(PROJECT_DIR, ".Renviron")
if (!nzchar(Sys.getenv("MARINE_P_DATA_DIR", unset = "")) &&
    file.exists(local_renviron)) {
  readRenviron(local_renviron)
}

marine_data_dir <- function(project_dir = PROJECT_DIR) {
  configured_dir <- Sys.getenv("MARINE_P_DATA_DIR", unset = "")

  if (nzchar(configured_dir)) {
    data_dir <- path.expand(configured_dir)
  } else {
    data_dir <- file.path(project_dir, "data")
  }

  if (!dir.exists(data_dir)) {
    stop(
      "Marine P data directory does not exist: ", data_dir,
      "\nSet MARINE_P_DATA_DIR in the project-local .Renviron file."
    )
  }

  normalizePath(data_dir, winslash = "/", mustWork = TRUE)
}

DATA_DIR <- marine_data_dir()
data_file <- function(...) file.path(DATA_DIR, ...)
