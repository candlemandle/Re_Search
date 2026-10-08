# Load Code A from the project root when tests run from tests/testthat.
code_a_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(code_a_root, "R", "code_a_config.R"))
source_code_a(code_a_root)
