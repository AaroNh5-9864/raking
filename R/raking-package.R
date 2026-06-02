#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import S7
## usethis namespace: end
NULL

# S7 methods defined inside a package must be registered when the package is
# loaded. Without this hook, method dispatch would fail.
.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}
