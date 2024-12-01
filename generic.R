# define fit generic
fit <- function(object, ...) {
  UseMethod("fit")
}
# define predict generic
predict <- function(object, ...) {
  UseMethod("predict")
}
