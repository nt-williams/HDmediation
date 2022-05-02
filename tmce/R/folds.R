make_folds <- function(data, V) {
    folds <- origami::make_folds(data, V = V)
    if (V == 1) {
        folds[[1]]$training_set <- folds[[1]]$validation_set
    }
    folds
}
