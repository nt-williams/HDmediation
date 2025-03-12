make_folds <- function(data, V, id = NULL) {
    folds <- origami::make_folds(data, V = V, cluster_ids = id)
    if (V == 1) {
        folds[[1]]$training_set <- folds[[1]]$validation_set
    }
    folds
}
