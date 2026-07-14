Impute_missing_covs_rf <- function(dat, v.fill, v.inform, v.factor = NULL) {
  require(randomForest)
  require(dplyr)
  require(stringr)

  # Define data objects and convert to data frames as needed.
  dat.mat <- is.matrix(dat)
  if(dat.mat) dat <- data.frame(dat)

  # Convert flagged columns to factors for classification-mode RF.
  # Record original types (and level order, for columns already stored as
  # factors) so they can be restored afterward.
  if(!is.null(v.factor)) {
    v.factor <- unique(v.factor[v.factor %in% names(dat)])
    bad.classes <- sapply(dat[, v.factor, drop = FALSE], function(x) length(class(x)) > 1)
    if(any(bad.classes)) {
      stop("v.factor columns must be simple atomic types (numeric/integer/character/factor); ",
           "unsupported type(s) found for: ", str_c(v.factor[bad.classes], collapse = ", "))
    }
    orig.class <- sapply(dat[, v.factor, drop = FALSE], function(x) class(x)[1])
    orig.levels <- lapply(dat[, v.factor, drop = FALSE], function(x) if(is.factor(x)) levels(x) else NULL)
    dat[, v.factor] <- lapply(dat[, v.factor, drop = FALSE], factor)
  }

  ind.missing <- which(is.na(dat[,v.fill]))
  ind.known <- which(!is.na(dat[,v.fill]))
  rf <- randomForest(as.formula(str_c(v.fill, "~",
                                      str_c(v.inform[which(v.inform != v.fill)],
                                            collapse = "+"))),
                     data = (dat %>% slice(ind.known)))
  predicted <- predict(rf, newdata = (dat %>% slice(ind.missing)))

  # predicted carries the same factor levels as dat[, v.fill] (whether v.fill
  # was already a factor, or converted via v.factor), so it assigns directly
  # with no type coercion needed.
  dat[ind.missing, v.fill] <- predicted

  # Restore v.factor columns to their original types (and level order, for
  # columns that were already factors).
  if(!is.null(v.factor)) {
    for(v in v.factor) {
      if(identical(orig.class[[v]], "factor")) {
        dat[, v] <- factor(as.character(dat[, v]), levels = orig.levels[[v]])
      } else {
        dat[, v] <- as(as.character(dat[, v]), orig.class[[v]])
      }
    }
  }

  if(dat.mat) dat <- data.matrix(dat)
  return(dat)
}
