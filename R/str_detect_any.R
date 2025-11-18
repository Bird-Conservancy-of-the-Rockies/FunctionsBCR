str_detect_any <- function(string, pattern.vec, ignore_case = FALSE) {
  require(stringr)
  require(dplyr)
  
  p <- pattern.vec[1]
  out <- str_detect(string, ifelse(ignore_case, paste0("(?i)", p), p))
  if(length(pattern.vec) > 1) {
    for(p in pattern.vec[-1]) {
      out2 <- str_detect(string, ifelse(ignore_case, paste0("(?i)", p), p))
      out <- cbind(out, out2) %>%
        apply(1, any)
    }
  }
  return(out)
}