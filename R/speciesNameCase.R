speciesNameCase <- function(rawNames) {
  require(stringr)
  
  # ---- helper functions ----
  
  # First-word rule: hyphen → lower
  fix_first <- function(w) {
    w <- str_replace(w, "^[a-z]", toupper)
    w <- str_replace_all(w, "(?<=-)[a-z]", tolower)
    w <- str_replace_all(w, "(?<=['])[a-z]", tolower)
    w
  }
  
  # Middle-word rule (for >2 words): also hyphen → lower
  fix_middle <- function(w) {
    w <- str_replace(w, "^[a-z]", toupper)
    w <- str_replace_all(w, "(?<=-)[a-z]", tolower)
    w <- str_replace_all(w, "(?<=['])[a-z]", tolower)
    w
  }
  
  # Last-word rule: hyphen → UPPER
  fix_last <- function(w) {
    w <- str_replace(w, "^[a-z]", toupper)
    w <- str_replace_all(w, "(?<=-)[a-z]", toupper)
    w <- str_replace_all(w, "(?<=['])[a-z]", tolower)
    w
  }
  
  # ---- main loop ----
  out <- character(length(rawNames))
  
  for (i in seq_along(rawNames)) {
    s <- rawNames[i]
    s <- str_to_lower(s)
    
    parts <- str_split(s, "\\s+")[[1]]
    n <- length(parts)
    
    if (n == 1) {
      # Single word → treat like first word
      fixed <- fix_first(parts[1])
      
    } else if (n == 2) {
      # Two words → first + last rules
      fixed <- paste(fix_first(parts[1]),
                     fix_last(parts[2]))
      
    } else {
      # Three or more words
      first <- fix_first(parts[1])
      middles <- vapply(parts[2:(n-1)], fix_middle, FUN.VALUE = character(1))
      last <- fix_last(parts[n])
      
      fixed <- paste(c(first, middles, last), collapse = " ")
    }
    
    out[i] <- fixed
  }
  
  out
}
