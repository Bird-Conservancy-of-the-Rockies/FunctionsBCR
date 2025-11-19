speciesNameCase <- function(rawNames) {
  require(stringr)
  
  # output vector
  out <- character(length(rawNames))
  
  for (i in seq_along(rawNames)) {
    s <- rawNames[i]
    
    # lower-case everything first
    s <- str_to_lower(s)
    
    # split into up to two parts
    parts <- str_split(s, "\\s+", n = 2)[[1]]
    
    # --- First word rules ---
    fix_first <- function(w) {
      w <- str_replace(w, "^[a-z]", toupper)             # capitalize first letter
      w <- str_replace_all(w, "(?<=-)[a-z]", tolower)    # after hyphen → lower
      w <- str_replace_all(w, "(?<=['])[a-z]", tolower)  # after apostrophe → lower
      w
    }
    
    # --- Second word rules ---
    fix_second <- function(w) {
      w <- str_replace(w, "^[a-z]", toupper)             # capitalize first letter
      w <- str_replace_all(w, "(?<=-)[a-z]", toupper)    # after hyphen → UPPER
      w <- str_replace_all(w, "(?<=['])[a-z]", tolower)  # after apostrophe → lower
      w
    }
    
    # Apply rules depending on number of words
    if (length(parts) == 1) {
      fixed <- fix_first(parts[1])
    } else {
      fixed <- paste(fix_first(parts[1]), fix_second(parts[2]))
    }
    
    out[i] <- fixed
  }
  
  out
}
