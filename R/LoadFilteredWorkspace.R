LoadFilteredWorkspace <- function(workspace, objects.exclude = c(), target.envir) {
  
  # Create a temporary environment
  temp_env <- new.env()
  
  # Load the .RData file into the temporary environment
  load(workspace, envir = temp_env)
  
  # Remove unwanted objects from the temporary environment
  rm(list = objects.exclude, envir = temp_env)
  
  # Copy remaining objects to the global environment
  list2env(as.list(temp_env), envir = target.envir)
  
  # Clean up
  rm(temp_env)
}
