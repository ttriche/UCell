#' Check signature names and add standard names if missing
#' 
#' @param   features    the signature feature names
#' 
#' @return              cleaned names
check_signature_names <- function(features) {
  defaultSigName <- paste0(rep("signature_", length(features)),
                           seq_along(features))
  if(is.null(names(features))){
    names(features) <- defaultSigName
  } else {
    invalidNames <- names(features) == "" | duplicated(names(features))
    names(features)[invalidNames] <- defaultSigName[invalidNames]
  }
  return(features)
}


