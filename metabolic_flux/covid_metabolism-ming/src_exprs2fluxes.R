
exprs2fluxes <- function(model, x, na.before=NA, na.after=0) {
  # map a numeric vector x of expression levels of genes to the flux levels of reactions in the model
  # x should contain either continuous expression values (then the output will also be continuous) or discrete values from {-1,0,1,NA} representing low/medium/high expression levels
  # x should be named by gene symbols as used in model$genes, or if it's unnamed and length being length(model$genes), assume it's already in the same order as model$genes
  # NA's in x will be kept and propagated during the conversion to flux values by default; or specify na.before and NA will be changed to this value
  # NA's in the result flux values will be changed to 0 by default; or specify na.after and NA will be changed to that value (or na.after=NA to keep NA)

  if (is.null(names(x))) {
    if (length(x)==length(model$genes)) {
      message("exprs2fluxes(): Assuming the input vector is in the same order as model$genes.")
    } else stop("Input vector and model$genes have different lengths!")
  } else {
    x <- x[model$genes]
    if (all(is.na(x))) stop("Input doesn't contain any of the model genes!")
  }
  x[is.na(x)] <- na.before

  if (all(x %in% c(-1,0,1,NA))) {
    `&` <- function(a,b) {
      if (isTRUE((is.na(a) && b==-1) || (is.na(b) && a==-1))) return(-1) # if one is NA and the other is -1, for sure the result is -1; all other NA cases are undetermined and NA will be returned
      min(a,b)
    }
    `|` <- function(a,b) {
      if (isTRUE((is.na(a) && b==1) || (is.na(b) && a==1))) return(1) # if one is NA and the other is 1, for sure the result is 1; all other NA cases are undetermined and NA will be returned
      max(a,b)
    }
  } else {
    `&` <- function(a,b) min(a,b)
    `|` <- function(a,b) max(a,b)
  }

  res <- sapply(1:length(model$rules), function(i) {
    ri=unlist(model$rules[[i]]) #FIX
    ri=gsub('\\(','\\[',ri)  #FIX
    ri=gsub('\\)','\\]',ri) #FIX
    #if (ri=="") {
    if(length(ri)==0){ #FIX
      return(NA)
    } else {
      return(tryCatch(eval(parse(text=ri)),
        error=function(e) {
          warning("Failed parsing rule #", i)
          NA
      }))
    }
  })
  res[is.na(res)] <- na.after
  res[model$rules==""] <- NA # still NA for rxns w/o genes
  unname(res)
}

