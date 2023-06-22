# https://rdrr.io/bioc/compEpiTools/src/R/simplifyGOterms.R

#source('src_simplifyGOterms.R')
#simple.gos<-simplifyGOterms(goterms=sig.go.terms$GO, maxOverlap= 0.01, ontology='BP', go2allEGs= org.Mm.egGO2ALLEGS)
#length(simple.gos)
#sapply(simple.gos,function(i)GOTERM[[i]]@Term)

simplifyGOterms <- function(goterms, maxOverlap=0.8, ontology, go2allEGs) {
    if(!is.character(goterms))
        stop('goterms has to be of class character ...')
    if(!is.numeric(maxOverlap))
        stop('maxOverlap has to be of class numeric ...')
    if(maxOverlap < 0 || maxOverlap > 1)
        stop('maxOverlap is a percentage and has to range in [0,1] ...')
    if(!all(ontology %in% c('BP','CC','MF')))
        stop('ontology has to be one of: CC, BP, MF ...')
    if(!is(go2allEGs,"AnnDbBimap"))
        stop('go2allEGs has to be of class AnnDbBimap ..')

    if(ontology == 'CC') go2parents <- as.list(GOCCPARENTS)
    if(ontology == 'BP') go2parents <- as.list(GOBPPARENTS)
    if(ontology == 'MF') go2parents <- as.list(GOMFPARENTS)
    go2discard <- NULL
    for(goterm in goterms) {
        parents <- go2parents[[goterm]]
        parents <- intersect(parents, goterms)
                                        # no parents are found for a given GO term, check the others
        if(length(parents) == 0) next
        ##gotermEGs <- go2allEGs[[goterm]] # EGs associated to a given GO term
        # EGs associated to a given GO term
        gotermEGs <- unlist(mget(goterm, go2allEGs))
        for(parent in parents) {
            parentEGs <- go2allEGs[[parent]] # EGs associated to its parent
            # EGs associated to its parent
            parentEGs <- unlist(mget(parent, go2allEGs))
            commonEGs <- intersect(gotermEGs, parentEGs)
            if(length(commonEGs) / length(parentEGs) > maxOverlap)
                go2discard <- c(go2discard, parent)
        }
    }
                                        # discard redundant parents
    if(length(go2discard) > 0)
        goterms <- goterms[-which(goterms %in% go2discard)]

    return(goterms)
}