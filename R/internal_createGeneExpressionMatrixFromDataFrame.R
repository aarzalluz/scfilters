.createGeneExpressionMatrixFromDataFrame <- function(myDataFrame) {
    
    if(!is.data.frame(myDataFrame)) {
        stop(paste(
            deparse(substitute(myDataFrame)),
            "is not a data.frame."
        ))
    }
    
    firstColumn <-  select(myDataFrame, 1) %>% unlist(use.names = FALSE)
    
    if(is.character(firstColumn) | is.factor(firstColumn)) {
        if(any(duplicated(firstColumn))) {
            warning(paste(
                "Non unique gene names in first column of",
                deparse(substitute(myDataFrame))
            ))
        }
        outputMatrix <- as.matrix(myDataFrame[, -1, drop = FALSE])
        row.names(outputMatrix) <- make.names(firstColumn, unique = TRUE)
    } else {
        outputMatrix <- as.matrix(myDataFrame)
    }
    
    if(!is.numeric(outputMatrix)) {
        stop(paste(
            "Could not tranform",
            deparse(substitute(myDataFrame)),
            "into a numeric matrix."
        ))
    }
    
    return(outputMatrix)
}

