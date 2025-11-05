summary2 <- function(input_data, rotate=FALSE, digits=3) {

  get_summary <- function(input, colname='summary', digits=3) {
    cs1 <- summary(input)
    names(cs1) <- NULL
    cs2 <- data.frame('summary'=c(
      'count' = round(length(input), 0), 
      'distinct' = round(length(unique(input)), digits),
      #'mean' = round(mean(input), digits),
      'mean' = round(cs1[4], digits), 
      'std' = round(sd(input), digits), 
      'min' = round(cs1[1], digits),
      'q1' = round(cs1[2], digits),
      'median' = round(cs1[3], digits),
      'q3' = round(cs1[5], digits),
      'max' = round(cs1[6], digits)
    ))
    colnames(cs2) <- c(colname)
    cs2
  }
    
  # Initialize report
  report <- data.frame()

  if ( is.null(colnames(input_data)) == FALSE ) {
    for (column_name in colnames(input_data)){
      # Only generate a summary for integer or numeric values
      if (class(input_data[, column_name]) %in% c('integer', 'numeric')){
        cs2 <- get_summary(input_data[, column_name], colname=column_name, digits=digits)
        # For the first variable, define the report dataframe
        if (dim(report)[1]==0){ 
          report <- cs2 
          colnames(report) <- column_name
        } else { 
          report[, column_name] <- cs2 
        }
      }
    } 
  } else {
    c2 <- get_summary(input_data, digits=digits)
    # For the first variable, define the report dataframe
    report <- c2
  }
  
  if ( rotate ){ t(report) } else { report }

}
