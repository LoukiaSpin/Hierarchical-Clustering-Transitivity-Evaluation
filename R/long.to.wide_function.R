#*******************************************************************************
#*
#*
#*            Turn long 'netmeta' format into wide 'gemtc' format                                                                           
#*
#*
#*******************************************************************************

long_to_wide <- function (input) {
  
  
  ## Check defaults
  # Dataset
  input <- if (any(sapply(input, typeof)[1:3] != "character")) {
    stop("The first three columns (trial and arms) must be 'characters'.", 
         call. = FALSE)
  } else if (any(sapply(input, typeof)[-c(1:3)] == "character")) {
    stop("The characteristics must be 'double'.", call. = FALSE)
  } else {
    input
  }
  colnames(input) <- c("trial", "arm1", "arm2", "sample1", "sample2")

  
  # Sort dataset by table names
  data_network_new <- input[order(input$trial),]
    #input[match(input$trial, names(table(input$trial))),]
  
  # And an indicator of the number of arms
  data_network_new$na_arm <- rep(table(input$trial), table(input$trial))
  
  # Split dataset by the indicator
  split_data <- split(data_network_new, f = data_network_new$na_arm)
  
  # Turn arms and samples columns into one vector each
  arm_long <- 
    do.call(rbind, 
            lapply(split_data, function(x) melt(x, id = "trial", measure.vars = c("arm1", "arm2"))))
  sample_long <- 
    do.call(rbind, 
            lapply(split_data, function(x) melt(x, id = "trial", measure.vars = c("sample1", "sample2"))))
  
  # Bring both in a data-frame
  combine <- cbind(arm_long, sample_long)[, -c(2, 4:5)]
  colnames(combine) <- c("study", "t", "n")
  
  # Remove duplicates according to indicator
  ind <- paste(combine$study, combine$t)
  combine_new0 <- combine[!duplicated(ind), ]
  
  # Order the trial so that pairs come sequentially
  combine_new <- combine_new0[order(combine_new0$study), ]
  
  # Add an indicator of number of arms per trial
  combine_new$id <- unlist(lapply(data.frame(table(combine_new$study))[, 2], function(x) 1:x))
  
  # Turn into wide
  network_data0 <- reshape(combine_new, 
                           v.names = c("t", "n"), 
                           idvar = "study", 
                           timevar = "id", 
                           direction = "wide")
  
  # Select columns by specific prefix
  treat_cols <- network_data0[, startsWith(colnames(network_data0), "t.")]
  sample_cols <- network_data0[, startsWith(colnames(network_data0), "n.")]
  
  # Get the finalised dataset
  network_data <- cbind(study = network_data0$study, treat_cols, sample_cols)
  
  return(network_data)
}