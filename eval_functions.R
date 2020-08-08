# Model analysis functions ####
# Robert Schmidt
# Description: helpful functions for model evaluation & analysis

# Libraries
library(tidyverse) # dplyr & ggplot
library(corrplot) # correlation plots
library(polycor) # for non-continuous variable correlations
library(caret) # pre-defined evaluation functions


# Variable correlations ####
make_corrPlot = function(df, col_labels = FALSE, text_vert_angle = 60, text_size = 0.5,
                         var_clusters = 5, insig.option = "n", sig.level = 0.05,
                         spearman_rank = FALSE) {
  
  # Input: data frame
  # Output: correlation plot (visualization of the correlation matrix)
  
  # Input explanations:
  # Text: text_vert_angle = angle of column text; text_size = size of text in [0, 1]
  # Variable clusters: var_clusters = number of variable clusters in corrplot hierarchical clustering
  # Significance: insig.option = "n" (no indication), "blank" (no plot if not significant); sig.level = threshold
  # Spearman rank: test for monotonicity rather than linearity; helpful for percent/bounded predictors
  
  if (spearman_rank == FALSE) {
    hetcor.df = hetcor(df) # hetcor adapts correlation matrix if not cts vs. cts variables
    cors = hetcor.df$correlations
    
    if (col_labels == FALSE) {
      colnames(cors) = c() # useful if the variable names are long
    }
    
    p.vals = hetcor.df$tests
    
    
  }
  
  
}