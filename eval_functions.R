# Model analysis functions ####
# Robert Schmidt
# Description: helpful functions for model evaluation & analysis

# Libraries
library(tidyverse) # dplyr & ggplot
library(corrplot) # correlation plots
library(polycor) # for non-continuous variable correlations
library(caret) # pre-defined evaluation functions
library(data.table) # for fast df manipulation


# Variable correlations ####

# Function: make correlation plot based on useful presets, rather than manual entry into corrplot
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
    
    corrplot(cors,
             p.mat = p.vals,
             sig.level = sig.level,
             insig = insig.option,
             method = "ellipse",
             tl.cex = text_size,
             tl.col = "black",
             tl.srt = text_vert_angle,
             order = "hclust",
             addrect = var_clusters,
             diag = FALSE)
  } else if (spearman_rank == TRUE) {
    cors = cor(df,
               use = "complete.obs",
               method = "spearman") # using only non-NA obs
    
    if (col_labels == FALSE) {
      colnames(cors) = c()
    }
    
    corrplot(cors,
             method = "ellipse",
             tl.cex = text_size,
             tl.col = "black",
             tl.srt = text_vert_angle,
             order = "hclust",
             addrect = var_clusters,
             diag = FALSE)
  }
}

# Function: find correlations of input and "flatten"/rank them in descending order
find_corr = function(input, method.input = "pearson") {
  # Input: df and method of correlation
  # Output: pairwise correlations, sorted in descending order by absolute value
  
  corr_matrix = cor(input, method = method.input)
  corr_matrix[lower.tri(corr_matrix, diag = TRUE)] = NA
  corr_matrix = as.data.frame(as.table(corr_matrix))
  corr_matrix = na.omit(corr_matrix)
  corr_matrix = corr_matrix[order(-abs(corr_matrix$Freq)), ]
  names(corr_matrix = c("predictor_1", "predictor_2", "correlation"))
  corr_matrix
}


# Model evaluation ####

# Function: plot percentile-scored values of predictions vs. actual response for binomial data & return correlation
ntile_eval = function(pred_rate, actual_counts, sample_weights, method = "preds") {
  # Note: this evaluation is for binomially distributed output (as used in binomial GLM)
  # Input: predictions, actual response, sample weights (the "n"/sample size for binomial data)
  # Option: bin by percentile of actual value or predicted value
  
  # Plot output: plot of values aggregated by percentile
  # Numeric output: correlation between predictions and actual response when aggregated by percentile
  
  df = data.frame(pred_rate, actual_counts, sample_weights)
  df$pred_counts = pred_rate*sample_weights
  
  if (method == "preds") {
    df$percentile = ntile(df$pred_rate, 100)
  } else if (method == "actual") {
    df$percentile = ntile(df$actual/df$sample_weights, 100)
  } else {
    print("Not a valid method!")
  }
  
  # Aggregate by percentile
  df_agg = aggregate(cbind(pred_counts, actual_counts, sample_weights) ~ percentile,
                     data = df,
                     sum)
  
  # Calculate rates when binned by percentile
  df_agg$actual_rate = df_agg$actual_counts/df_agg$sample_weights
  df_agg$pred_rate = df_agg$pred_counts/df_agg$sample_weights
  
  out.plot = ggplot(df_agg, aes(x = percentile)) +
    geom_point(aes(y = actual_rate, col = "actual")) + 
    geom_line(aes(y = actual_rate, col = "actual"), size = 1) +
    geom_point(aes(y = pred_rate, col = "preds")) +
    geom_line(aes(y = pred_rate, col = "preds"), size = 1) +
    xlab("Percentile") + ylab("Rate in percentiles") +
    ggtitle(case_when(method == "preds" ~ paste0("Rate by predicted percentile"),
                      method == "actual" ~ paste0("Rate by actual percentile"),
                      TRUE ~ paste0("Not a valid method!")))
  
  
  print(paste0("Correlation, predicted vs. actual rate: ",
               round(cor(df_agg$actual_rate, df_agg$pred_rate, method = "spearman"))))
  
  return(out.plot)
}

# Function: return RMSE performance by percentile
RMSE_ntile = function(preds, actual, n_percentile = 4) {
  # Input: predictions and true values
  # Parameter: n_percentile = number of percentiles to consider (by default, quartiles)
  # Output: RMSE by percentiles
  
  df = data.frame(preds = preds, actual = actual, percentile = ntile(actual, n_percentile))
  df.perf = data.frame(percentile = seq(1:n_percentile), test_RMSE = rep(NA, n_percentile))
  
  for (i in 1:n_percentile) {
    df.pctile = df %>% filter(percentile == i)
    df.perf[i, ]$test_RMSE = RMSE(df.pctile$preds, df.pctile$actual)
  }
  
  return(df.perf)
}

# Function: plot Kolmogorov-Smirmov test & return test statistic
ks_eval = function(sample_1, sample_2) {
  # Input: two samples for KS-test
  # Output: plot of KS test and the KS test statistic
  # Note: can be useful to see the max distributional difference between predicted and actual model
  
  group = c(rep("sample_1"), length(preds),
            rep("sample_2"), length(actual))
  
  df = data.frame(KSD = c(sample_1, sample_2), group = group)
  
  cdf_1 = ecdf(sample_1); cdf_2 = ecdf(sample_2) # create data ecdfs
  
  # Find min amd max to draw line between point of max distance for KS-test
  min_max = seq(min(sample_1, sample_2),
                max(sample_1, sample_2),
                length.out = length(sample_1))
  
  cdf_range = cdf_1(min_max) - cdf_2(min_max)
  
  x_crit = min_max[which(abs(cdf_range) == max(cdf_range))]
  y_1 = cdf_1(x_crit)
  y_2 = cdf_2(x_crit)
  
  # Plot the ecdfs and this difference
  out.plot = ggplot(df, aes(x = KSD, group = group, color = group)) +
    stat_ecdf(size = 1) +
    theme_bw() +
    xlab("Sample range") + ylab("ecdf")  +
    geom_segment(aes(x = x_crit[1], y = y_1[1],
                     xend = x_crit[1], yend = y_2[1]),
                 linetype = "dashed", color = "red") +
    geom_point(aes(x = x_crit[1], y = y_1[1]), color = "red", size = 2) +
    geom_point(aes(x = x_crit[1], y = y_2[1]), color = "red", size = 2) +
    ggtitle("K-S Test Plot")
  
  ks_stat = ks.test(sample_1, sample_2)$statistic
  print(paste0("K-S test statistic: ", round(ks_stat, 5)))
  
  return(out.plot)
}

# Function: plot the distributions of prediction vs. actual response
dist_eval = function(preds, actual, dist_names = c("preds, actual"),
                     x_scale = range(actual), line_size = 0.7, plot_alpha = 0.2, use_hist = FALSE,
                     x_label = "Response variable", plot_title = "Difference in distribution",
                     hist_bins = 100) {
  
  # Input: two vectors of values to compare
  # Output: comparison plot of smoothed distributions
  
  pred.df = data.frame(y_val = preds)
  pred.df$classification = dist_names[1]
  
  actual.df = data.frame(y_val = actual)
  actual.df$classification = dist_names[2]
  
  test.df = rbind(pred.df, actual.df)
  
  if (use_hist == TRUE) {
    ggplot(test.df, aes(x = y_val, fill = classification)) +
      geom_histogram(alpha = plot_alpha, size = line_size, position = "identity", bins = hist_bins) +
      xlab(x_label) + coord_cartesian(x_scale) +
      ggtitle(plot_title)
  } else {
    ggplot(test.df, aes(x = y_val, fill = classification)) +
      geom_density(aes(y = ..density..), alpha = plot_alpha, size = line_size) +
      xlab(x_label) + ylab("Density") + coord_cartesian(x_scale) +
      ggtitle(plot_title)
  }
}

# Variable significance ####

# Function: sort variables by significance from supplied fit
find_sig_vars = function(fit, sig.cutoff = 1) {
  # Input: fit and cutoff for significance (if wanted)
  # Output: variable significance as df
  
  sig.df = data.frame(summary(fit)$coef)
  
  sig.df$var_names = rownames(sig.df)
  rownames(sig.df) = seq(1:nrow(df))
  
  sig.df = sig.df[, c(5, 1, 2, 3, 4)]
  sig.df = format(sig.df, scientific = FALSE)
  
  colnames(sig.df) = c("var_names", "coef", "sd", "z", "p_val")
  sig.df$var_names = as.character(sig.df$var_names)
  sig.df[, -1] = sig.df[, -1] %>% apply(2, as.numeric) %>% apply(2, function(x) round(x, 7))
  
  sig.df = sig.df[sig.df$p_val < sig.cutoff, ]
  return(sig.df)
}


# Function: find nonzero lasso coefficients
find_lasso_coefs = function(lasso.fit, bestlambda) {
  # Input: lasso model and regularization parameter
  # Output: non-zero lasso model coefficients
  
  lasso.coefs = predict(lasso.fit, s = bestlambda, type = "coefficients")
  lasso.coefs = as.matrix(lasso.coefs)
  lasso.coefs = as.data.frame(lasso.coefs)
  lasso.coefs = setDT(lasso.coefs, keep.rownames = TRUE)[]
  names(lasso.coefs) = c("predictor", "coefficient")
  lasso.coefs %>% filter(coefficient != 0) %>% arrange(desc(abs(coefficient)))
}


# VIF: analysis & plots ####

# Function: make sorted df of model VIF
find_VIF = function(fit) {
  # Input: function fit
  # Output: vif as df, sorted in descending order of VIF
  
  vif.df = data.frame(vif(fit))
  colnames(vif.df) = "GVIF"
  vif.df$predictor = rownames(vif.df)
  rownames(vif.df) = seq(1:nrow(vif.df))
  
  vif.df = vif.df %>% select(predictor, GVIF) %>% arrange(desc(GVIF))
  return(vif.df)
}