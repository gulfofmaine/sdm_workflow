#####
## Taylor Diagram Function
#####

# Helper functions: correlation coefficient and bias ----------------------
corrcoeff_func_simp<- function(df){
  df_use<- df %>%
    drop_na(obs, mod)
  mean_obs<- mean(df_use$obs)
  mean_mod<- mean(df_use$mod)
  sd_obs<- sd(df_use$obs)
  sd_mod<- sd(df_use$mod)
  samps<- nrow(df_use)
  
  out<- ((1/samps)*(sum((df_use$mod - mean_mod)*(df_use$obs - mean_obs))))/(sd_obs*sd_mod)
  return(out)
}

bias_func_simp<- function(df){
  df_use<- df %>%
    drop_na(obs, mod)
  out<- sd(df_use$mod)/sd(df_use$obs)
  return(out)
}


# Main Taylor Diagram function --------------------------------------------
taylor_diagram<- function(dat, obs = "obs", mod = "mod", group = NULL, out_dir = NULL, grad_corr_lines = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1), pcex = 1, cex_axis = 1, normalize = TRUE, mar = c(5, 4, 6, 6), sd_r = 1, pt_col = NULL, pt_cols = NULL, shapes = NULL, example = FALSE) {
  
  ## Details
  # This function plots a Taylor Diagram of model prediction accuracy, sumamrizing the root mean square error, the coefficient of determination, and the ratio of standard deviations. 
  
  # Args:
  # dat = data frame with observations and model predictions, as well as group if necessary
  # obs = Column name for the observation response
  # mod = Column name for the modeled response
  # group = Grouping variable, used for comparing different species/ages or stages/models, etc
  # out.dir = Directory where to save the Taylor Diagram plot
  # ... All these other things correspond to some of the aesthetics of the plot. pt.col gives color if just plotting one point (group, model), pt.cols is a vector of colors for plotting multiple points (groups, models) on one plot.
  
  # Returns: NULL; saves plot to output directory
  
  ## Start function
  # Install libraries
  library(tidyverse)
  
  # Set arguments for debugging -- this will NOT run when you call the function. Though, you can run each line inside the {} and then you will have everything you need to walk through the rest of the function.
  if(FALSE){
    dat = mod_res$Pred_Baseline_Biomass[[1]]
    obs = "SUM_BIOMASS"
    mod = "Pred"
    group<- NULL
    out_dir<- NULL
    grad_corr_lines = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
    pcex = 1
    cex_axis = 1
    normalize = TRUE
    mar = c(5, 4, 6, 6)
    sd_r = 1
    pt_col<- "#006d2c"
  }
  if(example){
    # Create a data set with observations and predictions
    data(trees)
    tree_mod<- lm(Volume ~ Girth, data = trees)
    trees$Volume_pred<- as.numeric(predict(tree_mod))
    dat<- trees
    obs<- "Volume"
    mod<- "Volume_pred"
    group<- NULL
    out_dir<- "~/Desktop/"
    grad_corr_lines = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
    pcex = 1
    cex_axis = 1
    normalize = TRUE
    mar = c(5, 4, 6, 6)
    sd_r = 1
    pt_col<- "#006d2c"
  }
  
  # Some house keeping -- rename the obs and mod columns to work with generic functions
  old_names<- c(obs, mod)
  new_names<- c("obs", "mod")
  dat<- dat %>%
    rename_at(vars(all_of(old_names)), ~ new_names)
  
  # Calculate the correlation coefficient and bias, two stats needed for Taylor Diagram. Flexibility to group and then calculate stats by group (e.g., species, model, etc).
  if(is.null(group)){
    mod_stats<- dat %>%
      nest(data = everything()) %>%
      mutate(., "CorrCoeff" = as.numeric(map(data, corrcoeff_func_simp)),
             "Bias" = as.numeric(map(data, bias_func_simp)))
  } else {
    # Group by group and calculate stats
    mod_stats<- dat %>%
      group_by_at(group) %>%
      nest() %>%
      mutate(., "CorrCoeff" = as.numeric(map(data, corrcoeff_func_simp)),
             "Bias" = as.numeric(map(data, bias_func_simp)))
  }
  
  # Now plot creation....
  # Getting maxSD for plotting
  maxsd<- max(mod_stats$Bias, 1)
  
  # Empty plot first
  # Creating empty plot first
  plot_base<- ggplot() + 
    scale_x_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd+0.01), breaks = seq(from = 0, to = maxsd, by = 0.5), expand = c(0, 0)) +
    scale_y_continuous(name = "Standard deviation (normalized)", limits = c(-0.015, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5), expand = c(0, 0)) +
    theme_classic()
  
  # Coeff D rays 
  for(i in 1:length(grad_corr_lines)){
    x_vec<- c(0, maxsd*grad_corr_lines[i])
    y_vec<- c(0, maxsd*sqrt(1 - grad_corr_lines[i]^2))
    
    if(i ==1){
      coeffd_rays_df<- data.frame("Ray" = rep(1, length(x_vec)), "x" = x_vec, "y" = y_vec)
    } else {
      temp<- data.frame("Ray" = rep(i, length(x_vec)), "x" = x_vec, "y" = y_vec)
      coeffd_rays_df<- bind_rows(coeffd_rays_df, temp)
    }
  }
  
  # Add rays
  plot_coeffd<- plot_base +
    geom_line(data = coeffd_rays_df, aes(x = x, y = y, group = Ray), lty = "longdash", col = "lightgray") 
  
  coeffd_labs<- coeffd_rays_df %>%
    group_by(Ray) %>%
    summarize(., 
              "x" = max(x, na.rm = TRUE), 
              "y" = max(y, na.rm = TRUE)) %>%
    data.frame()
  
  coeffd_labs$Label<- grad_corr_lines
  
  plot_coeffd<- plot_coeffd +
    geom_label(data = coeffd_labs, aes(x = x, y = y, label = Label), fill = "white", label.size = NA)
  
  # SD arcs
  # Need to add in SD arcs
  sd_arcs<- seq(from = 0, to = maxsd, by = 0.5)
  
  for(i in 1:length(sd_arcs)){
    x_vec<- sd_arcs[i]*cos(seq(0, pi/2, by = 0.03))
    y_vec<- sd_arcs[i]*sin(seq(0, pi/2, by = 0.03))
    
    if(i ==1){
      sd_arcs_df<- data.frame("Arc" = rep(sd_arcs[1], length(x_vec)), "x" = x_vec, "y" = y_vec)
    } else {
      temp<- data.frame("Arc" = rep(sd_arcs[i], length(x_vec)), "x" = x_vec, "y" = y_vec)
      sd_arcs_df<- bind_rows(sd_arcs_df, temp)
    }
  }
  
  # Add arcs to plot.base
  plot_sd<- plot_coeffd +
    geom_line(data = sd_arcs_df, aes(x = x, y = y, group = Arc), lty = "dotted", color = "lightgray") 
  
  # Now gamma? -- Standard deviation arcs around the reference point
  gamma<- pretty(c(0, maxsd), n = 4)[-1]
  gamma<- gamma[-length(gamma)]
  labelpos<- seq(45, 70, length.out = length(gamma))
  
  for(gindex in 1:length(gamma)) {
    xcurve<- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + sd_r
    endcurve<- which(xcurve < 0)
    endcurve<- ifelse(length(endcurve), min(endcurve) - 1, 105)
    ycurve<- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
    maxcurve<- xcurve * xcurve + ycurve * ycurve
    startcurve<- which(maxcurve > maxsd * maxsd)
    startcurve<- ifelse(length(startcurve), max(startcurve) + 1, 0)
    x_vec<- xcurve[startcurve:endcurve]
    y_vec<- ycurve[startcurve:endcurve]
    
    if(gindex ==1){
      gamma_df<- data.frame("Gamma" = rep(gamma[1], length(x_vec)), "x" = x_vec, "y" = y_vec)
    } else {
      temp<- data.frame("Gamma" = rep(gamma[gindex], length(x_vec)), "x" = x_vec, "y" = y_vec)
      gamma_df<- bind_rows(gamma_df, temp)
    }
  }
  
  gamma_df$Gamma<- factor(gamma_df$Gamma, levels = unique(gamma_df$Gamma))
  
  # Add em
  plot_gamma<- plot_sd +
    geom_line(data = gamma_df, aes(x = x, y = y, group = Gamma), lty = "solid", col = "lightgray")
  
  # Label...
  gamma_labs<- gamma_df %>%
    group_by(Gamma) %>%
    summarize("x" = mean(x, na.rm = TRUE), 
              "y" = median(y, na.rm = TRUE))
  
  inflection_func<- function(df){
    d1<- diff(df$y)/diff(df$x)    
    pt_id<- which.max(d1)
    pt_out<- df[pt_id,]
    pt_out$y<- rep(0, nrow(pt_out))
    return(pt_out)
  }
  
  gamma_labs<- gamma_df %>%
    group_by(Gamma) %>%
    nest() %>%
    summarize("pt" = map(data, inflection_func)) %>%
    unnest(cols = c(pt))
  
  plot_gamma<- plot_gamma +
    geom_label(data = gamma_labs, aes(x = x, y = y, label = Gamma), fill = "white", label.size = NA)
  
  # Add in reference point
  plot_all<- plot_gamma +
    geom_point(aes(x = sd_r, y = 0), color = "black", size = 4) 
  
  # Add in reference points
  mod_td<- mod_stats %>%
    mutate(., "TD_X" = Bias * CorrCoeff,
           "TD_Y" = Bias * sin(acos(CorrCoeff)))
  
  if(is.null(group)){
    plot_td<- plot_all +
      geom_point(data = mod_td, aes(x = TD_X, y = TD_Y), color = pt_col, size = 3.5) +
      geom_text(aes(label = "Correlation coefficient", x = 0.8, y = 0.75), angle = -38)
  } else {
    plot_td<- plot_all +
      geom_point(data = mod_td, aes(x = TD_X, y = TD_Y, color = pt_cols), size = 3.5) +
      scale_color_manual(name = "Group", values = pt_cols) +
      scale_shape_manual(name = "Group", values = shapes) +
      geom_text(aes(label = "Correlation coefficient", x = 0.8, y = 0.75), angle = -38)
  }
  if(is.null(out_dir)){
    return(plot_td)
  } else {
    ggsave(paste(out_dir, "_TaylorDiagram.jpg", sep = ""), plot_td)
    return(plot_td)
  }
}