rm(list = ls())

library(rEDM)
library(plyr)
library(ggplot2)
library(parallel)
library(tidyr)
library(rgl)

source("model_descriptions.R")
source("entropy_calculations.R")
source("analysis_functions.R")
source("mesocosm_functions.R")
source("figures.R")
# source("lorenz_demo.R")
#source("make_figures.R")
#source("supplemental_analysis.R")

if(FALSE) # generate data
{
    make_model_1()
    make_model_2()
    make_model_3()
}

if(FALSE) # compare univariate forecasts with multiembed forecasts (Figure 2)
{
    lib_start <- 501
    lib_size <- 25
    compute_forecast_errors(lib = c(lib_start, lib_start + lib_size - 1))
    make_error_plot(to_plot = c(3001, 4000),
                    paste("figure_2_panels_", lib_start, "_to_", lib_start + lib_size - 1, ".pdf", sep = ""))
}

if(FALSE) # test convergence of entropy measures
{
    x <- rnorm(1000)
    cat("Estimated entropy by histogram was: ", entropy_hist(x), "\n", sep = "")
    cat("Estimated entropy by kernel was: ", entropy_kernel(x), "\n", sep = "")
    
    compute_entropy_convergence("model_1.Rdata", "entropy_results_1.Rdata")
    compute_entropy_convergence("model_2.Rdata", "entropy_results_2.Rdata")
    compute_entropy_convergence("model_3.Rdata", "entropy_results_3.Rdata")
    
    plot_entropy_results("entropy_calculations.pdf")
}

if(FALSE) # multiembed analysis
{
    num_trials <- 100
    set.seed(1234)
    lib_locations <- floor(runif(100, min = 501, max = 2001))
    mclapply(1:num_trials, function(i) {
        run_multiembed("model_1.Rdata", paste("multiembed_results_1_", i, ".Rdata", sep = ""), 
                       lib_start = lib_locations[i], num_to_merge = 8)
    }, mc.cores = 8)
    
    mclapply(1:num_trials, function(i) {
        run_multiembed("model_2.Rdata", paste("multiembed_results_2_", i, ".Rdata", sep = ""), 
                       lib_start = lib_locations[i], num_to_merge = 8)
    }, mc.cores = 8)
    
    mclapply(1:num_trials, function(i) {
        run_multiembed("model_3.Rdata", paste("multiembed_results_3_", i, ".Rdata", sep = ""), 
                       lib_start = lib_locations[i], num_to_merge = 8)
    }, mc.cores = 8)
    
    combine_multiembed_results(1)
    combine_multiembed_results(2)
    combine_multiembed_results(3)
    
    plot_full_model_results("multiembed_all_results_1.Rdata", out_file = "ed_figure_1.pdf")
    plot_full_model_results("multiembed_all_results_2.Rdata", out_file = "ed_figure_2.pdf")
    plot_full_model_results("multiembed_all_results_3.Rdata", out_file = "ed_figure_3.pdf")
    
    plot_combined_model_results(out_file = "figure_3.pdf")
}

if(FALSE) # mesocosm stuff
{
    process_mesocosm_data()
    
    mesocosm_ccm()
    plot_mesocosm_ccm(out_file = "figure_4b.pdf")
    
    mesocosm_multiembed()
    mesocosm_multiembed_stats()
    plot_mesocosm_multiembed("figure_4c.pdf", var = "rho")
    plot_mesocosm_multiembed_full("ed_figure_4.pdf")
}

if(FALSE) # see how forecast skill changes with number of embeddings to merge
{
    num_trials <- 100
    set.seed(1234)
    lib_locations <- floor(runif(100, min = 501, max = 2001))
    mclapply(1:num_trials, function(i) {
        run_multiembed_preds("model_1.Rdata", lib_start = lib_locations[i], 
                             paste("multiembed_preds_1_", i, ".Rdata", sep = ""))}, 
        mc.cores = 4)
    
    mclapply(1:num_trials, function(i) {
        run_multiembed_preds("model_2.Rdata", lib_start = lib_locations[i], 
                             paste("multiembed_preds_2_", i, ".Rdata", sep = ""))}, 
        mc.cores = 4)
    
    mclapply(1:num_trials, function(i) {
        run_multiembed_preds("model_3.Rdata", lib_start = lib_locations[i], 
                             paste("multiembed_preds_3_", i, ".Rdata", sep = ""))}, 
        mc.cores = 4)
    mclapply(1:num_trials, function(i) {
        run_multiembed_preds("model_hw_downsampled.Rdata", lib_start = lib_locations[i], 
                             paste("multiembed_pred_files/multiembed_preds_hw_", i, ".Rdata", sep = ""))}, 
        mc.cores = 8)
    
    compute_pred_skill_and_merge()
    compute_pred_skill_and_merge(base_pred_file = "multiembed_preds_2_", 
                                 output_file = "merged_stats_2.Rdata")
    compute_pred_skill_and_merge(base_pred_file = "multiembed_preds_3_", 
                                 output_file = "merged_stats_3.Rdata")
    compute_pred_skill_and_merge(base_pred_file = "multiembed_preds_hw_", 
                                 output_file = "merged_stats_hw.Rdata")
    
    plot_multiembed_perf_vs_nm(plot_file = "multiembed_rho_vs_nm_1.pdf")
    plot_multiembed_perf_vs_nm(stats_file = "merged_stats_2.Rdata", 
                               plot_file = "multiembed_rho_vs_nm_2.pdf")
    plot_multiembed_perf_vs_nm(stats_file = "merged_stats_3.Rdata", 
                               plot_file = "multiembed_rho_vs_nm_3.pdf")
    plot_multiembed_perf_vs_nm(stats_file = "merged_stats_hw.Rdata", 
                               plot_file = "multiembed_rho_vs_nm_hw.pdf", width = 10)
}

if(FALSE) # huisman-weissing model analysis
{
    make_HW_model()
    downsample_HW_model()
    make_HW_univariate_plots()
    
    # run multiembed trials
    num_trials <- 100
    set.seed(1234)
    lib_locations <- floor(runif(100, min = 501, max = 2001))
    mclapply(1:num_trials, function(i) {
        run_multiembed("model_hw_downsampled.Rdata", paste("multiembed_results_hw_", i, ".Rdata", sep = ""), 
                       lib_start = lib_locations[i], num_to_merge = 21)
    }, mc.cores = 8)
    
    combine_multiembed_results("hw")

    plot_full_model_results("multiembed_all_results_hw.Rdata", out_file = "ed_figure_hw.pdf", 
                            width = 10.5)

    plot_combined_model_results(out_file = "figure_3.pdf")    
}


if(FALSE) # summary statistics
{
    output_some_numbers()
    load("merged_multiembed_stats.Rdata")
    
    cat("=== ALL LIBS ===\n")
    
    my_range <- range(aggregated_stats$delta_univar_rho)
    cat("range of delta_univar_rho: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    my_range <- range(aggregated_stats$delta_multivar_rho)
    cat("range of delta_multivar_rho: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    
    my_range <- range(aggregated_stats$delta_univar_mae)
    cat("range of delta_univar_mae: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    my_range <- range(aggregated_stats$delta_multivar_mae)
    cat("range of delta_multivar_mae: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    
    my_range <- range(aggregated_stats$delta_univar_rmse)
    cat("range of delta_univar_rmse: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    my_range <- range(aggregated_stats$delta_multivar_rmse)
    cat("range of delta_multivar_rmse: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    
    aggregated_stats <- aggregated_stats[aggregated_stats$lib == 25,]
    cat("=== LIB = 25 ===\n")
    
    my_range <- range(aggregated_stats$delta_univar_rho)
    cat("range of delta_univar_rho: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    my_range <- range(aggregated_stats$delta_multivar_rho)
    cat("range of delta_multivar_rho: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    
    my_range <- range(aggregated_stats$delta_univar_mae)
    cat("range of delta_univar_mae: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    my_range <- range(aggregated_stats$delta_multivar_mae)
    cat("range of delta_multivar_mae: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    
    my_range <- range(aggregated_stats$delta_univar_rmse)
    cat("range of delta_univar_rmse: ", my_range[1], " - ", my_range[2], "\n", sep = "")
    my_range <- range(aggregated_stats$delta_multivar_rmse)
    cat("range of delta_multivar_rmse: ", my_range[1], " - ", my_range[2], "\n", sep = "")
}

if(FALSE) # make figures
{
    # raw time series
    make_figure_0(file = "model_1.Rdata", output = "figure_0a.pdf")
    make_figure_0(file = "model_2.Rdata", output = "figure_0b.pdf")
    make_figure_0(file = "model_3.Rdata", output = "figure_0c.pdf")
    
    # make_figures_univar() # figure S1 (univariate E vs. Rho)
    make_figure_1("model_1.Rdata") # cross-ts plots
    make_figure_2("model_1.Rdata") # cross-ts plots with lags  
    make_figure_3("multiembed_results_1.Rdata", lib = 25) # in-sample vs. out-of-sample forecast skill
    make_figure_4("multiembed_results_1.Rdata", lib = 25) # total time lag vs. out-of-sample forecast skill
    
    # comparison across prediction methods & library size
    make_figure_5(file = "multiembed_all_results_1.Rdata", 
                  title = "Model 1 (obs_err = 0.1, proc_err = 0.05)", output = "figure_5a.pdf")
    make_figure_5(file = "multiembed_all_results_2.Rdata", 
                  title = "Model 2 (obs_err = 0.1, proc_err = 0)", output = "figure_5b.pdf")
    make_figure_5(file = "multiembed_all_results_3.Rdata", 
                  title = "Model 3 (obs_err = 0.1, proc_err = 0.05)", output = "figure_5c.pdf")
    
    make_figure_6("multiembed_results_1.Rdata", lib = 25)
    
    # comparison across prediction methods & library size
    make_figure_8(output = "figure_8.pdf")
    
    # comparison w/ concatenation method data available:
    make_figure_9(outout = "figure_9.pdf")
}

# if(FALSE) # univariate analysis
# {
#     start_time <- proc.time()
#     message("Starting univariate calculations", appendLF = FALSE)
#     compute_univariate("model_1.Rdata", "univariate_results_1.Rdata", lib_sizes)
#     compute_univariate("model_2.Rdata", "univariate_results_2.Rdata", lib_sizes)
#     compute_univariate("model_3.Rdata", "univariate_results_3.Rdata", lib_sizes)
#     message("done! (", round((proc.time()-start_time)[3], 3), " seconds elapsed.)")
# }

# if(FALSE) # supplemental analysis
# {
#     run_in_out_analysis()
#     make_in_out_plots()
#     make_time_lag_rho_plots()
# }