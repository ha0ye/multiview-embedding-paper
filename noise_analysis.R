rm(list = ls())

library(rEDM)
library(plyr)
library(ggplot2)
library(parallel)
library(tidyr)

source("model_descriptions.R")
source("entropy_calculations.R")
source("analysis_functions.R")
source("figures.R")

if(FALSE) # generate data
{
    obs_err_cv <- 0.2
    make_model_1(obs_err_cv = obs_err_cv)
    make_model_2(obs_err_cv = obs_err_cv)
    make_HW_model(obs_err_cv = obs_err_cv)
    downsample_HW_model()

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
        run_multiembed("model_hw_downsampled.Rdata", paste("multiembed_results_hw_", i, ".Rdata", sep = ""), 
                       lib_start = lib_locations[i], num_to_merge = 21)
    }, mc.cores = 8)
    
    combine_multiembed_results(1, my_path = "")
    combine_multiembed_results(2, my_path = "")
    combine_multiembed_results("hw", my_path = "")
    
    plot_full_model_results("multiembed_all_results_1.Rdata", out_file = "ed_figure_1.pdf")
    plot_full_model_results("multiembed_all_results_2.Rdata", out_file = "ed_figure_2.pdf")
    plot_full_model_results("multiembed_all_results_hw.Rdata", out_file = "ed_figure_hw.pdf", 
                            width = 10.5)
    
    plot_combined_model_results(out_file = "figure_3.pdf")
    
    
    
}
