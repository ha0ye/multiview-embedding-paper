rm(list = ls())

library(rEDM)
library(parallel)
library(dplyr)
library(ggplot2)
#library(plyr)
#library(tidyr)
library(rgl)

source("model_descriptions.R")
source("analysis_functions.R")
source("figures.R")
source("science_rev_analysis.R")
source("mesocosm_functions.R")
#source("supplemental_analysis.R")

obs_err <- c(0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
systems <- c("1", "2", "3")

if(FALSE) # generate data
{
    make_model_1()
    make_model_2()
    make_model_3()
    make_model_4()
    make_HW_model()
}

if(FALSE) # run using different values of observation error and for different systems
{
    obs_err <- c(0.0, 0.10, 0.20, 0.30)
    num_trials <- 100
    set.seed(1234)
    lib_locations <- floor(runif(100, min = 501, max = 2001))
    
    for(obs_err_var in obs_err)
    {
        for(var in "4")
        {
            data_file <- paste0("model_", var, ".Rdata")
            pred_file <- paste0("output-nosync/model_", var, "_output_", obs_err_var, ".Rdata")
            stat_file <- paste0("model_", var, "_stats_", obs_err_var, ".Rdata")
            run_multiembed(data_file = data_file, 
                           out_file = pred_file, 
                           lib_start = lib_locations[61:100], 
                           obs_err_var = obs_err_var, 
                           lib_sizes = c(25, 50, 75, 100))
            compute_multiembed_stats(in_file = pred_file, 
                                     out_file = stat_file)
        }
    }
    
    merge_stats_files(obs_err)
}

if(FALSE) # compute stats using different values of k
{
    obs_err_var <- 0.1
    for(var in systems)
    {
        pred_file <- paste0("model_", var, "_output_", obs_err_var, ".Rdata")
        stat_file <- paste0("model_", var, "_stats_", obs_err_var, "_k.Rdata")
        compute_multiembed_stats_k(in_file = pred_file, 
                                   out_file = stat_file)
    }
    
    merge_stats_files_k()
}

if(FALSE) # plot rho vs. lib size
{
    my_plots <- lapply(obs_err, plot_metric_vs_lib_size)
    
    pdf("mve_plot_rho_vs_lib_size.pdf", width = 10, height = 8)
    print(my_plots)
    dev.off()
}

if(FALSE) # plot performance vs. lib size
{
    my_plots <- lapply(systems, plot_perf_vs_lib_size)
    
    pdf("mve_plot_perf_vs_lib_size.pdf", width = 7, height = 8)
    print(my_plots)
    dev.off()
}

if(FALSE) # plot performance vs. noise
{
    my_plots <- lapply(systems, plot_metric_vs_noise)
    
    pdf("mve_plot_mae_vs_noise.pdf", width = 8, height = 8)
    print(my_plots)
    dev.off()
    
    my_plots <- lapply(systems, plot_metric_vs_noise, "rmse")
    
    pdf("mve_plot_rmse_vs_noise.pdf", width = 8, height = 8)
    print(my_plots)
    dev.off()
    
    my_plots <- lapply(systems, plot_metric_vs_noise, "rho")
    
    pdf("mve_plot_rho_vs_noise.pdf", width = 8, height = 8)
    print(my_plots)
    dev.off()
}

if(FALSE) # plot rho vs. k
{
    my_plots <- lapply(systems, plot_perf_vs_k)
    
    pdf("mve_plot_rho_vs_k.pdf", width = 10, height = 8)
    print(my_plots)
    dev.off()
}

if(FALSE) # show embedding-unique errors
{
    compute_model_2_forecast_errors(lib = c(501, 525))
    make_error_plot(to_plot = c(2001, 3000))
    make_error_attractor_ts_plot(lib = c(501, 525))
}

if(FALSE) # show sparseness and noise
{
    # make more finely sampled data for model 2 for prettiness
    make_model_2("model_2_fine.Rdata", downsample = 50)
    
    plot_sparseness_model_2(pred_file = "figure_1a.pdf")
    
    plot_noise_model_2(pred_file = "figure_1b.pdf")
    
    plot_multiview_model_2()
}

if(FALSE) # mesocosm stuff
{
    process_mesocosm_data(fourth_root = FALSE)
    
    mesocosm_ccm()
    plot_mesocosm_ccm(out_file = "figure_4b.pdf")
    
    mesocosm_multiembed(ranking = "mae")
    mesocosm_multiembed_stats()
    # plot_mesocosm_multiembed("figure_4c.pdf")
    # plot_mesocosm_multiembed_full("ed_figure_4.pdf")
}

if(FALSE) # make combined figure 1
{
    pred_target <- 438
    lib <- c(501, 525)
    spline_t <- c(401, 525)
    native_view <- matrix(c(0.406, -0.914, -0.008, 0,
                            0.485, 0.208, 0.849, 0,
                            -0.775, -0.349, 0.528, 0,
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    native_view <- matrix(c(0.514, -0.857, -0.034, 0,
                            0.367, 0.184, 0.912, 0,
                            -0.775, -0.481, 0.409, 0,
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    native_view <- matrix(c(0.611, -0.791, -0.023, 0,
                            0.194, 0.121, 0.974, 0,
                            -0.767, -0.600, 0.227, 0,
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    y_view <- matrix(c(0.939, -0.029, -0.343, 0,
                       -0.053, 0.972, -0.227, 0,
                       0.340, 0.231, 0.912, 0,
                       0, 0, 0, 1), nrow = 4, byrow = TRUE)
    z_view <- matrix(c(0.942, 0.3333, -0.033, 0,
                       -0.054, 0.246, 0.968, 0,
                       0.330, -0.910, 0.250, 0,
                       0, 0, 0, 1), nrow = 4, byrow = TRUE)
    trans_gray <- rgb(204, 204, 204, alpha = 128, maxColorValue = 256)
    
    # panel A (sparse neighbors when lib_size = 25)
    plot_sparseness_model_2(pred_target = pred_target, 
                            pred = spline_t, lib = lib, 
                            user_matrix = native_view, 
                            plot_file = "figure_rev_1a.pdf")
    
    # panel B (time series and some attractors)
    plot_raw_time_series(lib = spline_t, 
                         plot_file = "figure_rev_1bi.pdf")
    plot_attractor_model_2(spline_t = spline_t, spline_col = "red", 
                           embedding = c(1, 0, 1, 1, 1, 2),
                           plot_file = "figure_rev_1bii.pdf")
    plot_attractor_model_2(spline_t = spline_t, spline_col = "green",
                           user_matrix = y_view,
                           embedding = c(2, 0, 2, 1, 2, 2),
                           plot_file = "figure_rev_1biii.pdf")
    plot_attractor_model_2(spline_t = spline_t, spline_col = "blue",
                           user_matrix = z_view,
                           embedding = c(3, 0, 3, 1, 3, 2),
                           plot_file = "figure_rev_1biiii.pdf")
    
    # panel C (univariate attractors)
    plot_attractor_model_2(spline_t = lib, spline_col = trans_gray, 
                           point_t = lib, point_col = "red", 
                           embedding = c(1, 0, 1, 1, 1, 2),
                           plot_file = "figure_rev_1ci.pdf")
    plot_attractor_model_2(spline_t = lib, spline_col = trans_gray,
                           point_t = lib, point_col = "green", 
                           user_matrix = y_view,
                           embedding = c(2, 0, 2, 1, 2, 2),
                           plot_file = "figure_rev_1cii.pdf")
    plot_attractor_model_2(spline_t = lib, spline_col = trans_gray,
                           point_t = lib, point_col = "blue", 
                           user_matrix = z_view,
                           embedding = c(3, 0, 3, 1, 3, 2),
                           plot_file = "figure_rev_1ciii.pdf")
    
    # oanel d (forecasts from univariate attractors)
    load("pred_errors.Rdata")
    plot_attractor_and_forecasts_model_2(mx_forecasts[2001:2999,], point_col = "red", 
                                         spline_t = spline_t, spline_col = trans_gray, 
                                         user_matrix = native_view, 
                                         plot_file = "figure_rev_1di.pdf")
    plot_attractor_and_forecasts_model_2(my_forecasts[2001:2999,], point_col = "green", 
                                         spline_t = spline_t, spline_col = trans_gray, 
                                         user_matrix = native_view, 
                                         plot_file = "figure_rev_1dii.pdf")
    plot_attractor_and_forecasts_model_2(mz_forecasts[2001:2999,], point_col = "blue", 
                                         spline_t = spline_t, spline_col = trans_gray, 
                                         user_matrix = native_view, 
                                         plot_file = "figure_rev_1diii.pdf")
    
    # panel E (multiview forecasts)
    ### find top 8 embeddings
    best_embeddings <- unique(do.call(rbind, best_embeddings))
    for(i in 1:NROW(best_embeddings))
    {
        plot_file <- paste0("figure_rev_1e", strrep("i", i), ".pdf")
        curr_embedding <- best_embeddings[i,]
        embedding <- c((curr_embedding[1]+2) %/% 3, (curr_embedding[1]-1) %% 3, 
                       (curr_embedding[2]+2) %/% 3, (curr_embedding[2]-1) %% 3, 
                       (curr_embedding[3]+2) %/% 3, (curr_embedding[3]-1) %% 3)
        plot_attractor_model_2(spline_t = lib, spline_col = trans_gray,
                                             point_t = lib, point_col = "black", 
                                             user_matrix = native_view,
                                             embedding = embedding,
                                             plot_file = plot_file)
    }
    
    ### plot MVE forecasts
    plot_attractor_and_forecasts_model_2(multi_forecasts[2001:2999,], point_col = "purple", 
                                         spline_t = spline_t, spline_col = trans_gray, 
                                         user_matrix = native_view, 
                                         plot_file = "figure_rev_1e.pdf")
    
}

if(FALSE) # make new panel for figure 2a
{
    pred_target <- 438
    lib <- c(501, 525)
    spline_t <- c(401, 525)
    native_view <- matrix(c(0.406, -0.914, -0.008, 0,
                            0.485, 0.208, 0.849, 0,
                            -0.775, -0.349, 0.528, 0,
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    
    plot_sparseness_model_2(pred_target = pred_target, pred = spline_t, 
                            lib = c(501, 600), 
                            user_matrix = native_view, 
                            make_prediction = TRUE, 
                            plot_file = "ED_Figure_0a.pdf")
}
