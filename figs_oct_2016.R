rm(list = ls())

library(rEDM)
library(rgl)

source("model_descriptions.R")
source("analysis_functions.R")
source("figures.R")

make_data <- function(RNGseed = 1234, noise_level = 0.2)
{
    set.seed(1234)
    
    # load and normalize time series
    load("model_2.Rdata")
    mu <- colMeans(data)
    sigma <- apply(data, 2, sd)
    orig_data <- normalize(data)
    noisy_data <- orig_data + rnorm(length(orig_data), sd = sqrt(noise_level))
    
    # load and normalize fine-scale data (for drawing attractors)
    load("model_2_fine.Rdata")
    fine_data <- data - matrix(rep(mu, each = NROW(data)), nrow = NROW(data))
    fine_data <- fine_data / matrix(rep(sigma, each = NROW(data)), nrow = NROW(data))
    
    save(fine_data, orig_data, noisy_data, file = "oct_2016_data.Rdata")
}

make_forecasts <- function(lib = c(551, 575), pred = c(2001, 2100))
{
    load("oct_2016_data.Rdata")
    
    # make block
    n <- NROW(noisy_data)
    block <- cbind(noisy_data[, 1], 
                   c(NA, noisy_data[1:(n-1), 1]), 
                   c(NA, NA, noisy_data[1:(n-2), 1]), 
                   noisy_data[, 2], 
                   c(NA, noisy_data[1:(n-1), 2]), 
                   c(NA, NA, noisy_data[1:(n-2), 2]), 
                   noisy_data[, 3], 
                   c(NA, noisy_data[1:(n-1), 3]), 
                   c(NA, NA, noisy_data[1:(n-2), 3]))
    
    target_cols <- c(1, 4, 7)
    
    # get obs
    obs <- block[(pred[1]+1):pred[2], target_cols]
    
    # get univariate forecasts
    univar_x_pred <- do.call(cbind, lapply(target_cols, function(target_column) {
        return(block_lnlp(block, lib = lib, pred = pred, 
               target_column = target_column, columns = 1:3, 
               short_output = TRUE, stats_only = FALSE)[[1]]$model_output$pred)
    }))
    univar_y_pred <- do.call(cbind, lapply(target_cols, function(target_column) {
        return(block_lnlp(block, lib = lib, pred = pred, 
                          target_column = target_column, columns = 4:6, 
                          short_output = TRUE, stats_only = FALSE)[[1]]$model_output$pred)
    }))
    univar_z_pred <- do.call(cbind, lapply(target_cols, function(target_column) {
        return(block_lnlp(block, lib = lib, pred = pred, 
                          target_column = target_column, columns = 7:9, 
                          short_output = TRUE, stats_only = FALSE)[[1]]$model_output$pred)
    }))
    
    # get multivariate forecasts (using x, y, z)
    multivar_pred <- do.call(cbind, lapply(target_cols, function(target_column) {
        return(block_lnlp(block, lib = lib, pred = pred, 
                          target_column = target_column, columns = target_cols, 
                          short_output = TRUE, stats_only = FALSE)[[1]]$model_output$pred)
    }))

    # do multiview forecasts
    embeddings_list <- t(combn(9, 3, simplify = TRUE))
    valid_embeddings <- apply(embeddings_list %% 3, 1, function (x) {1 %in% x})
    embeddings_list <- embeddings_list[valid_embeddings,]
    mve_pred <- matrix(NA, nrow = pred[2]-pred[1], ncol = 3)
    for(j in 1:3)
    {
        target_column <- target_cols[j]
        in_stats <- block_lnlp(block, lib = lib, pred = lib, num_neighbors = "e+1",
                               target_column = target_column, columns = embeddings_list,
                               silent = TRUE)
        
        in_sample_ranking <- order(in_stats$rho, decreasing = TRUE)
        multivar_embeddings <- embeddings_list[in_sample_ranking[1:8],]
        
        temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = 1,
                           target_column = target_column, columns = multivar_embeddings,
                           stats_only = FALSE, short_output = TRUE)
        
        top_embeddings_preds <- do.call(cbind, lapply(temp, function(x) x$model_output$pred))
        mve_pred[,j] <- rowMeans(top_embeddings_preds)
    }
    
    forecasts <- list(obs = obs, 
                      univar_x_pred = univar_x_pred, 
                      univar_y_pred = univar_y_pred, 
                      univar_z_pred = univar_z_pred, 
                      multivar_pred = multivar_pred, 
                      mve_pred = mve_pred)
    save(forecasts, lib, pred, file = "oct_2016_preds.Rdata")
    return()
}

make_forecast_plot <- function(plot_file = NULL, width = 6, height = 5)
{
    load("oct_2016_data.Rdata")
    load("oct_2016_preds.Rdata")
    fine_t <- (pred[1]*16-15):(pred[2]*16-15)
    x <- (fine_t+15)/16
    
    if(!is.null(plot_file))
        pdf(file = plot_file, width = width, height = height)
    
    #### plot raw time series ####
    par(mfrow = c(3,1), mar = c(1,6,1,1), oma = c(3,0,0,0), mgp = c(4,1,0), 
        las = 1)
    for(j in 1:3)
    {
        plot(x - min(x), fine_data[fine_t, j], type = "l", xaxt = ifelse(j == 3, "s", "n"), 
             xlab = "", ylab = c("x", "y", "z")[j], las = 1, lwd = 2, cex.axis = 1.5, 
             col = c("#FF0000", "#00FF00", "#0000FF")[j])
        if(j != 3)
            axis(side = 1, labels = FALSE)
    }
    # mtext("raw time series", 
    #       side = 3, outer = TRUE)
    
    #### plot attractor ####
    plot3d(fine_data[fine_t,], type = "l", lwd = 2, col = "#888888", 
           xlab = "X", ylab = "Y", zlab = "", axes = FALSE)
    box3d()
    axes3d(edges = c("x--", "y--", "z-+"), cex = 1.5)
    par3d(windowRect = c(0, 70, 600, 670), 
          userMatrix = matrix(c(0.399, -0.917, -0.008, 0, 
                                0.405, 0.168, 0.899, 0, 
                                -0.823, -0.362, 0.438, 0, 
                                0, 0, 0, 1), nrow = 4, byrow = TRUE), 
          scale = c(1.176, 0.898, 0.982))
    observer3d(0, 0, 15.400)
    rgl.postscript("oct_2016_plots.3d.pdf", fmt = "pdf")

    #### for each target ####
    par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0), mgp = c(4,1,0), 
        las = 1)
    my_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FF8800", "#8800FF")
    
    for(j in 1)
    {
        # univariate plots
        obs <- fine_data[fine_t, j]
        preds <- cbind(forecasts$univar_x_pred[,j], 
                              forecasts$univar_y_pred[,j], 
                              forecasts$univar_z_pred[,j], 
                       forecasts$multivar_pred[,j], 
                       forecasts$mve_pred[,j])
        for(i in 1:5)
        {
            x1 <- x - min(x)
            x2 <- seq(from = 1, length.out = NROW(preds))
            plot(x1, obs, type = "l", col = "#888888", 
                 xlab = "Time", ylab = c("X", "Y", "Z")[j], 
                 ylim = range(obs, preds, na.rm = TRUE), 
                 cex.axis = 1.5, lwd = 2)
            lines(x2, preds[,i], col = my_colors[i], lwd = 2)
            rho <- round(cor(obs[match(x2, x1)], preds[,i], use = "pairwise"), 2)
            legend("topleft", lwd = 2, col = my_colors[i], bty = "n", 
                   legend = paste0("rho = ", rho), cex = 1.5)
        }
        
        # 
        # make_panel <- function(x1, y1, x2, y2, color = "#FF0000", 
        #                        x_label = "Time", y_label = "y", ...)
        # {
        #     plot(x1, y1, type = "l", col = "#888888", 
        #          xlab = x_label, ylab = y_label, 
        #          ylim = range(y1, y2, na.rm = TRUE), 
        #          xaxs = "i", ...)
        #     lines(x2, y2, col = color, lwd = 2)
        #     
        #     rho <- round(cor(y1[match(x2, x1)], y2, use = "pairwise"), 2)
        #     legend("topleft", legend = paste0("rho = ", rho))
        # }
        # 
        # par(mfrow = c(5,1), mar = c(1,4,1,1), oma = c(3,0,3,0))    
        # # univariate plots
        # make_panel((fine_t+15)/16, fine_data[fine_t, j], 
        #            (pred[1]+1):pred[2], forecasts$univar_x_pred[,j], 
        #            "#FF0000", y_label = c("x", "y", "z")[j], xaxt = "n")
        # make_panel((fine_t+15)/16, fine_data[fine_t, j], 
        #            (pred[1]+1):pred[2], forecasts$univar_y_pred[,j], 
        #            "#00FF00", y_label = c("x", "y", "z")[j], xaxt = "n")
        # make_panel((fine_t+15)/16, fine_data[fine_t, j], 
        #            (pred[1]+1):pred[2], forecasts$univar_z_pred[,j], 
        #            "#0000FF", y_label = c("x", "y", "z")[j], xaxt = "n")
        # 
        # # multivariate_plot
        # make_panel((fine_t+15)/16, fine_data[fine_t, j], 
        #            (pred[1]+1):pred[2], forecasts$multivar_pred[,j], 
        #            "#FF00FF", y_label = c("x", "y", "z")[j], xaxt = "n")
        # 
        # # multiview plot
        # make_panel((fine_t+15)/16, fine_data[fine_t, j], 
        #            (pred[1]+1):pred[2], forecasts$mve_pred[,j], 
        #            "#FF8800", y_label = c("x", "y", "z")[j])
        # 
        # mtext(paste0("predicting ", c("x", "y", "z")[j]), 
        #       side = 3, outer = TRUE)
    }
    if(!is.null(plot_file))
        dev.off()
    
    return()
}

if(FALSE)
{
    make_data(noise_level = 0.25)
    make_forecasts(lib = c(556, 580), pred = c(2011, 2110))
    make_forecast_plot()
    make_forecast_plot("oct_2016_plots.pdf", width = 8)
}

# plot pred no noise
# fine_pred <- (pred[1]*16-15):(pred[2]*16-15)
# plot3d(fine_data[fine_pred,], type = "l", col = "#888888")s