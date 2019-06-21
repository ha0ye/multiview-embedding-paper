process_mesocosm_data <- function(fourth_root = TRUE)
{
    data <- read.table("mesocosm_data.txt", header = TRUE)
    lib <- as.matrix(read.table("mesocosm_data.lib"))
    
    # check for long stretches of 0s
    names(data) <- c("day", "cyclopoids", "calanoid_copepods", "rotifers", "protozoa", 
                     "nanophytoplankton", "picophytoplankton", "filamentous_diatoms", 
                     "ostracods", "harpacticoids", "bacteria", "nitrite", "nitrate", 
                     "ammonium", "nitrogen", "phosphorus")
    
    obs_data <- data
    if(fourth_root)
        obs_data[,2:NCOL(obs_data)] <- sqrt(sqrt(obs_data[,2:NCOL(obs_data)]))
    normed_data <- normalize(obs_data)
    normed_data[,1] <- data[,1]
    
    save(obs_data, normed_data, lib, file = "mesocosm_data.Rdata")
    return()
}

mesocosm_ccm <- function(lib_sizes = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650))
{
    load("mesocosm_data.Rdata")
    data <- normed_data
    
    cross_map_var_list <- data.frame(from = c("calanoid_copepods", "calanoid_copepods", "rotifers", "rotifers"), 
                                     to = c("nanophytoplankton", "picophytoplankton", "nanophytoplankton", "picophytoplankton"))
    
    ccm_results <- data.frame()
    start_time <- proc.time()
    message("Computing CCM for mesocosm data", appendLF = FALSE)
    for(i in 1:NROW(cross_map_var_list))
    {
        message(".", appendLF = FALSE)
        temp <- ccm(data, lib = lib, pred = lib, 
                    lib_sizes = lib_sizes, E = 3, random_libs = FALSE, 
                    lib_column = cross_map_var_list$from[i], 
                    target_column = cross_map_var_list$to[i], silent = TRUE)
        temp$lib_column <- cross_map_var_list$from[i]
        temp$target_column <- cross_map_var_list$to[i]
        ccm_results <- rbind(ccm_results, temp)
    }
    message("done! (", round((proc.time()-start_time)[3], 3), " seconds elapsed.)")
    save(ccm_results, file = "mesocosm_ccm_results.Rdata")
    return()
}

plot_mesocosm_ccm <- function(out_file = NULL)
{
    load("mesocosm_ccm_results.Rdata")
    
    to_plot <- summarySE(ccm_results, measurevar = "rho", 
                         groupvars = c("lib_size", "lib_column", "target_column"))
    my_plot <- ggplot(to_plot, aes(x = lib_size, y = rho, color = lib_column)) + 
        geom_errorbar(aes(ymin = rho - sd, ymax = rho + sd)) + 
        geom_line() + 
        facet_wrap(~ target_column, ncol = 1) + 
        scale_color_manual(values = c("green4", "darkorchid2")) + 
        xlab("library size (l)") + 
        ylab("cross map skill (rho)") + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              legend.background = element_rect(color = "black"), 
              legend.key = element_blank(), 
              legend.title = element_blank(), 
              legend.position = c(0, 1), 
              legend.justification = c(0, 1), 
              panel.background = element_rect(color = "black", fill = NA))
    print(my_plot)
    
    if(!is.null(out_file))
        pdf(file = out_file, width = 6, height = 4)
    
    print(my_plot)
    
    if(!is.null(out_file))
        dev.off()
    return()
}

mesocosm_multiembed <- function(...)
{
    run_multiembed <- function(data, lib, E = 3, target = 1, max_lag = 3, num_to_merge = 8, ...)
    {
        run_multiembed_lib <- function(block, embeddings_list, target_col, lib, pred, 
                                       num_to_merge, ranking = "rho")
        {
            pred_idx <- c()
            for(i in 1:NROW(pred))
            {
                pred_idx <- c(pred_idx, pred[i,1]:pred[i,2])
            }
            
            # run embeddings
            in_stats <- block_lnlp(block, lib = lib, pred = lib, num_neighbors = "e+1", 
                                   target_column = target_col, columns = embeddings_list, 
                                   first_column_time = TRUE, silent = TRUE)
            
            # ranking of embeddings
            if(ranking == "rho")
                in_sample_ranking <- order(in_stats$rho, decreasing = TRUE)
            else if(ranking == "mae")
                in_sample_ranking <- order(in_stats$mae, decreasing = FALSE)
            else
                in_sample_ranking <- order(in_stats$rmse, decreasing = FALSE)
            
            # get desired embeddings
            univar_embedding <- seq(from = target_col, to = target_col+E-1)
            multivar_embeddings <- embeddings_list[in_sample_ranking[1:num_to_merge],] # sorted
            
            # compute univariate and multivariate results
            temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = "e+1", 
                               target_column = target_col, first_column_time = TRUE, 
                               columns = rbind(univar_embedding, multivar_embeddings[1,]), 
                               stats_only = FALSE)
            
            results <- data.frame(time = temp[[1]]$model_output$time[pred_idx], 
                                  obs = temp[[1]]$model_output$obs[pred_idx], 
                                  univariate_pred = temp[[1]]$model_output$pred[pred_idx], 
                                  multivariate_pred = temp[[2]]$model_output$pred[pred_idx])
            
            # compute predictions for multiembed
            temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = 1, 
                               target_column = target_col, first_column_time = TRUE, 
                               columns =  multivar_embeddings[1:num_to_merge,], 
                               stats_only = FALSE)
            
            # calculate multiembed predictions
            for(j in 1:num_to_merge)
            {
                results[, paste("multiembed_", j, sep = "")] <- temp[[j]]$model_output$pred[pred_idx]
            }

            return(results)
        }
        
        n <- NROW(data)
        lib_splits <- list(1:3, 4:8, 9:12, 13:17)
        
        num_vars <- NCOL(data)
        block <- make_block(data, E, max_lag)
        embeddings_list <- make_embeddings(max_lag, E, num_vars)
        valid_embeddings <- apply(embeddings_list %% max_lag, 1, function (x) {1 %in% x})
        embeddings_list <- embeddings_list[valid_embeddings,]
        num_to_merge <- min(NROW(embeddings_list), num_to_merge)
        
        results <- do.call(rbind, lapply(1:length(lib_splits), function(k) {
            message(".", appendLF = FALSE)
            curr_pred <- lib[unlist(lib_splits[k]),]
            curr_lib <- lib[unlist(lib_splits[-k]),]
            return(run_multiembed_lib(block, embeddings_list, target, curr_lib, curr_pred, num_to_merge, ...))            
        }))

        return(results)
    }
    
    load("mesocosm_data.Rdata")
    
    start_time <- proc.time()
    message("Computing multiembed for mesocosm", appendLF = FALSE)
    data <- normed_data[, c("calanoid_copepods", "nanophytoplankton", "picophytoplankton")]
    calanoid_preds <- run_multiembed(data, lib, ...)
  
    data <- normed_data[, c("rotifers", "nanophytoplankton", "picophytoplankton")]
    rotifer_preds <- run_multiembed(data, lib, ...)
    message("done! (", round((proc.time()-start_time)[3], 3), " seconds elapsed.)")
        
    save(calanoid_preds, rotifer_preds, file = "mesocosm_multiembed_results.Rdata")
    return()
}

mesocosm_multiembed_stats <- function(denormalize = FALSE)
{
    process_mesocosm_stats <- function(preds, denormalize = FALSE)
    {
        if(denormalize) # denormalize time series
        {
            load("mesocosm_data.Rdata")
            mu <- mean(obs_data[,var_name], na.rm = TRUE)
            sigma <- sd(obs_data[,var_name], na.rm = TRUE)
            transform <- function(x) {(x * sigma + mu)^4}
            pred_cols <- names(preds)
            for(j in pred_cols[pred_cols != "time"])
            {
                preds[, j] <- transform(preds[, j])
            }
        }
        
        # average embeddings for multiembed
        cols_to_merge <- 5:NCOL(preds)
        preds$multiembed_pred <- if(length(cols_to_merge) > 1) rowMeans(preds[, cols_to_merge]) else preds[, cols_to_merge]
        
        return(list(univariate_stats = compute_stats(preds$obs, preds$univariate_pred), 
                    multivariate_stats = compute_stats(preds$obs, preds$multivariate_pred), 
                    multiembed_stats = compute_stats(preds$obs, preds$multiembed_pred)))
    }
    
    load("mesocosm_multiembed_results.Rdata")
    calanoid_stats <- process_mesocosm_stats(calanoid_preds, denormalize = denormalize)
    rotifer_stats <- process_mesocosm_stats(rotifer_preds, denormalize = denormalize)
        
    save(calanoid_stats, rotifer_stats, file = "mesocosm_multiembed_stats.Rdata")
    return()
}

plot_mesocosm_multiembed <- function(stats_file = "mesocosm_multiembed_stats.Rdata", 
                                     out_file = NULL, var = "rho", width = 2.5, height = 4)
{
    load(stats_file)
    
    if(!is.null(out_file))
        pdf(out_file, width = width, height = height)
    
    layout(matrix(c(1, 3, 2, 3), nrow = 2), heights = c(2.7, 1))
    par(mar = c(1, 2, 1, 1), oma = c(0.5, 2, 2, 0))
    
    for(stats in list(calanoid_stats, rotifer_stats))
    {
        to_plot <- c(stats$univariate_stats[, var], 
                     stats$multivariate_stats[, var], 
                     stats$multiembed_stats[, var])
        plot(to_plot[1], ylim = range(to_plot), xaxt = "n", col = "gray", pch = 8)
        points(to_plot[2], xaxt = "n", col = "red3", pch = 8)
        points(to_plot[3], xaxt = "n", col = "royalblue", pch = 8)
    }
    
    par(mar = c(0, 2, 0, 1))
    plot(c(0, 1), c(0, 1), type = "n", axes = FALSE)
    legend(0.5, 0.0, c("univariate", "best multivariate", "multiview embedding"), 
           col = c("gray", "red3", "royalblue"), pch = 8, xjust = 0.5, yjust = 0, horiz = FALSE)
    
    mtext(c("Calanoids", "Rotifers"), 3, line = 0.5, at = c(1/4, 3/4), outer = TRUE)
     mtext(var, 2, line = 0.5, at = 5/8, adj = 0.5, outer = TRUE)  
    
    if(!is.null(out_file))
        dev.off()
}

plot_mesocosm_multiembed_full <- function(stats_file = "mesocosm_multiembed_stats.Rdata", 
                                          out_file = NULL, width = 5, height = 5)
{
    load(stats_file)
    
    if(!is.null(out_file))
        pdf(out_file, width = width, height = height)
    
    layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7, 7), nrow = 3, byrow = TRUE), heights = c(3, 3, 1))
    par(mar = c(1, 2, 1, 1), oma = c(0.5, 2, 2, 0), las = 1)
    
    for(stats in list(calanoid_stats, rotifer_stats))
    {
        to_plot <- c(stats$univariate_stats$rho, 
                     stats$multivariate_stats$rho, 
                     stats$multiembed_stats$rho)
        plot(to_plot[1], ylim = range(to_plot), xaxt = "n", col = "gray", pch = 8)
        points(to_plot[2], xaxt = "n", col = "red3", pch = 8)
        points(to_plot[3], xaxt = "n", col = "royalblue", pch = 8)
        
        to_plot <- c(stats$univariate_stats$mae, 
                     stats$multivariate_stats$mae, 
                     stats$multiembed_stats$mae)
        plot(to_plot[1], ylim = range(to_plot), xaxt = "n", col = "gray", pch = 8)
        points(to_plot[2], xaxt = "n", col = "red3", pch = 8)
        points(to_plot[3], xaxt = "n", col = "royalblue", pch = 8)
        
        to_plot <- c(stats$univariate_stats$rmse, 
                     stats$multivariate_stats$rmse, 
                     stats$multiembed_stats$rmse)
        plot(to_plot[1], ylim = range(to_plot), xaxt = "n", col = "gray", pch = 8)
        points(to_plot[2], xaxt = "n", col = "red3", pch = 8)
        points(to_plot[3], xaxt = "n", col = "royalblue", pch = 8)
        
        # to_plot <- c(stats$univariate_stats$predictive_power, 
        #              stats$multivariate_stats$predictive_power, 
        #              stats$multiembed_stats$predictive_power)
        # plot(to_plot[1], ylim = range(to_plot), xaxt = "n", col = "gray", pch = 8)
        # points(to_plot[2], xaxt = "n", col = "red3", pch = 8)
        # points(to_plot[3], xaxt = "n", col = "royalblue", pch = 8)
    }
    
    par(mar = c(0, 2, 0, 1))
    plot(c(0, 1), c(0, 1), type = "n", axes = FALSE)
    legend(0.5, 0., c("univariate", "best multivariate", "multiview embedding"), 
           col = c("gray", "red3", "royalblue"), pch = 8, xjust = 0.5, yjust = 0, horiz = FALSE)
    
    mtext(c("Rotifers", "Calanoids"), 2, line = 0.5, at = c(5/14, 11/14), outer = TRUE)
    # mtext(c("accuracy (rho)", "error (MAE)", "error (RMSE)", "predictive power"), 
    #       3, line = 0.5, at = c(1/8, 3/8, 5/8, 7/8), adj = 0.5, outer = TRUE)  
    mtext(c("accuracy (rho)", "error (MAE)", "error (RMSE)"), 
          3, line = 0.5, at = c(1/6, 3/6, 5/6), adj = 0.5, outer = TRUE)  
    
    if(!is.null(out_file))
        dev.off()

    return()
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}