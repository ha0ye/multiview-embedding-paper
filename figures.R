plot_multiembed_perf_vs_nm <- function(stats_file = "merged_stats_1.Rdata", 
                                       plot_file = NULL, width = 8, height = 12)
{
    load(stats_file)
    
    stats <- output_stats
    stats$value = stats$rho
    to_plot <- summarySE(stats, measurevar = "value", groupvars = c("num_embeddings", "lib_size", "target"))
    my_plot <- ggplot(to_plot, aes(x = num_embeddings, y = value)) + 
#        geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
        geom_line(aes(y = value)) + 
        geom_line(aes(y = value - sd)) + 
        geom_line(aes(y = value + sd)) + 
        facet_wrap(~ lib_size + target, ncol = nlevels(as.factor(output_stats$target)), scales = "free") + 
#        scale_color_manual(values = c("gray40", "red3", "royalblue")) + 
        xlab("number of embeddings") + 
        ylab("forecast skill (rho)") + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              # legend.background = element_rect(color = "black"), 
              # legend.key = element_blank(), 
              # legend.title = element_blank(), 
              # legend.position = "bottom",  
              # legend.justification = c(0, 1), 
              panel.background = element_rect(color = "black", fill = NA))
    
    if(!is.null(plot_file))
    {
        pdf(plot_file, width = width, height = height)
    }
    print(my_plot)
    if(!is.null(plot_file))
    {
        dev.off()
    }
    return()
}

plot_combined_model_results <- function(var = "rho", out_file = NULL, width = 6.5, height = 6.5)
{
    # load results
    load("multiembed_all_results_1.Rdata")
    stats <- cbind(combined_stats, system = 1)
    load("multiembed_all_results_2.Rdata")
    stats <- rbind(stats, cbind(combined_stats, system = 2))
    load("multiembed_all_results_hw.Rdata")
    stats <- rbind(stats, cbind(combined_stats, system = 3))
    stats <- gather(stats, metric, value, univar_rho:multiembed_pp)
    metric_str <- strsplit(as.character(stats$metric), "_")
    stats$model <- sapply(metric_str, function(x) {x[1]})
    stats$metric <- sapply(metric_str, function(x) {x[2]})
    
    # fix levels
    stats$model <- factor(stats$model, levels = c("univar", "multivar", "multiembed"))
    levels(stats$model) <- c("univariate", "best multivariate", "multiview embedding")
    stats$metric <- factor(stats$metric, levels = c("rho", "pp", "mae", "rmse"))
    levels(stats$metric) <- c("rho", "PP", "MAE", "RMSE")
    
    # subset data and make plot
    stats <- stats[stats$metric == var,]
    to_plot <- summarySE(stats, measurevar = "value", groupvars = c("lib_size", "target", "model", "system"))
    my_plot <- ggplot(to_plot, aes(x = lib_size, y = value, color = model)) + 
        geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
        geom_line() + 
        facet_wrap(~ target + system, ncol = 3, scales = "free") + 
        scale_color_manual(values = c("gray40", "red3", "royalblue")) + 
        xlab("library size") + 
        ylab("forecast skill") + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              legend.background = element_rect(color = "black"), 
              legend.key = element_blank(), 
              legend.title = element_blank(), 
              legend.position = "bottom",  
              legend.justification = c(0, 1), 
              panel.background = element_rect(color = "black", fill = NA))
    
    if(!is.null(out_file))
        pdf(file = out_file, width = width, height = height)
    print(my_plot)
    if(!is.null(out_file))
        dev.off()
    return()   
}

plot_full_model_results <- function(data_file, out_file = NULL, width = 6.5, height = 9.5)
{
    load(data_file)
    
    stats <- gather(combined_stats, metric, value, univar_rho:multiembed_pp)
    metric_str <- strsplit(as.character(stats$metric), "_")
    stats$model <- sapply(metric_str, function(x) {x[1]})
    stats$metric <- sapply(metric_str, function(x) {x[2]})
    
    # fix levels
    stats$model <- factor(stats$model, levels = c("univar", "multivar", "multiembed"))
    levels(stats$model) <- c("univariate", "best multivariate", "multiview embedding")
    
    stats$metric <- factor(stats$metric, levels = c("rho", "mae", "rmse", "pp"))
    levels(stats$metric) <- c("rho", "MAE", "RMSE", "PP")
    
    to_plot <- summarySE(stats, measurevar = "value", groupvars = c("lib_size", "target", "metric", "model"))
    my_plot <- ggplot(to_plot, aes(x = lib_size, y = value, color = model)) + 
        geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
        geom_line() + 
        facet_wrap(metric ~ target, scales="free", nrow = 4) + 
        scale_color_manual(values = c("gray40", "red3", "royalblue")) + 
        xlab("library size") + 
        ylab("forecast skill") + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              legend.background = element_rect(color = "black"), 
              legend.key = element_blank(), 
              legend.title = element_blank(), 
              legend.position = "bottom",  
              legend.justification = c(0, 1), 
              panel.background = element_rect(color = "black", fill = NA))
    
    if(!is.null(out_file))
        pdf(file = out_file, width = width, height = height)
    print(my_plot)
    if(!is.null(out_file))
        dev.off()
    return()
}

plot_model_attractor <- function(data_file, t = 1000:3000, x = c(1,0), y = c(2,0), z = c(3,0), 
                                 x_label = "x", y_label = "y", z_label = "z", 
                                 plot_file = NULL)
{
    load(data_file)
    make_vec <- function(vec_definition)
    {
        return(data[t-vec_definition[2], vec_definition[1]])
    }
    t <- t + max(x[2], y[2], z[2])
    plot3d(make_vec(x), make_vec(y), make_vec(z), type = "l", 
           axes = FALSE, xlab = x_label, ylab = y_label, zlab = z_label)
    axes3d(edges = c("x--", "y+-", "z--"), tick = FALSE, labels = FALSE)
    m <- matrix(c(0.881, 0.473, 0.020, 0, 
                  -0.062, 0.071, 0.996, 0, 
                  0.470, -0.878, 0.092, 0, 
                  0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(userMatrix = m, windowRect = c(100, 100, 900, 800))
    if(!is.null(plot_file))
    {
        rgl.postscript(plot_file, fmt = "pdf")
    }
    return()
}

make_HW_plot_panel <- function(data_file, t = 500:1000, x = c(1,0), y = c(2,0), z = c(4,0), 
                               col = rgb(1.0, 0.0, 0.0), 
                               x_label = "x", y_label = "y", z_label = "z", 
                               plot_file = NULL, m = NULL, edges = c("x-+", "y--", "z--"))
{
    make_vec <- function(vec_definition)
    {
        return(data[t-vec_definition[2], vec_definition[1]])
    }
    t <- t + max(x[2], y[2], z[2])
    
    if(is.null(m))
    {
        m <- matrix(c(0.015, 0.598, -0.801, 0, 
                      0.930, 0.287, 0.231, 0, 
                      0.368, -0.748, -0.552, 0, 
                      0, 0, 0, 1), byrow = TRUE, nrow = 4)
    }
    
    load(data_file)
    xx <- make_vec(x)
    yy <- make_vec(y)
    zz <- make_vec(z)
    
    plot3d(NA, NA, NA, 
           xlim = range(xx), ylim = range(yy), zlim = range(zz), 
           xlab = x_label, ylab = y_label, zlab = z_label, axes = FALSE)
    axes3d(edges = edges)
    box3d()
    par3d(userMatrix = m, zoom = 1.25, windowRect = c(50, 50, 850, 850))
    lines3d(xx, yy, zz, col = col, lwd = 2)
    if(!is.null(plot_file))
    {
        rgl.postscript(plot_file, fmt = "pdf")
    }
    return()
}

# attractor plots of food chain model
if(FALSE)
{
    plot_model_attractor("model_2_alt.Rdata")
    plot_model_attractor("model_2_alt.Rdata", y = c(1,6), z = c(1,12))
    plot_model_attractor("model_2_alt.Rdata", x = c(2,0), y = c(2,6), z = c(2,12))
    plot_model_attractor("model_2_alt.Rdata", x = c(3,0), y = c(3,6), z = c(3,12))
    plot_model_attractor("model_2_alt.Rdata", x = c(1,0), y = c(2,0), z = c(1,6))
}

# time series plot for Huisman-Weissing model
if(FALSE)
{
    num_species <- length(grep("N[0-9]", names(data)))
    N_range <- range(data[,1:num_species])
    
    palette <- rainbow(num_species)
    t <- 1:NROW(data)
    plot(x = NA, type = "n", xlim = range(t), ylim = N_range, 
         xlab = "time", ylab = "abundance")
    for(i in 1:num_species)
    {
        lines(t, data[t,i], col = palette[i])
    }
}

# attractor plots for Huisman-Weissing model
if(FALSE)
{
    load("model_hw.Rdata")
    t <- 500:1000
    palette <- c(rgb(1.0, 0.0, 0.0), rgb(0.0, 0.0, 1.0))
    point_color <- palette[1+(data[t,4] > data[t,2])]
    
    make_HW_plot_panel("model_hw.Rdata", t, col = point_color, 
                    plot_file = "hw_native.pdf")
    
    m <- matrix(c(0.707, 0.707, 0, 0, 
                  -0.171, 0.181, 0.968, 0, 
                  0.684, -0.685, 0.249, 0, 
                  0, 0, 0, 1), byrow = TRUE, nrow = 4)
    
    make_HW_plot_panel("model_hw.Rdata", t, col = point_color, 
                    x = c(2, 0), y = c(2, 5), z = c(2, 10), 
                    m = m, edges = c("x--", "y+-", "z--"), 
                    x_label = "", y_label = "", z_label = "", 
                    plot_file = "hw_univar_2.pdf")
    
    make_HW_plot_panel("model_hw.Rdata", t, col = point_color, 
                    x = c(4, 0), y = c(4, 5), z = c(4, 10), 
                    m = m, edges = c("x--", "y+-", "z--"), 
                    x_label = "", y_label = "", z_label = "", 
                    plot_file = "hw_univar_4.pdf")
}


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
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

compute_forecast_errors <- function(lib = c(551, 600), pred = c(3001, 4000), 
                                    embedding_1 = 1:3, embedding_2 = 4:6, embedding_3 = 7:9)
{
    get_forecasts <- function(embedding, nn = (length(embedding)+1))
    {
        temp_x <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = nn, 
                             columns = embedding, target_column = 1, stats_only = FALSE)
        temp_y <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = nn, 
                             columns = embedding, target_column = 4, stats_only = FALSE)
        temp_z <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = nn, 
                             columns = embedding, target_column = 7, stats_only = FALSE)       
        return(cbind(temp_x[[1]]$model_output$pred, 
                     temp_y[[1]]$model_output$pred, 
                     temp_z[[1]]$model_output$pred))
    }
    
    get_multiembed_forecasts <- function()
    {
        embeddings_list <- t(combn(9, 3, simplify = TRUE))
        valid_embeddings <- apply(embeddings_list %% 3, 1, function (x) {1 %in% x})
        embeddings_list <- embeddings_list[valid_embeddings,]
        
        output_pred <- matrix(NA, nrow = NROW(block), ncol = 3)
        for(j in 1:3)
        {
            target_column <- c(1, 4, 7)[j]
            in_stats <- block_lnlp(block, lib = lib, pred = lib, num_neighbors = "e+1", 
                                   target_column = target_column, columns = embeddings_list, 
                                   silent = TRUE)
            
            in_sample_ranking <- order(in_stats$rho, decreasing = TRUE)
            multivar_embeddings <- embeddings_list[in_sample_ranking[1:8],]
            
            temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = 1, 
                               target_column = target_column, columns =  multivar_embeddings, 
                               stats_only = FALSE)
            
            # calculate multiembed predictions
            top_embeddings_preds <- do.call(cbind, lapply(temp, function(x) x$model_output$pred))
            output_pred[,j] <- rowMeans(top_embeddings_preds)
        }
        return(output_pred)
    }
    
    load("model_2.Rdata")
    orig_data <- normalize(data)
    data <- normed_data
    n <- NROW(data)
    block <- cbind(data[,1], c(NA, data[1:(n-1),1]), c(NA, NA, data[1:(n-2),1]), 
                   data[,2], c(NA, data[1:(n-1),2]), c(NA, NA, data[1:(n-2),2]), 
                   data[,3], c(NA, data[1:(n-1),3]), c(NA, NA, data[1:(n-2),3]))
    
    # compute errors for x reconstruction
    mx_forecasts <- get_forecasts(embedding_1)
    mx_errors <- mx_forecasts[pred[1]:(pred[2]-1),] - data[(pred[1]+1):pred[2],]
    mx_err_dist <- sqrt(rowSums(mx_errors * mx_errors))
    
    # compute errors for y reconstruction
    my_forecasts <- get_forecasts(embedding_2)
    my_errors <- my_forecasts[pred[1]:(pred[2]-1),] - data[(pred[1]+1):pred[2],]
    my_err_dist <- sqrt(rowSums(my_errors * my_errors))
    
    # compute errors for z reconstruction
    mz_forecasts <- get_forecasts(embedding_3)
    mz_errors <- mz_forecasts[pred[1]:(pred[2]-1),] - data[(pred[1]+1):pred[2],]
    mz_err_dist <- sqrt(rowSums(mz_errors * mz_errors))
    
    # compute errors for multiembed
    multi_forecasts <- get_multiembed_forecasts()
    multi_errors <- multi_forecasts[pred[1]:(pred[2]-1),] - data[(pred[1]+1):pred[2],]
    multi_err_dist <- sqrt(rowSums(multi_errors * multi_errors))
    
    #     mx_forecasts_1 <- get_forecasts(embedding_1, nn = 1)
    #     my_forecasts_1 <- get_forecasts(embedding_2, nn = 1)
    #     mz_forecasts_1 <- get_forecasts(embedding_3, nn = 1)
    #     mxy_forecasts <- (mx_forecasts_1 + my_forecasts_1 + mz_forecasts_1) / 3
    #     mxy_errors <- mxy_forecasts[pred[1]:(pred[2]-1),] - data[(pred[1]+1):pred[2],]
    #     mxy_err_dist <- sqrt(rowSums(mxy_errors * mxy_errors))
    
    save(mx_forecasts, mx_err_dist, my_forecasts, my_err_dist, 
         mz_forecasts, mz_err_dist, multi_forecasts, multi_err_dist, 
         pred, data, orig_data, file = "pred_errors.Rdata")
}

make_error_plot <- function(to_plot = NULL, pred_file = NULL)
{
    load("pred_errors.Rdata")
    if(is.null(to_plot))
    {
        to_plot <- pred
        layout3d(matrix(1:4, nrow = 2, byrow = TRUE))
    }
    t <- (to_plot[1]+1):to_plot[2]
    pred_t <- to_plot[1]:(to_plot[2]-1)
    
    user_matrix <- matrix(c(0.6, -0.8, 0, 0, 
                            0.25, 0.2, 1, 0, 
                            -0.75, -0.55, 0.3, 0, 
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(windowRect = c(0, 70, 600, 670))
    pt_size <- 4
    
    num_colors <- 30
    brks <- seq(length.out = num_colors, from = 0, 
                to = 3)#max(c(mx_err_dist, my_err_dist, mz_err_dist, multi_err_dist)))
    
    # plot M with Mx errors
    mx_pal <- colorRampPalette(c("red", "red"))
    mx_col <- mx_pal(num_colors)[as.numeric(cut(mx_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(mx_forecasts[pred_t,1], mx_forecasts[pred_t,2], mx_forecasts[pred_t,3], size = pt_size, col = mx_col)
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2a.pdf", "pdf")
    
    # plot M with My errors
    my_pal <- colorRampPalette(c("green", "white"))
    my_col <- my_pal(num_colors)[as.numeric(cut(my_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(my_forecasts[pred_t,1], my_forecasts[pred_t,2], my_forecasts[pred_t,3], size = pt_size, col = my_col)
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2d.pdf", "pdf")
    
    # plot M with Mz errors
    mz_pal <- colorRampPalette(c("blue", "blue"))
    mz_col <- mz_pal(num_colors)[as.numeric(cut(mz_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(mz_forecasts[pred_t,1], mz_forecasts[pred_t,2], mz_forecasts[pred_t,3], size = pt_size, col = mz_col)
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2b.pdf", "pdf")
    
    # plot M with multiembed errors
    multi_pal <- colorRampPalette(c("magenta", "magenta"))
    multi_col <- multi_pal(num_colors)[as.numeric(cut(multi_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(multi_forecasts[pred_t,1], multi_forecasts[pred_t,2], multi_forecasts[pred_t,3], size = pt_size, col = multi_col)
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2c.pdf", "pdf")
    
    if(!is.null(pred_file))
        rgl.postscript(pred_file, "pdf")
    return()
}