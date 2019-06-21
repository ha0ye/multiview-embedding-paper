plot_metric_vs_lib_size <- function(obs_err_var = 0.1, metric = "rho", stats_file = "combined_stats.Rdata")
{
    load(stats_file)
    stats_df <- stats %>% 
        filter(err == obs_err_var)
    stats_df$metric <- stats_df[, metric]
    
    to_plot <- stats_df %>% group_by(target, model, lib_size, method) %>% 
        summarize(avg = mean(metric), uq = quantile(metric, 0.75), lq = quantile(metric, 0.25))
    
    return(ggplot(to_plot, aes(x = lib_size, y = avg, color = method)) + 
               geom_line(size = 1) + 
               scale_x_continuous(limits = c(25, 100)) + 
               geom_line(aes(y = lq), linetype = 2, size = 0.25) + 
               geom_line(aes(y = uq), linetype = 2, size = 0.25) + 
#               geom_errorbar(aes(ymin = lq, ymax = uq)) + 
               facet_wrap(~ model + target, ncol = 5, scales = "free", dir = "v") + 
               scale_color_manual(values = c("gray40", "red3", "royalblue")) + 
               xlab("Library Size") + 
               ylab("Forecast Skill") + 
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(color = "black"), 
                     legend.background = element_rect(color = "black"), 
                     legend.key = element_blank(), 
                     legend.title = element_blank(), 
                     legend.position = "bottom",  
                     legend.justification = c(0, 1), 
                     panel.background = element_rect(color = "black", fill = NA)))
}

plot_perf_vs_lib_size <- function(system = 1, obs_err_var = 0.1, stats_file = "combined_stats.Rdata")
{
    load(stats_file)
    stats_df <- stats %>% 
        filter(model == system, err == obs_err_var) %>%
        tidyr::gather(metric, value, c(mae, rho, rmse))
    stats_df$metric <- as.factor(stats_df$metric)
    
    to_plot <- stats_df %>% group_by(target, lib_size, method, metric) %>% 
        summarize(avg = mean(value), uq = quantile(value, 0.75), lq = quantile(value, 0.25))
    to_plot$metric <- factor(to_plot$metric, levels(to_plot$metric)[c(2,1,3)])
    
    return(ggplot(to_plot, aes(x = lib_size, y = avg, color = method)) + 
               geom_line(size = 1) + 
               scale_x_continuous(limits = c(25, 100)) + 
               geom_line(aes(y = lq), linetype = 2, size = 0.25) + 
               geom_line(aes(y = uq), linetype = 2, size = 0.25) + 
#               geom_errorbar(aes(ymin = lq, ymax = uq)) + 
               facet_wrap(~ metric + target, ncol = 3, scales = "free", dir = "v") + 
               scale_color_manual(values = c("gray40", "red3", "royalblue")) + 
               xlab("Library Size") + 
               ylab("Forecast Skill") + 
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(color = "black"), 
                     legend.background = element_rect(color = "black"), 
                     legend.key = element_blank(), 
                     legend.title = element_blank(), 
                     legend.position = "bottom",  
                     legend.justification = c(0, 1), 
                     panel.background = element_rect(color = "black", fill = NA)))
}

plot_metric_vs_noise <- function(system = 1, metric = "mae", lib_sizes = c(25, 50, 100, 300), 
                                 stats_file = "combined_stats.Rdata")
{
    load(stats_file)
    stats_df <- stats %>%
        filter(model == system, lib_size %in% lib_sizes)
    stats_df$metric <- stats_df[, metric]
    
    univar_avg_errs <- stats_df %>% 
        filter(method == "univariate") %>%
        group_by(target, lib_size, err) %>%
        summarize(univar_scale = mean(metric))
    avg_errs <- stats_df %>% 
        group_by(target, lib_size, method, err) %>%
        summarize(avg = mean(metric), uq = quantile(metric, 0.75), lq = quantile(metric, 0.25))
    to_plot <- left_join(avg_errs, univar_avg_errs, 
                         by = c("target", "lib_size", "err"))
    to_plot$avg <- to_plot$avg / to_plot$univar_scale
    to_plot$lq <- to_plot$lq / to_plot$univar_scale
    to_plot$uq <- to_plot$uq / to_plot$univar_scale
    
    return(ggplot(to_plot, aes(x = err, y = avg, color = method)) + 
               geom_line(size = 1) + 
               geom_line(aes(y = lq, color = method), linetype = 2, size = 0.25) + 
               geom_line(aes(y = uq, color = method), linetype = 2, size = 0.25) + 
               #geom_errorbar(aes(ymin = lq, ymax = uq)) + 
               facet_wrap(~ lib_size + target, ncol = length(lib_sizes), scales = "free", dir = "v") + 
               scale_color_manual(values = c("gray40", "red3", "royalblue")) + 
               xlab("Noise (fraction added variance)") + 
               ylab("Forecast Error (normalized MAE)") + 
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(color = "black"), 
                     legend.background = element_rect(color = "black"), 
                     legend.key = element_blank(), 
                     legend.title = element_blank(), 
                     legend.position = "bottom",  
                     legend.justification = c(0, 1), 
                     panel.background = element_rect(color = "black", fill = NA)))
}

plot_perf_vs_k <- function(system = 1, lib_sizes = c(25, 50, 100, 300), 
                           stats_file = "combined_stats_k.Rdata", my_metric = "rho")
{
    load(stats_file)
    stats_df <- stats %>%
        filter(model == system, lib_size %in% lib_sizes, metric == my_metric)

    to_plot <- stats_df %>% 
        group_by(target, k, metric, lib_size) %>%
        summarize(avg = mean(value), uq = quantile(value, 0.75), lq = quantile(value, 0.25))

    return(ggplot(to_plot, aes(x = k, y = avg)) + 
               geom_line(size = 1, color = "royalblue") + 
               geom_line(aes(y = lq), linetype = 2, size = 0.25, color = "royalblue") + 
               geom_line(aes(y = uq), linetype = 2, size = 0.25, color = "royalblue") + 
               facet_wrap(~ lib_size + target, ncol = length(lib_sizes), scales = "free", dir = "v") + 
               xlab("Number of MVE Embeddings Used (k)") + 
               ylab("Forecast Skill") + 
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(color = "black"), 
                     legend.background = element_rect(color = "black"), 
                     legend.key = element_blank(), 
                     legend.title = element_blank(), 
                     legend.position = "bottom",  
                     legend.justification = c(0, 1), 
                     panel.background = element_rect(color = "black", fill = NA)))
}

plot_sparseness_model_2 <- function(pred_target = 438, pred = c(401, 525), 
                                    lib = c(501, 525), 
                                    user_matrix = diag(nrow = 4), 
                                    make_prediction = FALSE, 
                                    plot_file = NULL)
{
    distance <- function(x, y)
    {
        return(sqrt(sum((x-y)^2)))
    }
    
    add_vector <- function(t, vec_length = 0.5, line_width = 2, col = "black")
    {
        start_pos <- fine_data[t,]
        vec <- fine_data[t+1,] - start_pos
        vec <- vec * vec_length / sqrt(sum(vec^2))
        lines3d(rbind(start_pos, start_pos + vec), 
                lwd = line_width, col = col)
        return()
    }
    
    load("model_2.Rdata")
    mu <- colMeans(data)
    sigma <- apply(data, 2, sd)
    orig_data <- normalize(data)
    
    load("model_2_fine.Rdata")
    fine_data <- data - matrix(rep(mu, each = NROW(data)), nrow = NROW(data))
    fine_data <- fine_data / matrix(rep(sigma, each = NROW(data)), nrow = NROW(data))
    
    pt_size <- 8
    
    # plot attractor
    line_color <- rgb(204, 204, 204, alpha = 128, maxColorValue = 256)
    plot3d(fine_data[(pred[1]*16-15):(pred[2]*16-15),], type = "l", col = line_color,
           xlab = "", ylab = "", zlab = "", axes = FALSE)

    # plot predictee and predictee motion
    points3d(orig_data[pred_target,1], orig_data[pred_target,2], orig_data[pred_target,3], size = pt_size)
    add_vector(pred_target*16-15)
    
    # find neighbors
    lib_idx <- lib[1]:lib[2]
    distances_to_lib <- sapply(lib_idx, function(t) {distance(orig_data[pred_target,], fine_data[16*t-12,])})
    neighbors <- lib_idx[rank(distances_to_lib) <= 4]
    neighbor_dist <- distances_to_lib[rank(distances_to_lib) <= 4]
    
    # plot neighbors and neighbors motion
    points3d(fine_data[neighbors*16-12,], size = pt_size, col = "red")
#    points3d(orig_data[setdiff(lib_idx, neighbors),], size = pt_size, col = "blue")
    points3d(fine_data[setdiff(lib_idx, neighbors)*16-12,], size = pt_size, col = "blue")
    sapply(neighbors*16-12, add_vector, col = "red")
    
    # plot prediction vector
    if(make_prediction)
    {
        weights <- exp(-neighbor_dist / neighbor_dist[1])
        prediction <- weights %*% fine_data[neighbors*16-11,] / sum(weights)
        points3d(prediction, size = pt_size, col = "green")
    }
    
    
    box3d()
    par3d(windowRect = c(0, 70, 600, 670), userMatrix = user_matrix)
    
    if(!is.null(plot_file))
        rgl.postscript(plot_file, "pdf")
    return()
}

plot_raw_time_series <- function(lib = c(501, 525), 
                                 plot_file = NULL, 
                                 width = 6, height = 6)
{
    t <- (lib[1]*16-15):(lib[2]*16-15)
    load("model_2_fine.Rdata")
    
    if(!is.null(plot_file))
        pdf(plot_file, width = width, height = height)
    
    par(mfrow = c(3,1), mar = c(0.5,4,0.5,0.5), mgp = c(2.5,1,0), las = 1)
    plot(data[t,1], type = "l", col = "red", xaxt = "n", xlab = "", ylab = "x")
    plot(data[t,2], type = "l", col = "green", xaxt = "n", xlab = "", ylab = "y")
    plot(data[t,3], type = "l", col = "blue", xaxt = "s", xlab = "", ylab = "z")
    
    if(!is.null(plot_file))
        dev.off()
    return()
}

plot_attractor_model_2 <- function(spline_t = c(401, 525), spline_col = "black", 
                                   point_t = spline_t, point_col = NULL, point_size = 8, 
                                   user_matrix = diag(nrow = 4), 
                                   embedding = c(1, 0, 2, 0, 3, 0), 
                                   plot_file = NULL)
{
    # setup columns
    load("model_2_fine.Rdata")
    t <- (spline_t[1]*16-15):(spline_t[2]*16-15)
    spline_block <- cbind(data[t - embedding[2]*16, embedding[1]], 
                          data[t - embedding[4]*16, embedding[3]], 
                          data[t - embedding[6]*16, embedding[5]])
    load("model_2.Rdata")
    t <- (point_t[1]):(point_t[2])
    point_block <- cbind(data[t - embedding[2], embedding[1]], 
                         data[t - embedding[4], embedding[3]], 
                         data[t - embedding[6], embedding[5]])
    
    # plot attractor
    plot3d(spline_block, type = "l", col = spline_col,
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    if(!is.null(point_col))
        points3d(point_block, col = point_col, size = point_size)
    
    box3d()
    par3d(windowRect = c(0, 70, 600, 670), userMatrix = user_matrix)
    
    if(!is.null(plot_file))
        rgl.postscript(plot_file, "pdf")
    return()
}

plot_attractor_and_forecasts_model_2 <- function(points, point_col = "black", 
                                                 spline_t = c(401, 525), spline_col = "black", 
                                                 point_size = 8, 
                                                 user_matrix = diag(nrow = 4), 
                                                 plot_file = NULL)
{
    load("model_2.Rdata")
    mu <- colMeans(data)
    sigma <- apply(data, 2, sd)
    orig_data <- normalize(data)
    
    load("model_2_fine.Rdata")
    fine_data <- data - matrix(rep(mu, each = NROW(data)), nrow = NROW(data))
    fine_data <- fine_data / matrix(rep(sigma, each = NROW(data)), nrow = NROW(data))
    
    t <- (spline_t[1]*16-15):(spline_t[2]*16-15)
    spline_block <- fine_data[t,]
    
    # plot attractor
    plot3d(spline_block, type = "l", col = spline_col,
           xlab = "", ylab = "", zlab = "", axes = FALSE)    
    
    # plot points
    points3d(points, size = point_size, col = point_col)
    
    box3d()
    par3d(windowRect = c(0, 70, 600, 670), userMatrix = user_matrix)
    
    if(!is.null(plot_file))
        rgl.postscript(plot_file, "pdf")
    return()
}

make_error_attractor_ts_plot <- function(lib = c(551, 600), pred = c(2001, 3000), plot_file = NULL)
{
    load("pred_errors.Rdata")
    
    t <- lib[1]:lib[2]
    t1 <- (pred[1]+1):pred[2]
    user_matrix <- matrix(c(0.6, -0.8, 0, 0,
                            0.25, 0.2, 1, 0,
                            -0.75, -0.55, 0.3, 0,
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(windowRect = c(0, 70, 600, 670), userMatrix = user_matrix)
    pt_size <- 4
    plot3d(orig_data[t1,1], orig_data[t1,2], orig_data[t1, 3], type = "l", col = "#CCCCCC",
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(orig_data[t,1], orig_data[t,2], orig_data[t,3], col = "#000000")
    box3d()
    rgl.postscript("figure_2z.pdf", "pdf")
    
    pdf("figure_2zz.pdf", width = 6, height = 6)
    par(mfrow = c(3,1), mar = c(4,4,1,1), mgp = c(2.5,1,0), las = 1)
    plot(orig_data[t,1], type = "l", col = "red", xaxt = "n", xlab = "", ylab = "x")
    plot(orig_data[t,2], type = "l", col = "green", xaxt = "n", xlab = "", ylab = "y")
    plot(orig_data[t,3], type = "l", col = "blue", xaxt = "s", xlab = "", ylab = "z")
    dev.off()
    
    return()   
}

make_error_plot <- function(to_plot = NULL, plot_file = NULL)
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
    
    # plot M with Mx errors
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC",
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(mx_forecasts[pred_t,1], mx_forecasts[pred_t,2], mx_forecasts[pred_t,3], 
             size = pt_size, col = "red")
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2a.pdf", "pdf")
    
    # plot M with My errors
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC",
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(my_forecasts[pred_t,1], my_forecasts[pred_t,2], my_forecasts[pred_t,3], 
             size = pt_size, col = "green")
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2d.pdf", "pdf")
    
    # plot M with Mz errors
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC",
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(mz_forecasts[pred_t,1], mz_forecasts[pred_t,2], mz_forecasts[pred_t,3], 
             size = pt_size, col = "blue")
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2b.pdf", "pdf")
    
    # plot M with multiembed errors
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC",
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(multi_forecasts[pred_t,1], multi_forecasts[pred_t,2], multi_forecasts[pred_t,3], 
             size = pt_size, col = "magenta")
    box3d()
    par3d(userMatrix = user_matrix)
    rgl.postscript("figure_2c.pdf", "pdf")
    
    if(!is.null(plot_file))
        rgl.postscript(plot_file, "pdf")
    return()
}

plot_noise_model_2 <- function(pred_target = 438, pred = c(401, 525), 
                               lib = c(238, 638), 
                               plot_file = NULL)
{
    distance <- function(x, y)
    {
        return(sqrt(sum((x-y)^2)))
    }
    
    add_vector <- function(t, vec_length = 0.5, line_width = 2, col = "black")
    {
        start_pos <- fine_data[t*16-15,]
        vec <- fine_data[t*16-14,] - start_pos
        vec <- vec * vec_length / sqrt(sum(vec^2))
        lines3d(rbind(start_pos, start_pos + vec), 
                lwd = line_width, col = col)
        return()
    }
    
    add_vector_noisy <- function(t, vec_length = 0.5, line_width = 2, col = "black")
    {
        start_pos <- noisy_data[t,]
        vec <- fine_data[t*16-14,] - start_pos + rnorm(3, sd = sqrt(0.1))
        vec <- vec * vec_length / sqrt(sum(vec^2))
        lines3d(rbind(start_pos, start_pos + vec), 
                lwd = line_width, col = col)
        return()
    }
    
    load("model_2.Rdata")
    mu <- colMeans(data)
    sigma <- apply(data, 2, sd)
    orig_data <- normalize(data)
    noisy_data <- orig_data + rnorm(length(orig_data), sd = sqrt(0.1))
    
    load("model_2_fine.Rdata")
    fine_data <- data - matrix(rep(mu, each = NROW(data)), nrow = NROW(data))
    fine_data <- fine_data / matrix(rep(sigma, each = NROW(data)), nrow = NROW(data))
    
    user_matrix <- matrix(c(0.406, -0.914, -0.008, 0,
                            0.485, 0.208, 0.849, 0,
                            -0.775, -0.349, 0.528, 0,
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    pt_size <- 4
    
    # plot attractor
    line_color <- rgb(204, 204, 204, alpha = 128, maxColorValue = 256)
    plot3d(fine_data[(pred[1]*16-15):(pred[2]*16-15),], type = "l", col = line_color,
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    
    # plot predictee and predictee motion
    points3d(orig_data[pred_target,1], orig_data[pred_target,2], orig_data[pred_target,3], 
             size = pt_size)
    add_vector(pred_target)
    
    # find neighbors
    lib_idx <- lib[1]:lib[2]
    distances_to_lib <- sapply(lib_idx, function(t) {distance(orig_data[pred_target,], noisy_data[t,])})
    neighbors <- lib_idx[rank(distances_to_lib) <= 20]
    
    # plot neighbors and neighbors motion
    points3d(noisy_data[neighbors,], size = pt_size, col = "red")
    points3d(noisy_data[setdiff(lib_idx, neighbors),], size = pt_size, col = "blue")
    sapply(neighbors, add_vector_noisy, col = "red")
    
    box3d()
    par3d(windowRect = c(0, 70, 600, 670), userMatrix = user_matrix)
    
    if(!is.null(plot_file))
        rgl.postscript(plot_file, "pdf")
    return()
}

plot_multiview_model_2 <- function(pred_target = 438, pred = c(401, 525), 
                                   lib = c(501, 600), pt_size = 4)
{
    make_panel <- function(columns = 1:3, 
                           user_matrix = matrix(c(0.406, -0.914, -0.008, 0,
                                                  0.485, 0.208, 0.849, 0,
                                                  -0.775, -0.349, 0.528, 0,
                                                  0, 0, 0, 1), nrow = 4, byrow = TRUE), 
                           neighbors = NULL, neighbor_colors = NULL)
    {
        add_vector <- function(t, vec_length = 0.3, line_width = 2, col = "black")
        {
            start_pos <- reconstruction[t,]
            vec <- reconstruction[t+1,] - start_pos
            vec <- vec * vec_length / sqrt(sum(vec^2))
            lines3d(rbind(start_pos, start_pos + vec), 
                    lwd = line_width, col = col)
            return()
        }
        
        pred_t <- (pred[1]*16-15):(pred[2]*16-15)
        pred_target_t <- pred_target*16-15
        line_color <- rgb(204, 204, 204, alpha = 128, maxColorValue = 256)
        reconstruction <- block[, columns]
        plot3d(reconstruction[pred_t,], type = "l", col = line_color,
               xlab = "", ylab = "", zlab = "", axes = FALSE)
        box3d()
        par3d(windowRect = c(0, 70, 600, 670), userMatrix = user_matrix)
        
        points3d(reconstruction[pred_target_t,1], 
                 reconstruction[pred_target_t,2], 
                 reconstruction[pred_target_t,3], 
                 size = pt_size)
        add_vector(pred_target_t)
        
        # find neighbors
        lib_idx <- lib[1]:lib[2]
        distances_to_lib <- sapply(lib_idx, function(t) {
            distance(reconstruction[pred_target_t,], reconstruction[t*16-12,])})
        if(is.null(neighbors))
            neighbors <- lib_idx[rank(distances_to_lib) == 1]
        
        # plot neighbors and neighbors motion
        if(is.null(neighbor_colors))
        {
            points3d(reconstruction[neighbors*16-12, 1], 
                     reconstruction[neighbors*16-12, 2], 
                     reconstruction[neighbors*16-12, 3], 
                     size = pt_size, col = "red")
            sapply(neighbors*16-12, add_vector, col = "red")
        } else {
            for(i in seq_along(neighbors))
            {
                points3d(reconstruction[neighbors[i]*16-12, 1], 
                         reconstruction[neighbors[i]*16-12, 2], 
                         reconstruction[neighbors[i]*16-12, 3], 
                         size = pt_size, col = neighbor_colors[i])
                add_vector(neighbors[i]*16-12, col = neighbor_colors[i])
            }
        }
        points3d(reconstruction[setdiff(lib_idx, neighbors)*16-12,], size = pt_size, col = "blue")
        
        return(neighbors)
    }
    
    distance <- function(x, y)
    {
        return(sqrt(sum((x-y)^2)))
    }
    
    load("model_2.Rdata")
    mu <- colMeans(data)
    sigma <- apply(data, 2, sd)
    orig_data <- normalize(data)
    noisy_data <- orig_data + rnorm(length(orig_data), sd = sqrt(0.1))
    
    load("model_2_fine.Rdata")
    fine_data <- data - matrix(rep(mu, each = NROW(data)), nrow = NROW(data))
    fine_data <- fine_data / matrix(rep(sigma, each = NROW(data)), nrow = NROW(data))
    n <- NROW(fine_data)
    block <- cbind(fine_data[,1], 
                   c(rep.int(NA, 16), fine_data[1:(n-16), 1]), 
                   c(rep.int(NA, 32), fine_data[1:(n-32), 1]), 
                   fine_data[,2], 
                   c(rep.int(NA, 16), fine_data[1:(n-16), 2]), 
                   c(rep.int(NA, 32), fine_data[1:(n-32), 2]), 
                   fine_data[,3], 
                   c(rep.int(NA, 16), fine_data[1:(n-16), 3]), 
                   c(rep.int(NA, 32), fine_data[1:(n-32), 3]))
    
    # pdf("model_2_ts.pdf", width = 6, height = 10)
    # par(mfrow = c(9,1), mar = c(1,1,0,0), lwd = 2)
    # my_cols <- rainbow(9)
    # idx <- sample(1:9,9)
    # for(i in 1:9)
    # {
    #     plot(block[401:1200,idx[i]], col = my_cols[i], type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    # }
    # dev.off()
    
    
    # nn_1 <- make_panel(columns = c(2, 5, 6), user_matrix = matrix(c(0.992, -0.025, -0.127, 0, 
    #                                                           -0.019, 0.942, -0.335, 0, 
    #                                                           0.128, 0.225, 0.933, 0, 
    #                                                           0, 0, 0, 1), nrow = 4, byrow = TRUE))
    # rgl.postscript("figure_cc.pdf", "pdf")
    # 
    # nn_1 <- make_panel(columns = c(3, 7, 9))
    # rgl.postscript("figure_dd.pdf", "pdf")
    
    
    nn_1 <- make_panel(columns = 1:3, user_matrix = matrix(c(0.992, -0.025, -0.127, 0, 
                                                     -0.019, 0.942, -0.335, 0, 
                                                     0.128, 0.225, 0.933, 0, 
                                                     0, 0, 0, 1), nrow = 4, byrow = TRUE))
    rgl.postscript("figure_1ci.pdf", "pdf")
    nn_2 <- make_panel(columns = 4:6, user_matrix = matrix(c(0.491, 0.064, -0.869, 0, 
                                                     0.417, -0.893, 0.170, 0, 
                                                     -0.766, -0.446, -0.464, 0, 
                                                     0, 0, 0, 1), nrow = 4, byrow = TRUE))
    rgl.postscript("figure_1cii.pdf", "pdf")
    nn_3 <- make_panel(columns = 7:9, user_matrix = matrix(c(0.904, 0.428, -0.031, 0, 
                                                     -0.094, 0.267, 0.959, 0, 
                                                     0.418, -0.864, 0.282, 0, 
                                                     0, 0, 0, 1), nrow = 4, byrow = TRUE))
    rgl.postscript("figure_1ciii.pdf", "pdf")
    
    make_panel(columns = c(1,4,7), 
               neighbors = c(nn_1, nn_2, nn_3), 
               neighbor_colors = c("red", "green", "blue"))
    rgl.postscript("figure_1civ.pdf", "pdf")
    
    return()
}

plot_multiview_model_2_new_fig_1 <- function(pred_target = 438, pred = c(401, 525), 
                                   lib = c(501, 600), pt_size = 4)
{
    make_panel <- function(columns = 1:3, 
                           user_matrix = matrix(c(0.406, -0.914, -0.008, 0,
                                                  0.485, 0.208, 0.849, 0,
                                                  -0.775, -0.349, 0.528, 0,
                                                  0, 0, 0, 1), nrow = 4, byrow = TRUE))
    {
        pred_t <- (pred[1]*16-15):(pred[2]*16-15)
        pred_target_t <- pred_target*16-15
        line_color <- rgb(204, 204, 204, alpha = 128, maxColorValue = 256)
        reconstruction <- block[, columns]
        plot3d(reconstruction[pred_t,], type = "l", col = line_color,
               xlab = "", ylab = "", zlab = "", axes = FALSE)
        box3d()
        par3d(windowRect = c(0, 70, 600, 670), userMatrix = user_matrix)
        
        # add_vector <- function(t, vec_length = 0.3, line_width = 2, col = "black")
        # {
        #     start_pos <- reconstruction[t,]
        #     vec <- reconstruction[t+1,] - start_pos
        #     vec <- vec * vec_length / sqrt(sum(vec^2))
        #     lines3d(rbind(start_pos, start_pos + vec), 
        #             lwd = line_width, col = col)
        #     return()
        # }
        # points3d(reconstruction[pred_target_t,1], 
        #          reconstruction[pred_target_t,2], 
        #          reconstruction[pred_target_t,3], 
        #          size = pt_size)
        # add_vector(pred_target_t)
        
        return()
    }

    load("model_2.Rdata")
    mu <- colMeans(data)
    sigma <- apply(data, 2, sd)
    orig_data <- normalize(data)
    noisy_data <- orig_data + rnorm(length(orig_data), sd = sqrt(0.1))
    
    load("model_2_fine.Rdata")
    fine_data <- data - matrix(rep(mu, each = NROW(data)), nrow = NROW(data))
    fine_data <- fine_data / matrix(rep(sigma, each = NROW(data)), nrow = NROW(data))
    n <- NROW(fine_data)
    block <- cbind(fine_data[,1], 
                   c(rep.int(NA, 16), fine_data[1:(n-16), 1]), 
                   c(rep.int(NA, 32), fine_data[1:(n-32), 1]), 
                   fine_data[,2], 
                   c(rep.int(NA, 16), fine_data[1:(n-16), 2]), 
                   c(rep.int(NA, 32), fine_data[1:(n-32), 2]), 
                   fine_data[,3], 
                   c(rep.int(NA, 16), fine_data[1:(n-16), 3]), 
                   c(rep.int(NA, 32), fine_data[1:(n-32), 3]))
    
    pdf("figure_1b0.pdf", width = 6, height = 6)
    par(mfrow = c(3,1), mar = c(1,1,0,0), lwd = 2)
    my_cols <- rainbow(9)
    pred_t <- (pred[1]*16-15):(pred[2]*16-15)
    for(i in c(1, 4, 7))
    {
        plot(-block[pred_t,i], col = my_cols[i], type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    }
    dev.off()
    
    # nn_1 <- make_panel(columns = c(2, 5, 6), user_matrix = matrix(c(0.992, -0.025, -0.127, 0, 
    #                                                           -0.019, 0.942, -0.335, 0, 
    #                                                           0.128, 0.225, 0.933, 0, 
    #                                                           0, 0, 0, 1), nrow = 4, byrow = TRUE))
    # rgl.postscript("figure_cc.pdf", "pdf")
    # 
    # nn_1 <- make_panel(columns = c(3, 7, 9))
    # rgl.postscript("figure_dd.pdf", "pdf")
    
    make_panel(columns = c(1, 3, 9), user_matrix = matrix(c(0.992, -0.025, -0.127, 0, 
                                                             -0.019, 0.942, -0.335, 0, 
                                                             0.128, 0.225, 0.933, 0, 
                                                             0, 0, 0, 1), nrow = 4, byrow = TRUE))
    rgl.postscript("figure_1bi.pdf", "pdf")
    make_panel(columns = c(2, 5, 6), user_matrix = matrix(c(0.491, 0.064, -0.869, 0, 
                                                             0.417, -0.893, 0.170, 0, 
                                                             -0.766, -0.446, -0.464, 0, 
                                                             0, 0, 0, 1), nrow = 4, byrow = TRUE))
    rgl.postscript("figure_1bii.pdf", "pdf")
    make_panel(columns = c(4, 8, 9), user_matrix = matrix(c(0.904, 0.428, -0.031, 0, 
                                                             -0.094, 0.267, 0.959, 0, 
                                                             0.418, -0.864, 0.282, 0, 
                                                             0, 0, 0, 1), nrow = 4, byrow = TRUE))
    rgl.postscript("figure_1biii.pdf", "pdf")
    
    
    user_matrix = matrix(c(0.406, -0.914, -0.008, 0,
                           0.485, 0.208, 0.849, 0,
                           -0.775, -0.349, 0.528, 0,
                           0, 0, 0, 1), nrow = 4, byrow = TRUE)
    pred_t <- (pred[1]*16-15):(pred[2]*16-15)
    plot3d(block[pred_t, c(1, 4, 7)], type = "l", col = "black",
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    box3d()
    n <- length(pred_t)
    x <- block[pred_t, 1]
    y <- block[pred_t, 4]
    z <- block[pred_t, 7]
    max_z <- max(z)
    min_z <- min(z)
    lines3d(rep.int(x[1], n), 
            y, 
            seq(from = max_z, to = max_z*2 - min_z, length.out = n), 
            col = "red")
    lines3d(c(x[1],x[1]), rep.int(y[1], y[1]), c(z[1], max_z), lty = 2)
    par3d(windowRect = c(0, 70, 800, 670))#, userMatrix = user_matrix)
    rgl.postscript("figure_1ai.pdf", "pdf")
    
    return()
}
