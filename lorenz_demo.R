library(rgl)
library(rEDM)

make_lorenz_data <- function()
{
    run_lorenz_model <- function(num_points = 11000)
    {
        delta_t <- 0.01
        data_skip <- 1
        sigma = 10
        beta = 8/3
        rho = 28
        
        # motion equations
        d <- function(x) {
            return(c(sigma * (x[2] - x[1]), 
                     x[1] * (rho - x[3]) - x[2], 
                     x[1]*x[2] - beta*x[3]))
        }
        
        x <- matrix(20, nrow = num_points, ncol = 3)
        
        for(i in 1:(num_points-1))
        {
            xx <- x[i,]
            for(k in 1:data_skip)
            {
                xx <- xx + delta_t * d(xx)
            }
            x[i+1,] <- xx
        }
        return(x)
    }
    
    # generate data
    data <- run_lorenz_model()
    normed_data <- normalize(data)
    obs_data <- data + matrix(rnorm(NROW(data) * NCOL(data), sd = 0.1), 
                              nrow = NROW(data), ncol = NCOL(data))
    
    save(normed_data, obs_data, data, file = "model_lorenz.Rdata")
    return()
}

plot_lorenz <- function(start_time = 500, lib_size = 50)
{
    t <- (start_time+1):(start_time+lib_size)
    load("model_lorenz.Rdata")
    load("model_1.Rdata")
    plot3d(obs_data[t,1], obs_data[t,2], obs_data[t,3], type = "p")
    
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
            multivar_embeddings <- embeddings_list[in_sample_ranking[1:16],]
            
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

make_error_plot <- function()
{
    load("pred_errors.Rdata")
    t <- (pred[1]+1):pred[2]
    pred_t <- pred[1]:(pred[2]-1)
    
    user_matrix <- matrix(c(0.6, -0.8, 0, 0, 
                            0.25, 0.2, 1, 0, 
                            -0.75, -0.55, 0.3, 0, 
                            0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(userMatrix = user_matrix, windowRect = c(0, 70, 800, 870))
    
    
    num_colors <- 30
    brks <- seq(length.out = num_colors, from = 0, 
                to = 3)#max(c(mx_err_dist, my_err_dist, mz_err_dist, multi_err_dist)))
    
    # plot M with Mx errors
    mx_pal <- colorRampPalette(c("red", "red"))
    mx_col <- mx_pal(num_colors)[as.numeric(cut(mx_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(mx_forecasts[pred_t,1], mx_forecasts[pred_t,2], mx_forecasts[pred_t,3], size = 3, col = mx_col)
    box3d()
    rgl.postscript("figure_2a.pdf", "pdf")
    
    # plot M with My errors
    my_pal <- colorRampPalette(c("green", "white"))
    my_col <- my_pal(num_colors)[as.numeric(cut(my_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC")
    points3d(my_forecasts[pred_t,1], my_forecasts[pred_t,2], my_forecasts[pred_t,3], size = 3, col = my_col)
    
    # plot M with Mz errors
    mz_pal <- colorRampPalette(c("blue", "blue"))
    mz_col <- mz_pal(num_colors)[as.numeric(cut(mz_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(mz_forecasts[pred_t,1], mz_forecasts[pred_t,2], mz_forecasts[pred_t,3], size = 3, col = mz_col)
    box3d()
    rgl.postscript("figure_2b.pdf", "pdf")
    
    # plot M with multiembed errors
    multi_pal <- colorRampPalette(c("magenta", "magenta"))
    multi_col <- multi_pal(num_colors)[as.numeric(cut(multi_err_dist, breaks = brks))]
    plot3d(orig_data[t,1], orig_data[t,2], orig_data[t, 3], type = "l", col = "#CCCCCC", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(multi_forecasts[pred_t,1], multi_forecasts[pred_t,2], multi_forecasts[pred_t,3], size = 3, col = multi_col)
    box3d()
    rgl.postscript("figure_2c.pdf", "pdf")
    
}

make_lorenz_plot <- function(t = 50:2000, x = c(1,0), y = c(2,0), z = c(3,0), 
                             x_label = "x", y_label = "y", z_label = "z", 
                             plot_file = NULL)
{
    load("model_lorenz.Rdata")
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

if(FALSE)
{
    make_lorenz_data()
    
    compute_forecast_errors()
    load("pred_errors.Rdata")
    boxplot(list(mx_err_dist, my_err_dist, mz_err_dist, multi_err_dist))
    
    make_error_plot()
    
    # plot native lorenz attractor
    make_lorenz_plot(plot_file = "lorenz_native.pdf")
    
    # plot univariate reconstruction (x?)
    make_lorenz_plot(y = c(1,8), z = c(1,16), plot_file = "lorenz_x.pdf")
    
    # plot univariate reconstruction (z?)
    make_lorenz_plot(x = c(3,0), y = c(3,8), z = c(3,16), plot_file = "lorenz_z.pdf")
}

