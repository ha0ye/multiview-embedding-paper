run_multiembed <- function(data_file = "model_1.Rdata", 
                           out_file = "model_1_output.Rdata", 
                           num_neighbors = "e+1", 
                           obs_err_var = 0.1, 
                           lib_start = 501, 
                           lib_sizes = c(25, 50, 75, 100, 200, 300))
{  
    start_time <- proc.time()
    message("Starting multiembed calculations", appendLF = FALSE)
    
    # load data, normalize, and apply noise
    load(data_file)
    data <- normalize(data)
    obs_data <- add_observation_error(data, obs_err_sd = sqrt(obs_err_var))
    
    # for each lib, run calculations
    temp <- lapply(lib_start, function(lib_start) {
        compute_multiembed_preds(obs_data, lib_sizes, lib_start = lib_start, 
                                 num_neighbors = num_neighbors)
    })
    multiembed_output <- lapply(temp, function(x) {x$output})
    multiembed_embeddings <- lapply(temp, function(x) {x$multivar_embeddings})
    
    # save output file
    save(multiembed_output, multiembed_embeddings, obs_data, data, file = out_file)
    message("done! (", round((proc.time()-start_time)[3], 3), " seconds elapsed.)")
    return()
}

compute_multiembed_preds <- function(data, lib_sizes, lib_start = 501, 
                                     num_neighbors = "e+1", 
                                     pred = c(2501, 3000), 
                                     E = 3, max_lag = 3, ranking = "rho", 
                                     max_num_embeddings_to_save = 150)
{
    message(".", appendLF = FALSE)
    num_vars <- NCOL(data)
    
    ### make block & generate embeddings
    block <- make_block(data, E, max_lag)
    embeddings_list <- make_embeddings(max_lag, E, num_vars)
    valid_embeddings <- apply(embeddings_list %% max_lag, 1, function (x) {1 %in% x})
    embeddings_list <- embeddings_list[valid_embeddings,]
    num_embeddings <- NROW(embeddings_list)

    # setup output data structures
    combined_preds <- data.frame()
    
    params <- expand.grid(lib_size = lib_sizes, target = 1:num_vars)
    out <- mclapply(1:NROW(params), function(i) {
        lib_size <- params$lib_size[i]
        target_col <- params$target[i]*max_lag - max_lag + 1
        lib = c(lib_start, lib_start-1+lib_size)
        
        # run embeddings
        in_stats <- block_lnlp(block, lib = lib, pred = lib, num_neighbors = num_neighbors, 
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
        k <- min(num_embeddings, max_num_embeddings_to_save)
        univar_embedding <- seq(from = target_col, to = target_col+E-1)
        multivar_embeddings <- embeddings_list[in_sample_ranking[1:k],] # sorted
        
        # compute univariate and multivariate results
        temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = num_neighbors, 
                           target_column = target_col, first_column_time = TRUE, 
                           columns = rbind(univar_embedding, multivar_embeddings[1,]), 
                           stats_only = FALSE, short_output = TRUE)
        univariate_results <- temp[[1]]$model_output
        multivariate_results <- temp[[2]]$model_output
        
        # compute predictions for multiembed
        temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = 1, 
                           target_column = target_col, first_column_time = TRUE, 
                           columns =  multivar_embeddings[1:k,], 
                           stats_only = FALSE, short_output = TRUE)
        
        # calculate multiembed predictions
        multiembed_pred <- do.call(cbind, lapply(1:k, 
                                                 function(j) temp[[j]]$model_output$pred))
        
        return(list(output = data.frame(lib_size = params$lib_size[i], 
                                        target = params$target[i], 
                                        obs = univariate_results$obs, 
                                        univariate_pred = univariate_results$pred, 
                                        multivariate_pred = multivariate_results$pred, 
                                        multiembed_pred = multiembed_pred), 
                    multivar_embeddings = multivar_embeddings))
    }, mc.cores = 7)
    
    return(out)
}

compute_multiembed_stats <- function(in_file = "model_1_output.Rdata", 
                                     out_file = "model_1_stats.Rdata", 
                                     num_to_merge = "sqrt")
{
    load(in_file)
    
    obs <- data[2502:3000,]
    stats <- do.call(rbind, lapply(seq_along(multiembed_output), function(i) {
        preds_list <- multiembed_output[[i]]
        df <- do.call(rbind, lapply(seq_along(preds_list), function(j) {
            preds <- preds_list[[j]]
            lib_size <- preds$lib_size[1]
            target <- preds$target[1]
            univar_stats <- compute_stats(obs[,target], preds$univariate_pred)
            multivar_stats <- compute_stats(obs[,target], preds$multivariate_pred)

            m <- grep("multiembed.+", names(preds))
            if(num_to_merge == "sqrt")
            {    
                k <- ceiling(sqrt(length(m)))
            } else {
                k <- num_to_merge
            }
            if(k > 1)
            {
                mve_pred <- rowMeans(preds[, m[1:k]])
            } else {
                mve_pred <- preds[, m[1]]
            }
            mve_stats <- compute_stats(obs[,target], mve_pred)
            return(data.frame(lib_size = lib_size, 
                              target = target, 
                              univar_rho = univar_stats$rho, 
                              univar_mae = univar_stats$mae, 
                              univar_rmse = univar_stats$rmse, 
                              multivar_rho = multivar_stats$rho, 
                              multivar_mae = multivar_stats$mae, 
                              multivar_rmse = multivar_stats$rmse, 
                              mve_rho = mve_stats$rho, 
                              mve_mae = mve_stats$mae, 
                              mve_rmse = mve_stats$rmse))
        }))
        df <- cbind(i, df)
        return(df)
    }))
    save(stats, file = out_file)
    return()
}

merge_stats_files <- function(obs_err = c(0, 0.1))
{
    stats_df <- data.frame()
    for(obs_err_var in obs_err)
    {    
        load(paste0("model_1_stats_", obs_err_var, ".Rdata"))
        temp <- cbind(model = 1, err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
        load(paste0("model_2_stats_", obs_err_var, ".Rdata"))
        temp <- cbind(model = 2, err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
        load(paste0("model_3_stats_", obs_err_var, ".Rdata"))
        temp <- cbind(model = 3, err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
        load(paste0("model_hw_stats_", obs_err_var, ".Rdata"))
        temp <- cbind(model = "hw", err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
    }
    
    stats_df$model <- as.factor(stats_df$model)
    stats_df <- stats_df %>%
        tidyr::gather(key, value, -model, -err, -i, -lib_size, -target) %>%
        tidyr::extract(key, c("method", "metric"), "([a-z]+)_([a-z]+)") %>%
        tidyr::spread(metric, value)
    stats_df$method <- as.factor(stats_df$method)
    stats_df$method <- plyr::revalue(stats_df$method, 
                                     c("multivar" = "multivariate", 
                                       "mve" = "multiview embedding", 
                                       "univar" = "univariate"))
    stats_df$method <- factor(stats_df$method, levels(stats_df$method)[c(3,1,2)])
    
    stats <- stats_df
    save(stats, file = "combined_stats.Rdata")
    return()
}

compute_multiembed_stats_k <- function(in_file = "model_1_output_0.1.Rdata", 
                                       out_file = "model_1_stats_0.1_k.Rdata")
{
    load(in_file)
    
    obs <- data[2502:3000,]
    stats <- do.call(rbind, lapply(seq_along(multiembed_output), function(i) {
        preds_list <- multiembed_output[[i]]
        df <- do.call(rbind, lapply(preds_list, function(preds) {
            lib_size <- preds$lib_size[1]
            target <- preds$target[1]
            
            m <- grep("multiembed.+", names(preds))
            mve_stats <- do.call(rbind, lapply(seq_along(m), function(k) {
                if(k > 1)
                {
                    mve_pred <- rowMeans(preds[, m[1:k]])
                } else {
                    mve_pred <- preds[, m[1]]
                }
                return(compute_stats(obs[,target], mve_pred))
            }))
            return(cbind(lib_size = lib_size, 
                         target = target, 
                         k = seq_along(m), 
                         mve_stats))
        }))
        df <- cbind(i, df)
        return(df)
    }))
    save(stats, file = out_file)
    return()
}

merge_stats_files_k <- function(obs_err = 0.1)
{
    stats_df <- data.frame()
    for(obs_err_var in obs_err)
    {    
        load(paste0("model_1_stats_", obs_err_var, "_k.Rdata"))
        temp <- cbind(model = 1, err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
        load(paste0("model_2_stats_", obs_err_var, "_k.Rdata"))
        temp <- cbind(model = 2, err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
        load(paste0("model_3_stats_", obs_err_var, "_k.Rdata"))
        temp <- cbind(model = 3, err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
        load(paste0("model_hw_stats_", obs_err_var, "_k.Rdata"))
        temp <- cbind(model = "hw", err = obs_err_var, stats)
        stats_df <- rbind(stats_df, temp)
    }
    
    stats_df$model <- as.factor(stats_df$model)
    stats_df <- stats_df %>%
        tidyr::gather("metric", value, c(mae, rho, rmse))
    stats_df$metric <- factor(stats_df$metric)
    stats_df$metric <- factor(stats_df$metric, levels(stats_df$metric)[c(2,1,3)])
    
    stats <- stats_df
    save(stats, file = "combined_stats_k.Rdata")
}


compute_model_2_forecast_errors <- function(lib = c(551, 600), pred = c(2001, 3000),
                                            embedding_1 = 1:3, embedding_2 = 4:6, embedding_3 = 7:9, 
                                            obs_err_var = 0.1)
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
        best_embeddings <- list()
        
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
            best_embeddings[[j]] <- multivar_embeddings
        }
        return(list(forecasts = output_pred, 
                    best_embeddings = best_embeddings))
    }
    
    load("model_2.Rdata")
    orig_data <- normalize(data)
    data <- add_observation_error(orig_data, obs_err_sd = sqrt(obs_err_var))
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
    temp <- get_multiembed_forecasts()
    multi_forecasts <- temp$forecasts
    best_embeddings <- temp$best_embeddings
    multi_errors <- multi_forecasts[pred[1]:(pred[2]-1),] - data[(pred[1]+1):pred[2],]
    multi_err_dist <- sqrt(rowSums(multi_errors * multi_errors))
    
    save(mx_forecasts, mx_err_dist, my_forecasts, my_err_dist,
         mz_forecasts, mz_err_dist, multi_forecasts, multi_err_dist,
         pred, data, orig_data, best_embeddings, file = "pred_errors.Rdata")
}



make_embeddings <- function(max_lag, E, num_vars)
{
    return(t(combn(num_vars*max_lag, E, simplify = TRUE)))
}

make_block <- function(data, E = 3, max_lag = 3, t = NULL, lib = NULL, tau = 1)
{
    num_vars <- NCOL(data)
    num_rows <- NROW(data)
    block <- matrix(NA, nrow = num_rows, ncol = 1+num_vars*max_lag)
    if(is.null(t))
        block[, 1] <- 1:num_rows
    else
        block[, 1] <- t

    # add max_lag lags for each column in data
    col_index <- 2
    for (j in 1:num_vars)
    {
        ts <- data[,j]
        block[, col_index] <- ts
        col_index <- col_index + 1

        if(max_lag > 1)
        {
            for(i in 1:(max_lag-1))
            {
                ts <- c(rep_len(NA, tau), ts[1:(num_rows-tau)])
                if(!is.null(lib))
                {
                    for(k in 1:NROW(lib))
                    {
                        ts[lib[k,1]] <- NA
                    }
                }
                block[, col_index] <- ts
                col_index <- col_index + 1
            }
        }
    }

    return(block)
}

compute_stats <- function(obs, pred)
{
    # computes performance metrics for how well predictions match observations
    # obs = vector of observations
    # pred = vector of prediction

    N = sum(is.finite(obs) & is.finite(pred))
    rho = cor(obs, pred, use = "pairwise.complete.obs")
    mae = mean(abs(obs-pred), na.rm = TRUE)
    rmse = sqrt(mean((obs-pred)^2, na.rm = TRUE))
    return(data.frame(N = N, rho = rho, mae = mae, rmse = rmse))
}