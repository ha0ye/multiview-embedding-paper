make_HW_univariate_plots <- function(plot_file = NULL)
{
    load("hw_model_downsampled.Rdata")
    
    if(!is.null(plot_file))
    {
        pdf(plot_file, width = 6, height = 6)
    }
    for(i in 1:5)
    {
        ts <- data[,i]
        lib <- c(500, 600)
        pred <- c(2501, 3000)
        
        layout(matrix(c(1,1,2,3), byrow = TRUE, nrow = 2))
        par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
        plot(ts[lib[1]:lib[2]], type = "l", xlab = "Time", ylab = paste("N_", i, sep = ""))
        
        simplex_out <- simplex(ts, lib = lib, pred = pred)
        plot(simplex_out$E, simplex_out$rho, type = "l", xlab = "E", ylab = "rho")
        
        best_E <- simplex_out$E[which.max(simplex_out$rho)]
        
        smap_out <- s_map(ts, lib = lib, pred = pred, E = best_E)
        plot(smap_out$theta, smap_out$rho, type = "l", xlab = "theta", ylab = "rho")
    }
    if(!is.null(plot_file))
    {
        dev.off()
    }
    return()
}

compute_pred_skill_and_merge <- function(base_pred_file = "multiembed_preds_1_", 
                               output_file = "merged_stats_1.Rdata")
{
    message("Computing", appendLF = FALSE)
    output_stats <- do.call(rbind, mclapply(1:100, mc.cores = 6, function(file_idx) {
        load(paste("multiembed_pred_files/", base_pred_file, file_idx, ".Rdata", sep = ""))
        message(".", appendLF = FALSE)
        output_stats <- do.call(rbind, lapply(multiembed_output, function(preds)
        {
            df <- data.frame()
            temp_stats <- compute_stats(preds$obs, preds$univariate_pred)
            df <- rbind(df, data.frame(model = "univariate", num_embeddings = 1, 
                                       rho = temp_stats$rho, mae = temp_stats$mae))
            temp_stats <- compute_stats(preds$obs, preds$multivariate_pred)
            df <- rbind(df, data.frame(model = "multivariate", num_embeddings = 1, 
                                       rho = temp_stats$rho, mae = temp_stats$mae))
            
            multiembed_cols <- grep("multiembed", names(preds))
            for(i in seq_along(multiembed_cols))
            {
                temp_preds <- if(i == 1) preds[, multiembed_cols[1]] else rowMeans(preds[,multiembed_cols[1:i]])
                temp_stats <- compute_stats(preds$obs, temp_preds)
                df <- rbind(df, data.frame(model = "multiview", num_embeddings = i, 
                                           rho = temp_stats$rho, mae = temp_stats$mae))
            }
            
            df$lib_size <- preds$lib_size[1]
            df$target <- preds$target[1]
            return(df)
        }))
        output_stats$idx <- file_idx
        return(output_stats)
    }))
    save(output_stats, file = output_file)
    message("done!")
    return()
}

run_multiembed <- function(datafile, resultsfile, lib_start = 501, num_to_merge = 16, 
                           lib_sizes = c(25, 50, 75, 100, 200, 300))
{  
    start_time <- proc.time()
    message("Starting multiembed calculations", appendLF = FALSE)
    
    multiembed_output <- compute_multiembed(datafile, lib_sizes, lib_start = lib_start, num_to_merge = num_to_merge)
    combined_stats <- multiembed_output$combined_stats
    embeddings_list <- multiembed_output$embeddings_list
    num_vars <- multiembed_output$num_vars
    message("done! (", round((proc.time()-start_time)[3], 3), " seconds elapsed.)")
    
    # save results
    save(combined_stats, embeddings_list, num_vars, file = resultsfile)
    return()
}

run_multiembed_preds <- function(datafile, resultsfile, lib_start = 501, 
                           lib_sizes = c(25, 50, 75, 100, 200, 300))
{  
    start_time <- proc.time()
    message("Starting multiembed calculations", appendLF = FALSE)
    multiembed_output <- compute_multiembed_preds(datafile, lib_sizes, lib_start = lib_start)
    message("done! (", round((proc.time()-start_time)[3], 3), " seconds elapsed.)")
    
    # save results
    save(multiembed_output, file = resultsfile)
    return()
}


output_some_numbers <- function()
{
    combine_stats <- function()
    {
        stats <- data.frame()
        
        load("multiembed_all_results_1.Rdata")
        combined_stats$model <- 1
        stats <- rbind(stats, combined_stats)
        
        load("multiembed_all_results_2.Rdata")
        combined_stats$model <- 2
        stats <- rbind(stats, combined_stats)
        
        load("multiembed_all_results_3.Rdata")
        combined_stats$model <- 3
        stats <- rbind(stats, combined_stats)
        
        save(stats, file = "merged_multiembed_stats.Rdata")
        return()
    }
    
    get_multiembed_stats <- function(file)
    {
        
        load(file)
        lib_sizes <- sort(unique(combined_stats$lib))
        num_vars <- length(unique(combined_stats$target))
        
        for(j in 1:num_vars)
        {
            rows <- combined_stats$target == j
            rhos <- combined_stats[rows, c("lib", "univar_rho", "multivar_rho", "multiembed_rho")]
            maes <- combined_stats[rows, c("lib", "univar_mae", "multivar_mae", "multiembed_mae")]
            rmses <- combined_stats[rows, c("lib", "univar_rmse", "multivar_rmse", "multiembed_rmse")]
        }
        
        vals <- aggregate(values_to_plot, list(values_to_plot$lib), mean)
        vals <- vals[,3:5]
        return()
    }
    
    combine_stats()
    
    load("merged_multiembed_stats.Rdata")
    
    aggregated_stats <- aggregate(stats, list(model = stats$model, target = stats$target, lib = stats$lib), mean)
    aggregated_stats$delta_univar_rho <- aggregated_stats$multiembed_rho - aggregated_stats$univar_rho
    aggregated_stats$delta_multivar_rho <- aggregated_stats$multiembed_rho - aggregated_stats$multivar_rho
    aggregated_stats$delta_univar_mae <- aggregated_stats$multiembed_mae / aggregated_stats$univar_mae - 1
    aggregated_stats$delta_multivar_mae <- aggregated_stats$multiembed_mae / aggregated_stats$multivar_mae - 1
    aggregated_stats$delta_univar_rmse <- aggregated_stats$multiembed_rmse / aggregated_stats$univar_rmse - 1
    aggregated_stats$delta_multivar_rmse <- aggregated_stats$multiembed_rmse / aggregated_stats$multivar_rmse - 1
    
    stats$delta_univar_rho <- stats$multiembed_rho - stats$univar_rho
    stats$delta_multivar_rho <- stats$multiembed_rho - stats$multivar_rho
    stats$delta_univar_mae <- stats$multiembed_mae / stats$univar_mae - 1
    stats$delta_multivar_mae <- stats$multiembed_mae / stats$multivar_mae - 1
    stats$delta_univar_rmse <- stats$multiembed_rmse / stats$univar_rmse - 1
    stats$delta_multivar_rmse <- stats$multiembed_rmse / stats$multivar_rmse - 1
    
    save(stats, aggregated_stats, file = "merged_multiembed_stats.Rdata")
    return()
}

compute_univariate <- function(data_file, save_file, 
                               lib_sizes, lib_start = 501, pred = c(2501, 3000))
{
    load(data_file)
    data <- normed_data
    
    results <- data.frame()
    for(i in 1:NCOL(data))
    {
        ts <- data[,i]
        for (lib_size in lib_sizes) # for each lib size
        {
            message(".", appendLF = FALSE)
            lib <- c(lib_start, lib_start + lib_size-1)
            
            lnlp_out <- simplex(ts, lib = lib, pred = pred, E = 1:10, 
                                stats_only = FALSE)
            
            stats <- do.call(rbind, lapply(lnlp_out, function(l) {
                obs <- l$model_output$obs
                pred <- l$model_output$pred
                return(compute_stats(obs, pred))
            }))
            stats$lib_size <- lib_size
            stats$var <- i
            results <- rbind(results, stats)
        }
    }
    save(results, file = save_file)
    return(results)
}

##
compute_multiembed_preds <- function(data_file, lib_sizes, lib_start = 501, pred = c(2501, 3000), 
                                     E = 3, max_lag = 3, ranking = "rho")
{
    load(data_file)
    data <- normed_data
    num_vars <- NCOL(data)
    
    ### make block & generate embeddings
    block <- make_block(data, E, max_lag)
    embeddings_list <- make_embeddings(max_lag, E, num_vars)
    valid_embeddings <- apply(embeddings_list %% max_lag, 1, function (x) {1 %in% x})
    embeddings_list <- embeddings_list[valid_embeddings,]
    num_embeddings <- NROW(embeddings_list)
    block <- block[1:max(pred),]
    
    # setup output data structures
    combined_stats <- data.frame()
    
    params <- expand.grid(lib_size = lib_sizes, target = 1:num_vars)
    return(lapply(1:NROW(params), function(i) {
        message(".", appendLF = FALSE)
        lib_size <- params$lib_size[i]
        target_col <- params$target[i]*max_lag - max_lag + 1
        lib = c(lib_start, lib_start-1+lib_size)
        
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
        multivar_embeddings <- embeddings_list[in_sample_ranking[1:num_embeddings],] # sorted
        
        # compute univariate and multivariate results
        temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = "e+1", 
                           target_column = target_col, first_column_time = TRUE, 
                           columns = rbind(univar_embedding, multivar_embeddings[1,]), 
                           stats_only = FALSE, short_output = TRUE)
        univariate_results <- temp[[1]]$model_output
        multivariate_results <- temp[[2]]$model_output

        # compute predictions for multiembed
        temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = 1, 
                           target_column = target_col, first_column_time = TRUE, 
                           columns =  multivar_embeddings[1:num_embeddings,], 
                           stats_only = FALSE, short_output = TRUE)
        
        # calculate multiembed predictions
        multiembed_pred <- do.call(cbind, lapply(1:num_embeddings, 
                                                 function(j) temp[[j]]$model_output$pred))
        
        return(data.frame(lib_size = params$lib_size[i], 
                          target = params$target[i], 
                          obs = univariate_results$obs, 
                          univariate_pred = univariate_results$pred, 
                          multivariate_pred = multivariate_results$pred, 
                          multiembed_pred = multiembed_pred))
    }))
}

## full calculation of stats
compute_multiembed <- function(data_file, lib_sizes, lib_start = 501, pred = c(2501, 3000), 
                               E = 3, max_lag = 3, num_to_merge = 16, 
                               ranking = "rho")
{
    load(data_file)
    data <- normed_data
    num_vars <- NCOL(data)
    
    ### make block & generate embeddings
    block <- make_block(data, E, max_lag)
    embeddings_list <- make_embeddings(max_lag, E, num_vars)
    valid_embeddings <- apply(embeddings_list %% max_lag, 1, function (x) {1 %in% x})
    embeddings_list <- embeddings_list[valid_embeddings,]
    num_to_merge <- min(NROW(embeddings_list), num_to_merge)
    
    # setup output data structures
    combined_stats <- data.frame()
    
    params <- expand.grid(lib_size = lib_sizes, target = 1:num_vars)
    output <- do.call(rbind, lapply(1:NROW(params), function(i) {
        message(".", appendLF = FALSE)
        lib_size <- params$lib_size[i]
        target_col <- params$target[i]*max_lag - max_lag + 1
        lib = c(lib_start, lib_start-1+lib_size)
        
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
        univariate_results <- temp[[1]]$model_output
        univariate_stats <- compute_entropy_stats(univariate_results$obs, univariate_results$pred)
        multivariate_results <- temp[[2]]$model_output
        multivariate_stats <- compute_entropy_stats(multivariate_results$obs, multivariate_results$pred)
        
        # compute predictions for multiembed
        temp <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = 1, 
                           target_column = target_col, first_column_time = TRUE, 
                           columns =  multivar_embeddings[1:num_to_merge,], 
                           stats_only = FALSE)
        
        # calculate multiembed predictions
        top_embeddings_preds <- do.call(cbind, lapply(1:num_to_merge, 
                                                      function(j) temp[[j]]$model_output$pred))
        multiembed_pred <- if(num_to_merge > 1) rowMeans(top_embeddings_preds) else top_embeddings_preds
        multiembed_stats <- compute_entropy_stats(temp[[1]]$model_output$obs, multiembed_pred)
        
        return(data.frame(univar_rho = univariate_stats$rho, 
               univar_mae = univariate_stats$mae, 
               univar_rmse = univariate_stats$rmse, 
               univar_pp = univariate_stats$predictive_power, 
               multivar_rho = multivariate_stats$rho, 
               multivar_mae = multivariate_stats$mae, 
               multivar_rmse = multivariate_stats$rmse, 
               multivar_pp = multivariate_stats$predictive_power, 
               multiembed_rho = multiembed_stats$rho, 
               multiembed_mae = multiembed_stats$mae, 
               multiembed_rmse = multiembed_stats$rmse, 
               multiembed_pp = multiembed_stats$predictive_power))
    }))
    combined_stats <- cbind(params, output)
    
    # return results
    return(list(combined_stats = combined_stats, 
                embeddings_list = embeddings_list, 
                num_vars = num_vars))
}

combine_multiembed_results <- function(model_index, num_trials = 100, 
                                       my_path = "multiembed_temp/")
{
    all_stats <- data.frame()
    for(trial in 1:num_trials)
    {
        file <- paste(my_path, "multiembed_results_", model_index, "_", trial, ".Rdata", sep = "")
        result <- tryCatch({
            load(file)
            all_stats <- rbind(all_stats, combined_stats)
        }, warning = function(w) {
            stop("elevating warning to error")
        }, error = function(e) {
            cat("stopped combining multiembed results at trial #", trial, "\n", sep = "")
            return(NULL)
        }) # END tryCatch
        if(is.null(result))
            break
    }
    combined_stats <- all_stats
    save(combined_stats, file = paste("multiembed_all_results_", model_index, ".Rdata", sep = ""))
    
    return()
}

make_embeddings <- function(max_lag, E, num_vars)
{
    return(t(combn(num_vars*max_lag, E, simplify = TRUE)))
}

row.matches <- function(vec, mat)
{
    i <- seq(NROW(mat))
    j <- 0
    while(length(i) && (j <- j + 1) <= NCOL(mat))
        i <- i[mat[i,j] == vec[j]]
    i
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
    
    # add max_lag lags for each column in normed_data
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

join_preds <- function(preds)
{
    time <- sort(unique(do.call(c, sapply(preds, '[[', "time"))))
    output <- data.frame(time)
    output$obs <- rep_len(NA, length(time))
    
    for(i in seq_along(preds))
    {
        index <- match(preds[[i]]$time, time)
        
        # combine obs
        if(i > 1 && any(abs(output$obs[index] - preds[[i]]$obs) > 0.0001, na.rm = TRUE))
            warning("some observed points do not line up.")
        output$obs[index] <- preds[[i]]$obs
        
        # combine preds
        for(j in 3:length(preds[[i]]))
        {
            col_name <- paste(names(preds)[i], ".", j-2, sep = "")
            output[, col_name] <- rep_len(NaN, length(time))
            output[index, col_name] <- preds[[i]][,j]
        }
    }
    
    return(output)
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

compute_entropy_stats <- function(obs, pred)
{
    # computes performance metrics for how well predictions match observations
    # obs = vector of observations
    # pred = vector of prediction
    
    err <- pred - obs
    idx <- is.finite(obs) & is.finite(pred)
    prior_entropy <- entropy_kernel(obs[idx])
    posterior_entropy <- entropy_kernel(err[idx])
    
    return(data.frame(N = sum(idx), rho = cor(obs[idx], pred[idx]), 
                      mae = mean(abs(err[idx])), rmse = sqrt(mean(err[idx]^2)), 
                      predictive_power = 1 - exp(posterior_entropy - prior_entropy)))
}