run_in_out_analysis <- function()
{
    set.seed(1234)
    lib_locations <- floor(runif(100, min = 101, max = 1101))
    
    cat("running in-out comparison for model 1... ", sep = "")
    start_time <- proc.time()
    model_1_in_out <- lapply(c(1:100), function(i) {
        cat(".", sep = "")
        do_in_out_comparison("model_1.Rdata", lib_sizes = c(25, 50, 100), 
                             lib_start = lib_locations[i])
    })
    save(model_1_in_out, file = "multiembed_in_out_1.Rdata")
    elapsed_time <- proc.time() - start_time
    cat("(", elapsed_time[3], " sec.)\n", sep = "")
    
    cat("running in-out comparison for model 2... ", sep = "")
    start_time <- proc.time()
    model_2_in_out <- lapply(c(1:100), function(i) {
        cat(".", sep = "")
        do_in_out_comparison("model_2.Rdata", lib_sizes = c(25, 50, 100), 
                             lib_start = lib_locations[i])
    })
    save(model_2_in_out, file = "multiembed_in_out_2.Rdata")
    elapsed_time <- proc.time() - start_time
    cat("(", elapsed_time[3], " sec.)\n", sep = "")
    
    cat("running in-out comparison for model 3... ", sep = "")
    start_time <- proc.time()
    model_3_in_out <- lapply(c(1:100), function(i) {
        cat(".", sep = "")
        do_in_out_comparison("model_3.Rdata", lib_sizes = c(25, 50, 100), 
                             lib_start = lib_locations[i])
    })
    save(model_3_in_out, file = "multiembed_in_out_3.Rdata")
    elapsed_time <- proc.time() - start_time
    cat("(", elapsed_time[3], " sec.)\n", sep = "")
    return()
}

make_in_out_plots <- function()
{
    load("multiembed_in_out_1.Rdata")
    make_multiembed_in_out_plot(data = model_1_in_out, file = "fig_s1.pdf")
    
    load("multiembed_in_out_2.Rdata")
    make_multiembed_in_out_plot(data = model_2_in_out, file = "fig_s2.pdf")
    
    load("multiembed_in_out_3.Rdata")
    make_multiembed_in_out_plot(data = model_3_in_out, file = "fig_s3.pdf", 
                                x_lab = "predicting larvae", 
                                y_lab = "predicting pupae", 
                                z_lab = "predicting adults")
    return()
}

make_time_lag_rho_plots <- function()
{
    load("multiembed_in_out_1.Rdata")
    make_time_lag_plot(data = model_1_in_out, file = "fig_s4.pdf")
    
    load("multiembed_in_out_2.Rdata")
    make_time_lag_plot(data = model_2_in_out, file = "fig_s5.pdf")
    
    load("multiembed_in_out_3.Rdata")
    make_time_lag_plot(data = model_3_in_out, file = "fig_s6.pdf", 
                       x_lab = "predicting larvae", 
                       y_lab = "predicting pupae", 
                       z_lab = "predicting adults")
    
    return()
}

do_in_out_comparison <- function(file, lib_sizes, lib_start = 101, pred = c(1501, 2000), 
                                 E = 3, max_lag = 3)
{
    load(file)
    data <- normed_data
    num_vars <- NCOL(data)
    
    ### make block & generate embeddings
    block <- make_block(data, E, max_lag)
    embeddings_list <- make_embeddings(max_lag, E, num_vars)
    valid_embeddings <- apply(embeddings_list %% max_lag, 1, function (x) {1 %in% x})
    embeddings_list <- embeddings_list[valid_embeddings,]
    
    # setup output data structures
    params <- expand.grid(lib_size = lib_sizes, target = 1:num_vars)
    output <- do.call(rbind, lapply(1:NROW(params), function(i) {
        lib_size <- params$lib_size[i]
        target_col <- params$target[i]*max_lag - max_lag + 1
        lib = c(lib_start, lib_start-1+lib_size)
        
        # run embeddings
        in_stats <- block_lnlp(block, lib = lib, pred = lib, num_neighbors = "e+1", 
                               target_column = target_col, columns = embeddings_list, 
                               first_column_time = TRUE, silent = TRUE)
        out_stats <- block_lnlp(block, lib = lib, pred = pred, num_neighbors = "e+1", 
                                target_column = target_col, columns = embeddings_list, 
                                first_column_time = TRUE, silent = TRUE)
        
        stats <- data.frame(embedding = in_stats$embedding, lib_size = lib_size, 
                            target = target_col, 
                            in_rho = in_stats$rho, in_mae = in_stats$mae, in_rmse = in_stats$rmse, 
                            out_rho = out_stats$rho, out_mae = out_stats$mae, out_rmse = out_stats$rmse)
        
        return(stats)
    }))
    
    # return results
    return(output)
}

make_multiembed_in_out_plot <- function(data, file = NULL, width = 6, height = 5, 
                                        x_lab = "predicting x", 
                                        y_lab = "predicting y", 
                                        z_lab = "predicting z")
{
    data <- do.call(rbind, data)
    model_summary <- ddply(data, .(embedding, lib_size, target), 
                           summarize, in_50 = median(in_rho), 
                           out_50 = median(out_rho), 
                           in_10 = quantile(in_rho, 0.1), 
                           in_90 = quantile(in_rho, 0.9), 
                           out_10 = quantile(out_rho, 0.1), 
                           out_90 = quantile(out_rho, 0.9), 
                           in_25 = quantile(in_rho, 0.25), 
                           in_75 = quantile(in_rho, 0.75), 
                           out_25 = quantile(out_rho, 0.25), 
                           out_75 = quantile(out_rho, 0.75))
    new_labels <- function(var, value) 
    {
        if (var == "target")
        {
            value[value == 1] <- x_lab
            value[value == 4] <- y_lab
            value[value == 7] <- z_lab
        }
        else
        {
            value[value == 25] <- "library size = 25"
            value[value == 50] <- "library size = 50"
            value[value == 100] <- "library size = 100"
        }
        return(value)
    }
    #    my_plot <- ggplot(model_summary, aes(in_50, out_50)) + 
    #    geom_point(size = 1, shape = 3) + 
    my_plot <- ggplot(model_summary, aes(in_50, out_50)) + 
        geom_linerange(size = 0.25, aes(ymin = out_25, ymax = out_75)) + 
        geom_errorbarh(size = 0.25, aes(xmin = in_25, xmax = in_75, height = 0)) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2, color = "blue") + 
        facet_grid(lib_size ~ target, labeller = new_labels) + 
        xlab(expression("In-Sample Forecast Skill ("*rho*")")) + 
        ylab(expression("Out-of-Sample Forecast Skill ("*rho*")")) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              panel.background = element_rect(color = "black", fill = NA))
    
    if(!is.null(file))
        pdf(file = file, width = width, height = height)
    print(my_plot)
    if(!is.null(file))
        dev.off()
}

make_time_lag_plot <- function(data, file = NULL, width = 6, height = 4, 
                               x_lab = "predicting x", 
                               y_lab = "predicting y", 
                               z_lab = "predicting z")
{
    new_labels <- function(var, value) 
    {
        if (var == "target")
        {
            value[value == 1] <- x_lab
            value[value == 4] <- y_lab
            value[value == 7] <- z_lab
        }
        return(value)
    }
    data <- do.call(rbind, data)
    model_summary <- ddply(data, .(embedding, lib_size, target), 
                           summarize, in_rho = median(in_rho), out_rho = median(out_rho))
    vars <- strsplit(as.character(model_summary$embedding), ", ")
    model_summary$total_lag <- sapply(vars, function(embedding) {
        cols <- as.numeric(embedding)
        sum((cols-1) %% 3)
    })
    
    model_summary <- model_summary[model_summary$lib_size == 100,]
    my_plot <- ggplot(model_summary, aes(total_lag, out_rho)) + 
        geom_point(size = 1, shape = 3) + 
        geom_smooth(method = lm, se = FALSE) + 
        facet_grid(. ~ target, labeller = new_labels) + 
        xlab(expression("Total Time Lag of Reconstruction")) + 
        ylab(expression("Out-of-Sample Forecast Skill ("*rho*")")) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"), 
              panel.background = element_rect(color = "black", fill = NA))
    
    if(!is.null(file))
        pdf(file = file, width = width, height = height)
    print(my_plot)
    if(!is.null(file))
        dev.off()
    return()
}