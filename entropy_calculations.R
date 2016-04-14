library(ggplot2)

entropy_hist <- function(x)
{
    I <- function(h)
    {
        b <- c(seq(from = mid, to = r[1]-h, by = -h), seq(from = mid, to = r[2]+h, by = h))
        N <- hist(x, breaks = sort(unique(b)), plot = FALSE)$counts
        I_hat <- 1/n * sum(N * log(N), na.rm = TRUE) - log(n * h)
        return(I_hat - sum(N > 0)/n)
    }
    
    n <- length(x)
    r <- range(x)
    iqr <- quantile(x, probs = c(0.25, 0.75))
    mid <- (iqr[1]+iqr[2])/2
    h_to_test <- seq(from = (iqr[2]-iqr[1])/sqrt(n), to = (r[2]-r[1])/4, length.out = 100)
    I_under <- sapply(h_to_test, I)
    I_meds <- runmed(I_under, 7)
    h <- h_to_test[which.max(I_meds)]
    
    #     h <- optim(2*(iqr[2]-iqr[1])/(n^(1/3)), f, 
    #                method = "L-BFGS-B", control = list(fnscale = -1), 
    #                lower = (r[2]-r[1])/(sqrt(n)), upper = r[2]-r[1])
    
    b <- c(seq(from = mid, to = r[1]-h, by = -h), seq(from = mid, to = r[2]+h, by = h))
    N <- hist(x, breaks = sort(unique(b)), plot = FALSE)$counts
    I_hat <- 1/n * sum(N * log(N), na.rm = TRUE)
    return(-I_hat)
}

entropy_kernel <- function(x)
{
    I <- function(h)
    {
        K <- as.matrix(dt(d/h, df = 4))
        f <- rowSums(K)
        I_hat <- sum(log(f))/n - log(n-1) - log(h)
        return(I_hat)
    }
    
    n <- length(x)
    r <- range(x)
    iqr <- quantile(x, probs = c(0.25, 0.75))
    mid <- (iqr[1]+iqr[2])/2
    h_to_test <- seq(from = (iqr[2]-iqr[1])/sqrt(n), to = (r[2]-r[1])/4, length.out = 100)
    d <- dist(x)
    I_hat <- sapply(h_to_test, I)
    return(-max(I_hat))
}

compute_entropy_convergence <- function(data_file, save_file)
{
    f <- function(x, method = "hist")
    {
        if(method == "hist")
            return(entropy_hist(x))
        return(entropy_kernel(x))
    }
    load(data_file)
    num_vars <- NCOL(data)
    
    pred_sizes <- c(25, 50, 100, 200, 400, 800)
    num_samples <- 100
    rand_starts <- sample(1001:(11000-max(pred_sizes)), num_samples)
    
    start_time <- proc.time()
    message("Calculating entropy", appendLF = FALSE)
    df <- data.frame()
    for(i in 1:num_samples)
    {
        message(".", appendLF = FALSE)
        for(pred_size in pred_sizes)
        {
            segment_idx <- seq(from = rand_starts[i], length.out = pred_size)
            bootstrap_idx <- sample(1001:11000, size = pred_size, replace = TRUE)
            for(j in 1:num_vars)
            {
                for(m in c("hist", "kernel"))
                {
                    df <- rbind(df, 
                                data.frame(index = i, 
                                           variable = j, 
                                           pred_size = pred_size, 
                                           entropy = f(data[segment_idx, j], m), 
                                           entropy_method = m, 
                                           sample_method = "segment"), 
                                data.frame(index = i, 
                                           variable = j, 
                                           pred_size = pred_size, 
                                           entropy = f(data[bootstrap_idx, j], m), 
                                           entropy_method = m, 
                                           sample_method = "bootstrap")
                    )
                }
            }
        }
    }
    message("done! (elapsed time = ", round((proc.time()-start_time)[3], 2), " seconds)")
    
    save(df, file = save_file)
    return()
}

plot_entropy_results <- function(plot_file = NULL, width = 6, height = 6)
{
    data_files <- c("entropy_results_1.Rdata", "entropy_results_2.Rdata", "entropy_results_3.Rdata")
    my_plots <- lapply(data_files,
                       function(file) {
                           load(file)
                           erc <- summarySE(df, measurevar = "entropy",
                                            groupvars = c("pred_size", "entropy_method", 
                                                          "variable", "sample_method"))
                           my_plot <- ggplot(erc, 
                                             aes(x = pred_size, y = entropy, 
                                                 color = entropy_method, linetype = sample_method, 
                                                 ymin = entropy-sd, ymax = entropy+sd)) + 
                               geom_errorbar(width = 0.1) + 
                               geom_line() + 
                               geom_point() +
                               facet_grid(variable ~ .) + 
                               theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     axis.text = element_text(color = "black"), 
                                     panel.background = element_rect(color = "black", fill = NA)) + 
                               xlab("Number of Points") + ylab("Estimated Differential Entropy")
                           
                           return(my_plot)
                       })
    
    if(!is.null(plot_file))
        pdf(file = plot_file, width = width, height = height, useDingbats = FALSE)
    
    lapply(my_plots, print)
    
    if(!is.null(plot_file))
        dev.off()
    
    return()
}

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


