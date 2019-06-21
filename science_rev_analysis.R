obs_err <- c(0.0, 0.10, 0.20, 0.30)
folders <- c("~/Desktop/bak/", "~/Desktop/bak2/", "~/Desktop/bak3/", 
             "~/Desktop/bak4/", "~/Desktop/bak5/")

# compute stats files
if(FALSE)
{
    for(folder_dir in folders)
    {
        for(obs_err_var in obs_err)
        {
            pred_file <- paste0(folder_dir, "model_4_output_", obs_err_var, ".Rdata")
            stat_file <- paste0(folder_dir, "model_4_stats_", obs_err_var, ".Rdata")
            compute_multiembed_stats(in_file = pred_file, 
                                     out_file = stat_file)
        }
    }
}

# merge stats files
if(FALSE)
{
    stats_df <- data.frame()
    cum_idx <- 0
    for(folder_dir in folders)
    {
        temp_df <- data.frame()
        for(obs_err_var in obs_err)
        {    
            load(paste0(folder_dir, "model_4_stats_", obs_err_var, ".Rdata"))
            temp <- cbind(model = 4, err = obs_err_var, stats)
            temp_df <- rbind(temp_df, temp)
        }
        temp_df$i <- temp_df$i + cum_idx
        cum_idx <- max(temp_df$i)
        stats_df <- rbind(stats_df, temp_df)        
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
    save(stats, file = "combined_stats_rev.Rdata")
}

# compute stats files (adjusting for k)
if(FALSE)
{
    for(folder_dir in folders)
    {
        for(obs_err_var in 0.3)
        {
            pred_file <- paste0(folder_dir, "model_4_output_", obs_err_var, ".Rdata")
            stat_file <- paste0(folder_dir, "model_4_stats_", obs_err_var, "_k.Rdata")
            compute_multiembed_stats_k(in_file = pred_file, 
                                     out_file = stat_file)
        }
    }
}

# merge stats files (adjusting for k)
if(FALSE)
{
    stats_df <- data.frame()
    cum_idx <- 0
    for(folder_dir in folders)
    {
        temp_df <- data.frame()
        for(obs_err_var in 0.3)
        {    
            load(paste0(folder_dir, "model_4_stats_", obs_err_var, "_k.Rdata"))
            temp <- cbind(model = 4, err = obs_err_var, stats)
            temp_df <- rbind(temp_df, temp)
        }
        temp_df$i <- temp_df$i + cum_idx
        cum_idx <- max(temp_df$i)
        stats_df <- rbind(stats_df, temp_df)
    }
    stats_df$model <- as.factor(stats_df$model)
    stats_df <- stats_df %>%
        tidyr::gather("metric", value, c(mae, rho, rmse))
    stats_df$metric <- factor(stats_df$metric)
    stats_df$metric <- factor(stats_df$metric, levels(stats_df$metric)[c(2,1,3)])
    
    stats <- stats_df
    save(stats, file = "combined_stats_rev_k.Rdata")
}

# produce plots for 12-species system
if(TRUE)
{
    my_plots <- lapply(obs_err, plot_metric_vs_lib_size, stats_file = "combined_stats_rev.Rdata")
    pdf("mve_plot_rho_vs_lib_size_rev.pdf", width = 10, height = 8)
    print(my_plots)
    dev.off()
    
    my_plots <- plot_perf_vs_lib_size(system = 4, obs_err_var = 0.3, stats_file = "combined_stats_rev.Rdata")
    pdf("mve_plot_perf_vs_lib_size_rev.pdf", width = 8, height = 32)
    print(my_plots)
    dev.off()
    
    my_plots <- plot_metric_vs_noise(system = 4,stats_file = "combined_stats_rev.Rdata", 
                                     lib_sizes = c(25, 50, 100))
    pdf("mve_plot_mae_vs_noise_rev.pdf", width = 8, height = 32)
    print(my_plots)
    dev.off()
    
    my_plots <- plot_perf_vs_k(system = 4, stats_file = "combined_stats_rev_k.Rdata", 
                               lib_sizes = c(25, 50, 100))
    pdf("mve_plot_rho_vs_k_rev.pdf", width = 7, height = 28)
    print(my_plots)
    dev.off()
}

# compare 4th-root-transformed predictions with de-transformed predictions
if(FALSE)
{
    make_plot <- function(df, raw_vals)
    {
        make_panel <- function(x, y, plot_limits = range(x, y, na.rm = TRUE))
        {
            plot(x, y, xlim = plot_limits, ylim = plot_limits)
            abline(a = 0, b = 1, lty = 2, col = "blue")
            legend("topleft", border = NA, 
                   legend = c(paste0("cor = ", round(c(cor(x, y, use = "pairwise")), 3)), 
                              paste0("mae = ", round(sum(abs(x - y), na.rm = TRUE), 3))))
            return()
        }
        df$multiembed_pred <- rowMeans(df[, 5:12])
        
        # plot pred vs. obs (transformed)
        par(mfrow = c(2,3), mar = c(4,4,1,1), pty = "s")
        plot_limits <- range(df$obs, df$univariate_pred, 
                             df$multivariate_pred, df$multiembed_pred, na.rm = TRUE)
        make_panel(df$obs, df$univariate_pred, plot_limits)
        make_panel(df$obs, df$multivariate_pred, plot_limits)
        make_panel(df$obs, df$multiembed_pred, plot_limits)
        
        # denorm and detransform
        denorm_lm <- lm(raw_vals[2:725] ~ df$obs[1:724])
        m <- denorm_lm$coefficients[2]
        b <- denorm_lm$coefficients[1]
        for(j in 2:NCOL(df))
        {
            df[,j] <- (df[,j] * m + b)^4
        }
        
        # plot pred vs. obs (raw)
        plot_limits <- range(df$obs, df$univariate_pred, 
                             df$multivariate_pred, df$multiembed_pred, na.rm = TRUE)
        make_panel(df$obs, df$univariate_pred, plot_limits)
        make_panel(df$obs, df$multivariate_pred, plot_limits)
        make_panel(df$obs, df$multiembed_pred, plot_limits)
        return()
    }
    
    load("base_results_apr_2016/mesocosm_data.Rdata")
    load("base_results_apr_2016/mesocosm_multiembed_results.Rdata")
    load("base_results_apr_2016/mesocosm_multiembed_stats.Rdata")
    
    make_plot(calanoid_preds, obs_data$calanoid_copepods)
    make_plot(rotifer_preds, obs_data$rotifers)
}

# de-transform and re-compute stats
if(FALSE)
{
    compute_raw_stats <- function(df, raw_vals)
    {
        df$multiembed_pred <- rowMeans(df[, 5:12])
        
        # denorm and detransform
        denorm_lm <- lm(raw_vals[2:725] ~ df$obs[1:724])
        m <- denorm_lm$coefficients[2]
        b <- denorm_lm$coefficients[1]
        for(j in 2:NCOL(df))
        {
            df[,j] <- (df[,j] * m + b)^4
        }
        
        output <- lapply(list(df$univariate_pred, 
                              df$multivariate_pred, 
                              df$multiembed_pred), function(preds)
                              {
                                  temp <- compute_stats(df$obs, preds)
                                  temp$N <- temp$num_pred
                                  temp <- temp[, -1]
                                  temp$predictive_power <- NA
                                  return(temp)
                              })
        names(output) <- c("univariate_stats", 
                           "multivariate_stats", 
                           "multiembed_stats")
        return(output)
    }
    
    load("base_results_apr_2016/mesocosm_data.Rdata")
    load("base_results_apr_2016/mesocosm_multiembed_results.Rdata")

    calanoid_stats <- compute_raw_stats(calanoid_preds, obs_data$calanoid_copepods)
    rotifer_stats <- compute_raw_stats(rotifer_preds, obs_data$rotifers)
    save(calanoid_stats, rotifer_stats, file = "base_results_apr_2016/mesocosm_multiembed_raw_stats.Rdata")
}

# make plots for de-transformed stats
if(FALSE)
{
    plot_mesocosm_multiembed(stats_file = "base_results_apr_2016/mesocosm_multiembed_raw_stats.Rdata", out_file = "figure_4c.pdf")
    plot_mesocosm_multiembed_full(stats_file = "base_results_apr_2016/mesocosm_multiembed_raw_stats.Rdata", out_file = "ed_figure_4.pdf")
}

# raw mesocosm time series plots
if(FALSE)
{
    load("base_results_apr_2016/mesocosm_data.Rdata")
    obs_data <- obs_data[order(obs_data$day),]
    for(j in 2:NCOL(obs_data))
        obs_data[, j] <- obs_data[, j] ^ 4

    pdf("mesocosm_ts.pdf", width = 6, height = 9)
    par(mfrow = c(4,1), mgp = c(2.5,1,0), 
        mar = c(0.5,4,0.5,0.5), oma = c(4,0,0,0))
    plot(obs_data$day, obs_data$calanoid_copepods, type = "l", 
         xlab = "", ylab = "Calanoid Copepods", xaxt = "n")
    axis(1, labels = FALSE)
    plot(obs_data$day, obs_data$rotifers, type = "l", 
         xlab = "", ylab = "Rotifers", xaxt = "n")
    axis(1, labels = FALSE)
    plot(obs_data$day, obs_data$picophytoplankton, type = "l", 
         xlab = "", ylab = "Picophytoplankton", xaxt = "n")
    axis(1, labels = FALSE)
    plot(obs_data$day, obs_data$nanophytoplankton, type = "l", 
         xlab = "Day", ylab = "Nanophytoplankton")
    mtext("Day", 1, line = 2, outer = TRUE)
    dev.off()
}

rho_comp <- function(x1, x2, y)
    # computes p-value for cor(x1, y) > cor(x2, y) using 
    # t-test with df = length(y) - 2
{
    if(identical(x1, x2))
        return(0.5)
    n <- sum(is.finite(x1) & is.finite(x2) & is.finite(y))
    x1y <- cor(x1, y, use = "pairwise")
    x2y <- cor(x2, y, use = "pairwise")
    err <- TWOpov_err(as.matrix(cbind(x1, x2)), y)
    p_value <- 1 - pt((x1y - x2y) / err, df = n-2, lower.tail = TRUE)
    return(data.frame(df = n-2, statistic = (x1y - x2y) / err, p.value = p_value))
}

TWOpov<-function(x,y,alpha=.05,CN=F)
{
    #
    # Comparing two dependent correlations: Overlapping case
    #
    # x is assumed to be a matrix with 2 columns
    #
    #  Compare correlation of x[,1] with y to x[,2] with y
    #
    if(!is.matrix(x))stop("x should be a matrix")
    if(ncol(x)!=2)stop("x should be a matrix with two columns")
    xy=elimna(cbind(x,y))
    x1=xy[,1]
    x2=xy[,2]
    y=xy[,3]
    r12=cor(x1,y)
    r13=cor(x2,y)
    r23=cor(x1,x2)
    ci12=pcorhc4(x1,y,alpha=alpha,CN=CN)$ci
    ci13=pcorhc4(x2,y,alpha=alpha,CN=CN)$ci
    corhat=((r23-.5*r12*r13)*(1-r12^2-r13^2-r23^2)+r23^3)/((1-r12^2)*(1-r13^2))
    term1=2*corhat*(r12-ci12[1])*(ci13[2]-r13)
    term2=2*corhat*(r12-ci12[2])*(ci13[1]-r13)
    L=r12-r13-sqrt((r12-ci12[1])^2+(ci13[2]-r13)^2-term1)
    U=r12-r13+sqrt((r12-ci12[2])^2+(ci13[1]-r13)^2-term2)
    c(L,U)
}

TWOpov_err <- function(x,y,CN=F)
{
    #
    # Comparing two dependent correlations: Overlapping case
    #
    # x is assumed to be a matrix with 2 columns
    #
    #  Compare correlation of x[,1] with y to x[,2] with y
    #
    # returns p-value
    if(!is.matrix(x))stop("x should be a matrix")
    if(ncol(x)!=2)stop("x should be a matrix with two columns")
    xy=elimna(cbind(x,y))
    x1=xy[,1]
    x2=xy[,2]
    y=xy[,3]
    r12=cor(x1,y)
    r13=cor(x2,y)
    r23=cor(x1,x2)
    err12 <- pcorhc4_err(x1,y,CN=CN)
    err13 <- pcorhc4_err(x2,y,CN=CN)
    corhat=((r23-.5*r12*r13)*(1-r12^2-r13^2-r23^2)+r23^3)/((1-r12^2)*(1-r13^2))
    err_correction_term = 2*corhat*(err12)*(err13)
    err_diff <- sqrt(err12^2 + err13^2 - err_correction_term)
    return(err_diff)
}

elimna<-function(m)
{
    #
    # remove any rows of data having missing values
    #
    if(is.null(dim(m)))m<-as.matrix(m)
    ikeep<-c(1:nrow(m))
    for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
    elimna<-m[ikeep[ikeep>=1],]
    elimna
}

pcorhc4_err <- function(x, y, CN = FALSE)
{
    z1 <- (x - mean(x)) / sd(x)
    z2 <- (y - mean(y)) / sd(y)
    ans <- olshc4(z1, z2, alpha = 0.05, CN = CN)
    return(ans$ci[2, 6])
}

pcorhc4<-function(x,y,alpha=.05,CN=F)
{
    #
    #   Compute a .95 confidence interval for Pearson's correlation coefficient.
    #   using the HC4 method
    #
    # CN=F, degrees of freedom are n-p; seems better for general use.
    # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
    #
    xy<-elimna(cbind(x,y))
    x<-xy[,1]
    y<-xy[,2]
    z1=(x-mean(x))/sqrt(var(x))
    z2=(y-mean(y))/sqrt(var(y))
    ans=olshc4(z1,z2,alpha=alpha,CN=CN)
    list(r=ans$r,ci=ans$ci[2,3:4],p.value=ans$ci[2,5])
}

olshc4<-function(x,y,alpha=.05,CN=FALSE,xout=FALSE,outfun=outpro,HC3=FALSE,...)
{
    #
    # Compute confidence for least squares
    # regression using heteroscedastic method
    # recommended by Cribari-Neto (2004).
    # CN=F, degrees of freedom are n-p
    # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
    # All indications are that CN=F is best for general use.
    #
    #  HC3=TRUE, will replace the HC4 estimator with the HC3 estimator.
    #
    x<-as.matrix(x)
    if(nrow(x) != length(y))stop("Length of y does not match number of x values")
    m<-cbind(x,y)
    m<-elimna(m)
    y<-m[,ncol(x)+1]
    x=m[,1:ncol(x)]
    n=length(y)
    nrem=n
    n.keep=length(y)
    x<-as.matrix(x)
    if(xout){
        flag<-outfun(x,...)$keep
        x<-as.matrix(x)
        x<-x[flag,]
        y<-y[flag]
        n.keep=length(y)
        x<-as.matrix(x)
    }
    temp<-lsfit(x,y)
    x<-cbind(rep(1,nrow(x)),x)
    xtx<-solve(t(x)%*%x)
    h<-diag(x%*%xtx%*%t(x))
    n<-length(h)
    d<-(n*h)/sum(h)
    for(i in 1:length(d)){
        d[i]<-min(4, d[i])
    }
    if(HC3)d=2
    hc4<-xtx%*%t(x)%*%diag(temp$res^2/(1-h)^d)%*%x%*%xtx
    df<-nrow(x)-ncol(x)
    crit<-qt(1-alpha/2,df)
    if(CN)crit=qnorm(1-alpha/2)
    al<-ncol(x)
    p=al-1
    ci<-matrix(NA,nrow=al,ncol=6)
    lab.out=rep("Slope",p)
    dimnames(ci)<-list(c("(Intercept)",lab.out),c("Coef.","Estimates",
                                                  "ci.lower","ci.upper","p-value","Std.Error"))
    for(j in 1:al){
        ci[j,1]<-j-1
        ci[j,2]<-temp$coef[j]
        ci[j,3]<-temp$coef[j]-crit*sqrt(hc4[j,j])
        ci[j,4]<-temp$coef[j]+crit*sqrt(hc4[j,j])
        test<-temp$coef[j]/sqrt(hc4[j,j])
        ci[j,5]<-2*(1-pt(abs(test),df))
        if(CN)ci[j,5]<-2*(1-pnorm(abs(test),df))
    }
    ci[,6]=sqrt(diag(hc4))
    list(n=nrem,n.keep=n.keep,ci=ci, cov=hc4)
}




