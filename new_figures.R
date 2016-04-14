rm(list = ls())

if(FALSE)
{
    make_figure_1a("fig_1a.pdf")
    make_figure_1b("fig_1b.pdf")
    make_figure_1c("fig_1c.pdf")
    make_figure_1d("fig_1d.pdf")
    
    plot_lorenz_c("figure_1c.pdf")
}

make_figure_1a <- function(file = NULL)
{
    library(rgl)
    load("model_2_alt.Rdata")
    k <- 155*800
    t <- seq(from = 121*800, to = k, by = 80)
    t2 <- seq(from = k, to = 175.5*800, by = 80)
    plot3d(data[t,], type = "l", xlab = "", ylab = "", zlab = "", 
           lwd = 2, xlim = c(0, 1), ylim = c(0, 1.1), zlim = c(2, 6.2), axes = FALSE)
    axes3d(edges = c("x-+", "y++", "z--"))
    #box3d()
    
    points3d(data[k,1], data[k,2], data[k,3], size = 4)
    lines3d(c(data[k,1], 1, 1), 
            c(data[k,2], data[k,2], data[k,2]), 
            c(data[k,3], data[k,3], 6.2), col = "red")
    lines3d(data[t2,], col = "gray")
    
    p <- matrix(c(0.395, -0.041, 0.918, 0, 
                  0.210, 0.977, -0.047, 0, 
                  -0.894, 0.211, 0.395, 0, 
                  0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(userMatrix = p, windowRect = c(0, 70, 800, 870))
    
    if(!is.null(file))
    {
        rgl.postscript(file, "pdf")
    }
}

make_figure_1b <- function(file = NULL)
{
    load("model_2_alt.Rdata")
    k <- 155*800
    t <- seq(from = 121*800, to = k, by = 80)
    t2 <- seq(from = k, to = 175.5*800, by = 80)
    if(!is.null(file))
    {
        pdf(file = file, width = 6, height = 4)
    }
    plot(t, data[t,2], type = "l", col = "red", xlim = rev(range(t)), lwd = 2)
    if(!is.null(file))
    {
        dev.off()
    }
}

make_figure_1c <- function(file = NULL)
{
    load("model_2_alt.Rdata")
    k <- 155*800
    t <- seq(from = 121*800, to = k, by = 80)
    t2 <- seq(from = k, to = 175.5*800, by = 80)
    if(!is.null(file))
    {
        pdf(file = file, width = 6, height = 4)
    }
    par(mfrow = c(3, 1), mar = c(2, 2, 0, 0), oma = c(1, 1, 1, 1))
    plot(t, data[t,2], type = "l", col = "#FF0000", lwd = 2, xlim = range(c(t, t2)))
    lines(t2, data[t2,2], type = "l", col = "#FF0000", lwd = 2, lty = 2, xlim = range(c(t, t2)))

    plot(t, data[t-800,2], type = "l", col = "#BB0000", lwd = 2, xlim = range(c(t, t2)))
    lines(t2, data[t2-800,2], type = "l", col = "#BB0000", lwd = 2, lty = 2, xlim = range(c(t, t2)))
    
    plot(t, data[t-1600,2], type = "l", col = "#770000", lwd = 2, xlim = range(c(t, t2)))
    lines(t2, data[t2-1600,2], type = "l", col = "#770000", lwd = 2, lty = 2, xlim = range(c(t, t2)))
    if(!is.null(file))
    {
        dev.off()
    }
}

make_figure_1d <- function(file = NULL)
{
    library(rgl)
    load("model_2_alt.Rdata")
    k <- 155*800
    t <- seq(from = 121*800, to = k, by = 80)
    t2 <- seq(from = k, to = 175.5*800, by = 80)
    to_plot <- cbind(data[t, 2], 
                     data[t-800, 2], 
                     data[t-2*800, 2])
    
    plot3d(to_plot, type = "l", xlab = "x", ylab = "y", zlab = "z", 
           lwd = 2, xlim = c(0, 1.1), ylim = c(0, 1.1), zlim = c(0, 1.1), axes = FALSE)
    axes3d(edges = c("x+-", "y-+", "z-+"))
    
    points3d(data[k,2], data[k-800,2], data[k-2*800,2], size = 4)
    to_plot <- cbind(data[t2, 2], 
                 data[t2-800, 2], 
                 data[t2-2*800, 2])
    lines3d(to_plot, col = "gray")
    
    p <- matrix(c(0.035, -0.741, 0.672, 0, 
                  0.985, -0.087, -0.148, 0, 
                  0.168, 0.667, 0.726, 0, 
                  0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(userMatrix = p, windowRect = c(0, 70, 800, 870))
    
    if(!is.null(file))
    {
        rgl.postscript(file, "pdf")
    }
}

make_lorenz_data <- function()
{
    # initialize params
    num_points <- 2000
    dt <- 0.01
    half_dt <- dt / 2
    x <- matrix(20, nrow = num_points, ncol = 3)
    
    f <- function(x) {
        c(10 * (x[2] - x[1]), # sigma = 10
          28 * x[1] - x[1] * x[3] - x[2], # rho = 28
          x[1] * x[2] - 8/3 * x[3]) # beta = 8/3
    }
    
    # generate data
    for(i in 2:num_points)
    {
        xx <- x[i-1,]
        k1 <- f(xx)
        x[i,] <- xx + k1 * dt
        #        k2 <- f(xx + half_dt * k1)
        #        k3 <- f(xx + half_dt * k2)
        #        k4 <- f(xx + dt * k3)
        #        x[i,] <- xx + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)
    }
    save(x, file = "lorenz_data.RData")
}

plot_lorenz_a <- function(file = NULL)
{
    library(rgl)
    load("lorenz_data.Rdata")
    k <- 1120
    t <- 1:k
    t2 <- k:(k+40)
    xlim <- c(-16.9, 20.6)
    ylim <- c(-21.3, 22.1)
    zlim <- c(10.8, 46.5)
    
    plot3d(x[t,], type = "l", xlab = "", ylab = "", zlab = "", 
           lwd = 2, xlim = xlim, ylim = ylim, zlim = zlim, axes = FALSE)
    axes3d(edges = c('x+-', 'y--', 'z--'))
    box3d()
    points3d(x[k,1], x[k,2], x[k,3], size = 4)
    lines3d(c(x[k,1], x[k,1], x[k,1]), 
            c(x[k,2], 21, 21), 
            c(x[k,3], x[k,3], 11), col = "red")
    lines3d(x[t2,], col = "gray")
    
    p <- matrix(c(-0.033, 0.491, -0.870, 0, 
                  0.979, 0.195, 0.071, 0, 
                  0.204, -0.848, -0.488, 0, 
                  0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(userMatrix = p, windowRect = c(0, 70, 800, 870))
    if(!is.null(file))
    {
        rgl.postscript(file, "pdf")
    }
    return()
}

plot_lorenz_b <- function(file = NULL, width = 6, height = 4)
{
    load("lorenz_data.Rdata")
    tt <- 1120
    t <- 1:tt
    t2 <- tt:(tt+40)
    if(!is.null(file))
        pdf(file = file, width = width, height = height)
    par(mar = c(4, 4, 1, 1))
    plot(t, x[t,1], type = "l", lwd = 2, xlab = "time", ylab = "X")
    lines(t2, x[t2,1], type = "l", lty = 2, col = "gray")
    if(!is.null(file))
        dev.off()
    return()
}

plot_lorenz_c <- function(file = NULL)
{
    library(rgl)
    load("lorenz_data.Rdata")
    tt <- 1088
    t <- 17:tt
    xlim <- c(-12, 4)
    ylim <- c(-18, -2)
    zlim <- c(-18, 2)
    
    # find segments
    tau <- 8
    d <- cbind(x[t,1], x[t-tau,1], x[t-2*tau,1])
    in_box <- d[,1] > xlim[1] & d[,1] < xlim[2] & 
        d[,2] > ylim[1] & d[,2] < ylim[2] & 
        d[,3] > zlim[1] & d[,3] < zlim[2]
    diff <- in_box[2:length(in_box)] - in_box[1:(length(in_box)-1)]
    starts <- which(diff == 1)+1
    ends <- which(diff == -1)
    if(in_box[1])
        starts <- c(1, starts)
    if(in_box[length(in_box)])
        ends <- c(ends, length(in_box))
    segments <- cbind(starts, ends)
    
    # find nearest neighbors
    num_neighbors <- 4
    del_t <- 5
    index <- rev(seq(from = tt-del_t, to = 17, by = -del_t))
    lib <- cbind(x[index,1], 
                 x[index-tau,1], 
                 x[index-2*tau,1])
    s <- matrix(c(x[tt,1], x[tt-tau,1], x[tt-2*tau,1]), nrow = NROW(lib), ncol = 3, byrow = TRUE)
    dist <- sqrt(rowSums((lib - s)^2))
    nn <- which(rank(dist) <= num_neighbors)
    weights <- exp(-dist[nn]/dist[nn[1]])
    weights <- weights / sum(weights)
    nn <- index[nn]
    
    plot3d(d[segments[1,1]:segments[1,2],], type = "l", xlab = "", ylab = "", zlab = "", 
           xlim = xlim, ylim = ylim, zlim = zlim, 
           lwd = 1, axes = FALSE)
    
    # plot segments
    for(j in 2:NROW(segments))
    {
        lines3d(d[segments[j,1]:segments[j,2],])
    }
    
    axes3d(edges = c('x--', 'y+-', 'z--'))
    box3d()
    pred <- array(0, dim = c(9, 3))
    for(k in seq_along(nn))
    {
        i <- nn[k]
        points3d(x[i,1], x[i-tau,1], x[i-2*tau,1], size = 4, pch = 5, col = "blue")
        ii <- i:(i+tau)
        lines3d(x[ii,1], x[ii-tau,1], x[ii-2*tau,1], col = "cyan", lwd = 3)
        deltas <- cbind(x[ii,1], x[ii-tau,1], x[ii-2*tau,1]) - 
            matrix(c(x[i,1], x[i-tau,1], x[i-2*tau,1]), nrow = length(ii), ncol = 3, byrow = TRUE)
        pred <- pred + deltas*weights[k]
    }
    points3d(x[tt,1], x[tt-tau,1], x[tt-2*tau,1], size = 4)
    pred <- pred + matrix(c(x[tt,1], x[tt-tau,1], x[tt-2*tau,1]), nrow = NROW(pred), ncol = 3, byrow = TRUE)
    lines3d(pred, col = "gray", lwd = 3)
    
    p <- matrix(c(0.884, 0.467, -0.008, 0, 
                 -0.182, 0.362, 0.914, 0, 
                 0.430, -0.807, 0.405, 0, 
                 0, 0, 0, 1), nrow = 4, byrow = TRUE)
    
    #    lines3d(c(x[k,1], x[k,1], x[k,1]), 
    #            c(x[k,2], 21, 21), 
    #            c(x[k,3], x[k,3], 11), col = "red")
    #    lines3d(x[t2,1], x[t2-8,1], x[t2-16,1], col = "gray")
    
    #    p <- matrix(c(-0.033, 0.491, -0.870, 0, 
    #                  0.979, 0.195, 0.071, 0, 
    #                  0.204, -0.848, -0.488, 0, 
    #                  0, 0, 0, 1), nrow = 4, byrow = TRUE)
    #p <- matrix(c(0.949, 0.319, 0.0233, 0, 
    #              -0.073, 0.142, 0.987, 0, 
    #              0.310, -0.936, 0.158, 0, 
    #              0, 0, 0, 1), nrow = 4, byrow = TRUE)
    #par3d(userMatrix = p)
    par3d(windowRect = c(0, 70, 800, 870))
    if(!is.null(file))
    {
        rgl.postscript(file, "pdf")
    }
    return()
}

make_lorenz_data_noisy <- function()
{
    # initialize params
    num_points <- 20000
    dt <- 0.01
    x <- matrix(20, nrow = num_points, ncol = 3)
    
    f <- function(x) {
        c(10 * (x[2] - x[1]), # sigma = 10
          28 * x[1] - x[1] * x[3] - x[2], # rho = 28
          x[1] * x[2] - 8/3 * x[3]) # beta = 8/3
    }
    
    # generate data
    for(i in 2:num_points)
    {
        xx <- x[i-1,]
        k1 <- f(xx)
        xx <- xx + k1 * dt
        x[i,] <- xx
    }
    
    data <- x
    obs_data <- x
    # add observational error
    for(j in 1:NCOL(x))
    {
        obs_data[,j] <- data[,j] + rnorm(NROW(x), mean = 0, sd = sd(data[,j])*0.2)
    }

    save(data, obs_data, file = "lorenz_data_noisy.RData")
}

# returns 3D array (i, j, k)
# probability of observing lag coordinate j of observed vector k 
# given actual state = lib vector i
compute_lorenz_error_probs <- function(data, obs_data, pred, lib, tau = 1)
{
    obs_err <- c(NA, NA, NA)
    for(j in 1:NCOL(data))
    {
        obs_err[j] <- sd(data[,j])*0.2
    }
    obs_err <- rep(obs_err, each = 3)
    
    ni <- length(lib)
    nj <- 9
    nk <- length(pred)
    
    # setup matrices
    lib_vectors <- cbind(data[lib,1], data[lib-tau,1], data[lib-2*tau,1], 
                         data[lib,2], data[lib-tau,2], data[lib-2*tau,2], 
                         data[lib,3], data[lib-tau,3], data[lib-2*tau,3])
    pred_vectors <- cbind(obs_data[pred,1], obs_data[pred-tau,1], obs_data[pred-2*tau,1], 
                          obs_data[pred,2], obs_data[pred-tau,2], obs_data[pred-2*tau,2], 
                          obs_data[pred,3], obs_data[pred-tau,3], obs_data[pred-2*tau,3])
    error_probs <- array(NA, dim = c(ni, nj, nk))
    
    for(k in seq_along(pred))
    {
        errors <- lib_vectors - matrix(pred_vectors[k,], nrow = ni, ncol = nj, byrow = TRUE)
        for(j in 1:NCOL(errors))
            error_probs[,j,k] <- dnorm(errors[,j], sd = obs_err[j])
#        error_probs[,,k] <- dnorm(errors, sd = obs_err[k])
    }
    return(error_probs)
}

# draw univariate reconstructions
make_univar_figures <- function(coord = 1, block_coords = 1:3)
{
    file_1 <- paste("nearest_neighbors_SSR_coord=", coord, ".pdf", sep = "")
    file_2 <- paste("nearest_neighbors_orig_coord=", coord, ".pdf", sep = "")
    
    library(rgl)
    load("lorenz_data_noisy.RData")
    
    # make vectors
    tau <- 4
    pred <- 1088
    lib <- 2001:20000
    error_probs <- compute_lorenz_error_probs(obs_data, obs_data, pred, lib, tau = 4)
    lib_vectors <- cbind(obs_data[lib,coord], obs_data[lib-tau,coord], obs_data[lib-2*tau,coord])
    pred_vectors <- cbind(obs_data[pred,coord], obs_data[pred-tau,coord], obs_data[pred-2*tau,coord])
    
    # find most probable nearest neighbors
    prob <- apply(error_probs[, block_coords, 1], 1, prod)
    prob <- prob / max(prob)
    wt_sum <- sum(prob)
    ord <- order(-prob)
    wt <- 0
    for(i in ord)
    {
        wt <- wt + prob[i]
        if(wt > wt_sum*0.95)
        {
            prob_threshold <- prob[i]
            break
        }
    }
    index <- prob >= prob_threshold
    
    # make plot (reconstructed state space)
    t1 <- (2*tau+1):4000
    t2 <- t1 - tau
    t3 <- t1 - 2*tau
    plot3d(data[t1, coord], data[t2, coord], data[t3, coord], type = "l", col = "#BFBFBF", 
           xlab = "", ylab = "", zlab = "", axes = FALSE) # actual lagged reconstruction
    points3d(pred_vectors) # observed point
    intensity <- prob[index]
    cols <- rgb(1, 0, intensity)
    points3d(lib_vectors[index,], col = cols) # nearest neighbors
    box3d()
    p <- matrix(c(1, 0, 0, 0, 
                  0, 0, 1, 0, 
                  0, -1, 0, 0, 
                  0, 0, 0, 1), nrow = 4, byrow = TRUE)
    par3d(windowRect = c(0, 70, 800, 870), userMatrix = p)
    rgl.postscript(file_1, "pdf")
    
    # make plot (original state space)
    plot3d(data[1:4000,], type = "l", col = "#BFBFBF", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(t(data[pred,]))
    points3d(data[lib[index],], col = cols)
    box3d()
    par3d(windowRect = c(0, 70, 800, 870), userMatrix = p)
    rgl.postscript(file_2, "pdf")
    
    return()
}
make_univar_figures()
make_univar_figures(coord = 2, block_coords = 4:6)
make_univar_figures(coord = 3, block_coords = 1:6)

# do bayes calculations on lorenz
bayes_calc_lorenz <- function()
{
    library(rgl)
    load("lorenz_data_noisy.RData")
    pred <- 1088
    lib <- 2001:20000
    
    error_probs <- compute_lorenz_error_probs(data, obs_data, pred, lib, tau = 4)
    k <- 1 # index of observed vector
    
    # compute posterior for univariate x reconstruction
    coords <- c(1, 2, 3)
    prob <- apply(error_probs[, coords, k], 1, prod)
    prob <- prob / max(prob)
    wt_sum <- sum(prob)
    ord <- order(-prob)
    wt <- 0
    for(i in ord)
    {
        wt <- wt + prob[i]
        if(wt > wt_sum*0.99)
        {
            prob_threshold <- prob[i]
            break
        }
    }
    index <- prob >= prob_threshold
    intensity <- prob[index]
    cols <- rgb(1, 0, intensity)
    
    # make plot
    plot3d(data[1:4000,], type = "l", col = "#BFBFBF", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(data[lib[index],], col = cols)
    box3d()
    par3d(windowRect = c(0, 70, 800, 870))
    rgl.postscript("x_prob_dens.pdf", "pdf")
    
    # compute posterior for univariate y reconstruction
    coords <- c(4, 5, 6)
    prob <- apply(error_probs[, coords, k], 1, prod)
    prob <- prob / max(prob)
    wt_sum <- sum(prob)
    ord <- order(-prob)
    wt <- 0
    for(i in ord)
    {
        wt <- wt + prob[i]
        if(wt > wt_sum*0.99)
        {
            prob_threshold <- prob[i]
            break
        }
    }
    index <- prob >= prob_threshold
    intensity <- prob[index]
    cols <- rgb(0, 1, intensity)
    
    # make plot
    plot3d(data[1:4000,], type = "l", col = "#BFBFBF", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(data[lib[index],], col = cols)
    box3d()
    par3d(windowRect = c(0, 70, 800, 870))
    rgl.postscript("y_prob_dens.pdf", "pdf")
    
    # compute posterior for multiview with x and y reconstruction
    coords <- c(1, 2, 3, 4, 5, 6)
    prob <- apply(error_probs[, coords, k], 1, prod)
    prob <- prob / max(prob)
    wt_sum <- sum(prob)
    ord <- order(-prob)
    wt <- 0
    for(i in ord)
    {
        wt <- wt + prob[i]
        if(wt > wt_sum*0.99)
        {
            prob_threshold <- prob[i]
            break
        }
    }
    index <- prob >= prob_threshold
    intensity <- prob[index]
    cols <- rgb(0, 1, intensity)
    
    # make plot
    plot3d(data[1:4000,], type = "l", col = "#BFBFBF", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    points3d(data[lib[index],], col = cols)
    box3d()
    par3d(windowRect = c(0, 70, 800, 870))
    rgl.postscript("xy_prob_dens.pdf", "pdf")
    
#     
#     plot(data[lib,1], data[lib,2], col = x_cols)
#     points(data[lib,1], data[lib,2], col = y_cols)
#     
#     coords <- c(1, 2, 3, 4, 5, 6)
#     prob <- apply(error_probs[, coords, k], 1, prod)
#     cols <- rgb(0, 0, 1, alpha = prob / max(prob))
#     points(data[lib,1], data[lib,2], col = cols)
#     
    
    
    return()
}

