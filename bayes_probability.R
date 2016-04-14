# returns 3D array (i, j, k)
# probability of observing lag coordinate j of observed vector k 
# given actual state = lib vector i
compute_error_probs <- function(data, obs_data, pred, lib)
{
    obs_err_cv = 0.1
    obs_err_sd <- sqrt(log(1+obs_err_cv^2))
    obs_err_mean <- -(obs_err_sd^2/2)
    ni <- length(lib)
    nj <- 9
    nk <- length(pred)
    
    # setup matrices
    lib_vectors <- cbind(data[lib,1], data[lib-1,1], data[lib-2,1], 
                         data[lib,2], data[lib-1,2], data[lib-2,2], 
                         data[lib,3], data[lib-1,3], data[lib-2,3])
    pred_vectors <- cbind(obs_data[pred,1], obs_data[pred-1,1], obs_data[pred-2,1], 
                          obs_data[pred,2], obs_data[pred-1,2], obs_data[pred-2,2], 
                          obs_data[pred,3], obs_data[pred-1,3], obs_data[pred-2,3])
    error_probs <- array(NA, dim = c(ni, nj, nk))
    
    for(k in seq_along(pred))
    {
        errors <- log(lib_vectors) - matrix(log(pred_vectors[k,]), nrow = ni, ncol = nj, byrow = TRUE)
        error_probs[,,k] <- dnorm(errors, mean = obs_err_mean, sd = obs_err_sd)
    }
    return(error_probs)
}

make_prec_plot <- function(file, out_file, names = c("Univariate (x)", 
                                                     "Univariate (y)", 
                                                     "Univariate (z)", 
                                                     "Best Multivariate", 
                                                     "Multiview"))
{
    load(file)
    
    embeddings_matrix <- do.call(rbind, embeddings_list)
    univar_x_index <- match(data.frame(1:3), embeddings_list)
    univar_y_index <- match(data.frame(4:6), embeddings_list)
    univar_z_index <- match(data.frame(7:9), embeddings_list)
    
    prec_list <- lapply(prec_list, function(x) -log(x))
    best_embedding_index <- which.max(sapply(prec_list, median))
    multiview_prec <- -log(multiview_prec)
    
    prec_df <- data.frame(univar_x = prec_list[[univar_x_index]], 
                          univar_y = prec_list[[univar_y_index]], 
                          univar_z = prec_list[[univar_z_index]], 
                          best_multivar = prec_list[[best_embedding_index]], 
                          multiview = multiview_prec)
    
    pdf(file = out_file, width = 8, height = 6, useDingbats = FALSE)
    par(mar = c(4,4,1,1))
    boxplot(prec_df, ylab = "Precision (log-scale)", names = names)
    dev.off()
}

if(FALSE) # do calculations
{
    compute_prec <- function(coords)
    {
        sapply(seq_along(pred), function(layer) {
            p <- apply(error_probs[,coords,layer],1, prod)
            s <- data[pred[layer],]
            
            prec <- 1
            for(j in 1:3)
            {
                mse <- sum((data[lib,j] - s[j])^2 * p) / sum(p)
                null_mse <- mean((data[lib,j] - s[j])^2)
                prec <- prec * mse / null_mse
            }
            return(prec)
        })
    }
    
    file <- "model_1.Rdata"
    outfile <- "model_1_prec.Rdata"
    load(file)
    pred <- 501:2500
    lib <- 3001:11000
    error_probs <- compute_error_probs(data, obs_data, pred, lib)
    
    embeddings_list <- combn(1:9, 3, simplify = FALSE)
    prec_list <- lapply(embeddings_list, compute_prec)
    multiview_prec <- compute_prec(coords = 1:9)
    save(embeddings_list, prec_list, multiview_prec, file = outfile)
    
    file <- "model_2.Rdata"
    outfile <- "model_2_prec.Rdata"
    load(file)
    pred <- 501:2500
    lib <- 3001:11000
    error_probs <- compute_error_probs(data, obs_data, pred, lib)
    
    embeddings_list <- combn(1:9, 3, simplify = FALSE)
    prec_list <- lapply(embeddings_list, compute_prec)
    multiview_prec <- compute_prec(coords = 1:9)
    save(embeddings_list, prec_list, multiview_prec, file = outfile)
    
    file <- "model_3.Rdata"
    outfile <- "model_3_prec.Rdata"
    load(file)
    pred <- 501:2500
    lib <- 3001:11000
    error_probs <- compute_error_probs(data, obs_data, pred, lib)
    
    embeddings_list <- combn(1:9, 3, simplify = FALSE)
    prec_list <- lapply(embeddings_list, compute_prec)
    multiview_prec <- compute_prec(coords = 1:9)
    save(embeddings_list, prec_list, multiview_prec, file = outfile)
}

if(FALSE) # make plots
{
    make_prec_plot("model_1_prec.Rdata", "model_1_prec.pdf")
    make_prec_plot("model_2_prec.Rdata", "model_2_prec.pdf")
    make_prec_plot("model_3_prec.Rdata", "model_3_prec.pdf", 
                   names = c("Univariate\n(larvae)", 
                             "Univariate\n(pupae)", 
                             "Univariate\n(adults)", 
                             "Best Multivariate", 
                             "Multiview"))
}


# mean (weighted sum)
# mean_ws <- sum(data[lib,1] * p) / sum(p)
# variance (weighted sum)
# var_ws <- sum((data[lib,1])^2 * p) / sum(p) - mean_ws^2


# for each pred (layer)
# compute MSE, conditioned on 


#pdf(file = "state_probability_good.pdf", width = 8, height = 6)
#par(mfrow = c(3,1), mar = c(3, 5, 0, 1), oma = c(0, 0, 1, 0))
#plot(data[lib,1], p, xlab = "", ylab = "P(x)", pch = 3)
#plot(data[lib,2], p, xlab = "", ylab = "P(y)", pch = 3)
#plot(data[lib,3], p, xlab = "", ylab = "P(z)", pch = 3)
#dev.off()
# 
# df <- data.frame(x = data[lib,1], y = data[lib,2], z = data[lib,3], p = p)
# df <- df[order(df$x),]
# px <- loess(p ~ x, df, span = 0.1)
# plot(data[lib,1], p, xlab = "", ylab = "P(x)", pch = 3, cex = 0.5)
# lines(df$x, predict(px), col = "blue")
# 
# # total p
# const <- integrate(function(x) {predict(px, data.frame(x = x))}, lower = min(df$x), upper = max(df$x))
# 
# # mean (weighted sum)
# mean_ws <- sum(data[lib,1] * p) / sum(p)
# 
# # mean (integrated using loess)
# mean_func <- function(x) {x*predict(px, data.frame(x = x))}
# mean_int <- integrate(mean_func, lower = min(df$x), upper = max(df$x))
# mean_int <- mean_int$value / const$value
# 
# # variance (weighted sum)
# var_ws <- sum((data[lib,1]-mean_ws)^2 * p) / sum(p)
# 
# # variance (integrated)
# var_func <- function(x) {x * x * predict(px, data.frame(x = x))}
# var_int <- integrate(func, lower = min(df$x), upper = max(df$x))
# var_int <- var_int$value / const$value - mean_int^2