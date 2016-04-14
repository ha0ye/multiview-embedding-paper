library(scatterplot3d)

make_figure_0 <- function(file, output = NULL)
{
  load(file)
  if(!is.null(output))
    pdf(file = output, width = 6.5, height = 6.5)
  layout(matrix(c(1, 2, 3), nrow = 3))
  par(mar = c(2, 2, 1, 1), oma = c(1, 3, 0, 0))
  
  x <- 101:200
  for(j in 1:NCOL(data))
  {
    #plot(x, obs_data[x, j], type = "l")
    plot(x, data[x, j], type = "l", lty = 2)
  }
  mtext(c("x", "y", "z"), 2, line = 1, adj = 1-c(1/6, 3/6, 5/6), outer = TRUE)
  
  if(!is.null(output))
    dev.off()
  return()
}

make_figure_1 <- function(file)
{
  load(file)
  x <- normed_data[101:200,1]
  y <- normed_data[101:200,2]
  z <- normed_data[101:200,3]
 
  pdf("figure_1.pdf", width = 6.5, height = 2.4)
  par(mfrow = c(1, 3), mar = c(4, 3.5, 1, 1), oma = c(0, 1, 0, 0))
  plot.lim <- c(-3, 3)
  
  plot(x, y, cex.lab = 1.25, xlab = "", ylab = "", pch = 3, 
       xlim = plot.lim, ylim = plot.lim, asp = 1)
  mtext("x", side = 1, line = 2.5, cex = 5/6)
  mtext("y", side = 2, line = 2.5, cex = 5/6)
  
  plot(x, z, cex.lab = 1.25, xlab = "", ylab = "", pch = 3, 
       xlim = plot.lim, ylim = plot.lim, asp = 1)
  mtext("x", side = 1, line = 2.5, cex = 5/6)
  mtext("z", side = 2, line = 2.5, cex = 5/6)
  
  plot(y, z, cex.lab = 1.25, xlab = "", ylab = "", pch = 3, 
       xlim = plot.lim, ylim = plot.lim, asp = 1,)
  mtext("y", side = 1, line = 2.5, cex = 5/6)
  mtext("z", side = 2, line = 2.5, cex = 5/6)

  dev.off()
  
  return()
}

make_figure_2 <- function(file)
{
  splinify <- function(ts, spline_density = 20)
  {
    n <- length(ts)  
    times <- seq(1, n-2, len = (n-1)*spline_density)
  
    # make spline function
    xs <- splinefun(1:n, ts)
    
    return(xs(times))
  }
  
  load(file)
  
  x <- normed_data[101:200,1]
  y <- normed_data[101:200,2]
  z <- normed_data[101:200,3]
  
  n <- length(x)
  
  pdf("figure_2.pdf", width = 6.5, height = 4.5)
  par(mfrow = c(2, 3), mar = c(4, 3.5, 1, 1), oma = c(0, 1, 0, 0))
  plot.lim <- c(-3, 3)
  
  plot(splinify(x[2:n]), splinify(y[1:(n-1)]), type = "l", 
       xlim = plot.lim, ylim = plot.lim, asp = 1, 
       xlab = "", ylab = "", cex.lab = 1.25)
  mtext("x(t)", side = 1, line = 2.5, cex = 5/6)
  mtext("y(t+1)", side = 2, line = 2.5, cex = 5/6)
  
  plot(splinify(x[3:n]), splinify(y[1:(n-2)]), type = "l", 
       xlim = plot.lim, ylim = plot.lim, asp = 1, 
       xlab = "", ylab = "", cex.lab = 1.25)
  mtext("x(t)", side = 1, line = 2.5, cex = 5/6)
  mtext("y(t+2)", side = 2, line = 2.5, cex = 5/6)
  
  plot(splinify(x[4:n]), splinify(y[1:(n-3)]), type = "l", 
       xlim = plot.lim, ylim = plot.lim, asp = 1, 
       xlab = "", ylab = "", cex.lab = 1.25)
  mtext("x(t)", side = 1, line = 2.5, cex = 5/6)
  mtext("y(t+3)", side = 2, line = 2.5, cex = 5/6)
  
  plot(splinify(x[2:n]), splinify(z[1:(n-1)]), type = "l", 
       xlim = plot.lim, ylim = plot.lim, asp = 1, 
       xlab = "", ylab = "", cex.lab = 1.25)
  mtext("x(t)", side = 1, line = 2.5, cex = 5/6)
  mtext("z(t+1)", side = 2, line = 2.5, cex = 5/6)
  
  plot(splinify(x[3:n]), splinify(z[1:(n-2)]), type = "l", 
       xlim = plot.lim, ylim = plot.lim, asp = 1, 
       xlab = "", ylab = "", cex.lab = 1.25)
  mtext("x(t)", side = 1, line = 2.5, cex = 5/6)
  mtext("z(t+2)", side = 2, line = 2.5, cex = 5/6)
  
  plot(splinify(x[4:n]), splinify(z[1:(n-3)]), type = "l", 
       xlim = plot.lim, ylim = plot.lim, asp = 1, 
       xlab = "", ylab = "", cex.lab = 1.25)
  mtext("x(t)", side = 1, line = 2.5, cex = 5/6)
  mtext("z(t+3)", side = 2, line = 2.5, cex = 5/6)

  dev.off()
  
  return()
}

make_figure_3 <- function(file, lib = 25)
{
  load(file)
  
  index <- paste("lib=", lib, sep = "")
  results <- combined_results[[which(names(combined_results) == index)]]
  plot_pch <- 1:num_vars

  in_rhos <- do.call(c, lapply(results, function(x) {x$in_stats$rho}))
  out_rhos <- do.call(c, lapply(results, function(x) {x$out_stats$rho}))
  plot_symbols <- rep(plot_pch, each = NROW(embeddings_list))
  legend_text <- sapply(1:num_vars, function(x) {paste("predicting y_", x, sep = "")})
  
  pdf(file = "figure_3.pdf", width = 6.5, height = 6)
  par(mar = c(4.5, 4.5, 1, 1))
  plot(in_rhos, out_rhos, pch = plot_symbols, xlim = c(-0.8, 0.9), ylim = c(-0.2, 0.9), 
       xaxp = c(-0.8, 0.8, 4), 
       xlab = expression(rho ~ "(in-sample)"), 
       ylab = expression(rho ~ "(out-of-sample)"))
  legend(-0.8, 0.9, legend = legend_text, xjust = 0, yjust = 1, pch = plot_pch)
  abline(a = 0, b = 1, lty = 2)
  dev.off()
  
  return()
}

make_figure_4 <- function(file, lib = 25)
{
  load(file)
  
  index <- paste("lib=", lib, sep = "")
  results <- combined_results[[which(names(combined_results) == index)]]
  plot_pch <- 1:num_vars
  max_lag <- 3
  
  #other_embeddings <- embeddings_list[apply((embeddings_list-1) %% max_lag, 1, min) != 0,]
  #ordered_embeddings_list <- rbind(other_embeddings, embeddings_list)
  total_lag <- jitter(rep(rowSums((embeddings_list - 1) %% max_lag), times = num_vars))
  out_rhos <- do.call(c, lapply(results, function(x) {x$out_stats$rho}))
  plot_symbols <- rep(plot_pch, each = NROW(embeddings_list))
  legend_text <- sapply(1:num_vars, function(x) {paste("predicting y_", x, sep = "")})
  
  pdf(file = "figure_4.pdf", width = 6.5, height = 6)
  par(mar = c(4.5, 4.5, 1, 1))
  plot(total_lag, out_rhos, pch = plot_symbols, ylim = c(-0.2, 0.8), 
       xlab = "total time lag", 
       ylab = expression(rho ~ "(out-of-sample)"))
  legend(0, -0.2, legend = legend_text, xjust = 0, yjust = 0, pch = plot_pch)
  dev.off()
  
  return()
}

make_figure_5 <- function(file, title = file, output = NULL, line_width = 1.5)
{
  make_plots <- function(values_to_plot, x_axt, sub_label = "", line_width)
  {
    vals <- aggregate(values_to_plot, list(values_to_plot$lib), mean)
    vals <- vals[,3:5]
    y_limits <- range(vals)
    
    plot(vals[,1], col = "gray40", type = "l", ylim = y_limits, xaxt = "n", 
         xlab = "Library Size", ylab = "", lwd = line_width)
    lines(vals[,2], col = "red3", lwd = line_width)
    lines(vals[,3], col = "royalblue", lwd = line_width)
    mtext(sub_label, 3, at = 0.08, line = 0, adj = 0)
    if(x_axt)
      axis(1, at = seq_along(lib_sizes), lib_sizes)
    else
      axis(1, at = seq_along(lib_sizes), labels = FALSE)
    return()
  }
  
  load(file)
  lib_sizes <- sort(unique(combined_stats$lib))
  num_vars <- length(unique(combined_stats$target))
  if(!is.null(output))
    pdf(file = output, width = 6.5, height = 6.5)
  layout(matrix(c(1, 2, 3, 10, 4, 5, 6, 10, 7, 8, 9, 10), nrow = 4), heights = c(3, 3, 3, 1.3))
  par(mar = c(1, 2, 1, 1), oma = c(0.5, 2, 2, 0))
  
  letters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  
  for(j in 1:num_vars)
  {
    rows <- combined_stats$target == j
    make_plots(combined_stats[rows, c("lib_size", "univar_rho", "multivar_rho", "multiembed_rho")], x_axt = FALSE, letters[j], line_width)
    make_plots(combined_stats[rows, c("lib_size", "univar_mae", "multivar_mae", "multiembed_mae")], x_axt = FALSE, letters[j+3], line_width)
    make_plots(combined_stats[rows, c("lib_size", "univar_rmse", "multivar_rmse", "multiembed_rmse")], x_axt = TRUE, letters[j+6], line_width)
  }
  
  par(mar = c(0, 0, 2, 0))
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE)
  legend(0.5, 0, c("univariate", "best multivariate", "multiview embedding"), 
         col = c("gray40", "red3", "royalblue"), lwd = line_width, xjust = 0.5, yjust = 0, horiz = TRUE)
  #mtext(c("library size", "library size", "library size"), 3, at = c(1/6, 3/6, 5/6), adj = 0.5, line = -1, outer = FALSE)
  mtext("library size", 3, at = 0.5, adj = 0.5, line = -1, outer = FALSE)
  
  mtext(c("x", "y", "z"), 3, line = 0.5, at = c(1/6, 3/6, 5/6), outer = TRUE)
  # mtext(title, 3, line = 2, outer = TRUE)
  mtext(c(expression("accuracy (" ~ rho ~ ")"), "error (MAE)", "error (RMSE)"), 
        2, line = 0.5, at = 1-c(3/20, 9/20, 15/20), adj = 0.5, outer = TRUE)

  if(!is.null(output))
    dev.off()
  return()
}

make_figure_8 <- function(file, title = file, output = NULL, line_width = 1.5)
{
  make_plots <- function(values_to_plot, x_axt, sub_label = "")
  {
    vals <- aggregate(values_to_plot, list(values_to_plot$lib), mean)
    vals <- vals[,3:5]
    y_limits <- range(vals)
    
    plot(vals[,1], col = "gray40", type = "l", ylim = y_limits, xaxt = "n", 
         xlab = "", ylab = "", lwd = line_width)
    lines(vals[,2], col = "red3", lwd = line_width)
    lines(vals[,3], col = "royalblue", lwd = line_width)
    mtext(sub_label, 3, at = 0.08, line = -0.2, adj = 0)
    if(x_axt)
      axis(1, at = seq_along(lib_sizes), lib_sizes)
    else
      axis(1, at = seq_along(lib_sizes), labels = FALSE)
    return()
  }
  
  if(!is.null(output))
    pdf(file = output, width = 6.5, height = 6.5)
  layout(matrix(c(1, 2, 3, 10, 4, 5, 6, 10, 7, 8, 9, 10), nrow = 4), heights = c(3, 3, 3, 1.3))
  par(mar = c(1, 2, 1, 1), oma = c(0.5, 2, 2, 0))
  
  letters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  
  load("multiembed_all_results_1.Rdata")
  lib_sizes <- sort(unique(combined_stats$lib))
  num_vars <- length(unique(combined_stats$target))
  
  for(j in 1:num_vars)
  {
    rows <- combined_stats$target == j
    make_plots(combined_stats[rows, c("lib", "univar_rho", "multivar_rho", "multiembed_rho")], 
               x_axt = (j==num_vars), letters[j])
  }

  load("multiembed_all_results_2.Rdata")
  for(j in 1:num_vars)
  {
    rows <- combined_stats$target == j
    make_plots(combined_stats[rows, c("lib", "univar_rho", "multivar_rho", "multiembed_rho")], 
               x_axt = (j==num_vars), letters[j+3])
  }  
  
  load("multiembed_all_results_3.Rdata")
  for(j in 1:num_vars)
  {
    rows <- combined_stats$target == j
    make_plots(combined_stats[rows, c("lib", "univar_rho", "multivar_rho", "multiembed_rho")], 
               x_axt = (j==num_vars), letters[j+6])
  }
  
  par(mar = c(0, 0, 2, 0))
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE)
  legend(0.5, 0, c("univariate", "best multivariate", "multiview embedding"), 
         col = c("gray40", "red3", "royalblue"), lwd = line_width, xjust = 0.5, yjust = 0, horiz = TRUE)
  #mtext(c("library size", "library size", "library size"), 3, at = c(1/6, 3/6, 5/6), adj = 0.5, line = -1, outer = FALSE)
  mtext("library size", 3, at = 0.5, adj = 0.5, line = -1, outer = FALSE)
  
  mtext(c("coupled logistic", "food chain", "flour beetle"), 3, line = 0.5, at = c(1/6, 3/6, 5/6), outer = TRUE)
  # mtext(title, 3, line = 2, outer = TRUE)
  mtext(c(expression("accuracy (" ~ rho ~ ")"), expression("accuracy (" ~ rho ~ ")"), expression("accuracy (" ~ rho ~ ")")), 
        2, line = 0.5, at = 1-c(3/20, 9/20, 15/20), adj = 0.5, outer = TRUE)
  
  if(!is.null(output))
    dev.off()
  return()

}

make_figure_9 <- function(output = NULL)
{
  palette = c("red3", "royalblue", "green4")
  
  make_attractor_plot <- function(columns = c(1, 2, 3), time_delay = c(0, 1, 2), attractor_color)
  {
    var_labels <- c("x", "y", "z")
    labels <- ifelse(time_delay > 0, 
                     paste(var_labels[columns], "(t-", time_delay, ")", sep = ""), 
                     paste(var_labels[columns], "(t)", sep = ""))
    labels = c("", "", "")
    to_plot <- cbind(data[t-time_delay[1]*800, columns[1]], 
                     data[t-time_delay[2]*800, columns[2]], 
                     data[t-time_delay[3]*800, columns[3]])
    labels[2] <- ""
    scatterplot3d(to_plot, type = "l", mar = c(1.5, 1, 0, 0), angle = 60, 
                  tick.marks = FALSE, lwd = 1.5, color = attractor_color, 
                  xlab = labels[1], ylab = labels[2], zlab = labels[3])
    return()
  }
  
  load("model_2_alt.Rdata")
  palette = c("red3", "royalblue", "green4")
  
  if(!is.null(output))
    pdf(file = output, width = 6.5, height = 5)
  layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), nrow = 3, byrow = TRUE), heights = c(3, 5, 5))
  par(oma = c(0.5, 2, 2, 0))
  
  t <- seq(from = 121*800, to = 180*800, by = 80)
  y_labels <- c("x", "y", "z")
  x <- data[t,1]
  y <- data[t,2]
  z <- data[t,3]
  z <- (z - min(z)) / (max(z) - min(z))
  ### PANEL A
  par(mar = c(1, 2, 0, 1), xaxt = "n", yaxt = "n")
  plot(t, x, xlim = c(126*800, 295*800), type = "l", lty = 1, lwd = 1.5, col = palette[1])
  lines(t+60*800, y, type = "l", lty = 1, lwd = 1.5, col = palette[2])
  lines(t+120*800, z, type = "l", lty = 1, lwd = 1.5, col = palette[3])
  mtext("A", side = 3, line = 0, adj = -0.05, font = 2, cex = 14/12)
  
  abline(v = 180*800+400, lty = 2, col = "gray40")
  abline(v = 240*800+400, lty = 2, col = "gray40")
  
  ### PANEL B
  #
  make_attractor_plot(columns = c(1, 1, 1), attractor_color = palette[1])
  mtext("B", side = 3, line = 0, adj = -0.05, font = 2, cex = 14/12)
  make_attractor_plot(columns = c(2, 2, 2), attractor_color = palette[2])
  make_attractor_plot(columns = c(3, 3, 3), attractor_color = palette[3])
  
  make_attractor_plot(columns = c(1, 2, 2), attractor_color = "darkorchid2")
  mtext("C", side = 3, line = 0, adj = -0.05, font = 2, cex = 14/12)
  make_attractor_plot(columns = c(2, 3, 3), attractor_color = "darkorchid2")
  make_attractor_plot(columns = c(3, 1, 1), attractor_color = "darkorchid2")
  
  if(!is.null(output))
    dev.off()
  return()
}

#make_figure_9("figure_9.pdf")

# multiembed rho comparison between number of embedddings merged
make_figure_6 <- function(file, lib = 25)
{
  load(file)
  
  combined_rhos <- combined_output[combined_output$lib == lib, c(5:NCOL(combined_output))]
  combined_rhos <- sapply(combined_rhos, function(x) {pmax(0, x)})

  pdf(file = "figure_6.pdf", width = 6.5, height = 7.5)
  par(mfrow = c(num_vars, 1), mar = c(4.5, 4.5, 1, 1))
  
  for(j in 1:num_vars)
  {
    temp_str <- expression(rho ~ "(predicting y_" * j)
    plot(1:NCOL(combined_rhos), combined_rhos[j, ], type = "l", col = "royalblue", lwd = 2, 
            xlab = expression(kappa ~ "(number of embeddings)"), 
         ylab = bquote(rho ~ "(predicting y_" * .(j) * ")"))
  }
  
  dev.off()
  return()
}

make_figure_7 <- function(output = NULL)
{
  make_attractor_plot <- function(columns = c(1, 2, 3), time_delay = c(0, 0, 0))
  {
    var_labels <- c("x", "y", "z")
    labels <- ifelse(time_delay > 0, 
                     paste(var_labels[columns], "(t-", time_delay, ")", sep = ""), 
                     paste(var_labels[columns], "(t)", sep = ""))
    to_plot <- cbind(data[x-time_delay[1]*800, columns[1]], 
                     data[x-time_delay[2]*800, columns[2]], 
                     data[x-time_delay[3]*800, columns[3]])
    labels[2] <- ""
    scatterplot3d(to_plot, type = "l", mar = c(1.5, 1, 0, 0), angle = 60, 
                  tick.marks = FALSE, lwd = 1.5, 
                  xlab = labels[1], ylab = labels[2], zlab = labels[3])
    return()
  }
  
  load("model_2_alt.Rdata")
  
  if(!is.null(output))
    pdf(file = output, width = 6.5, height = 8.0)
  #layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7), nrow = 5, byrow = TRUE), heights = c(1, 1, 1, 3, 1))
  #layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(2, 1))
  layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 4, 4, 6), nrow = 5, byrow = TRUE), heights = c(1, 1, 1, 1.8, 1.8))
  par(oma = c(0.5, 2, 2, 0))
  
  x <- seq(from = 121*800, to = 180*800, by = 80)
  y_labels <- c("x", "y", "z")
  ### PANEL A
  par(mar = c(1, 2, 0, 1), xaxt = "n")
  for(j in 1:NCOL(data))
  {
    plot(x, data[x, j], type = "l", lty = 1, lwd = 1.5)
    mtext(y_labels[j], side = 2, line = 2.5, adj = 0.5)
    if(j == 1)
      mtext("A", side = 3, line = 0, adj = -0.05)
  }

  ### PANEL B
  make_attractor_plot(columns = c(2, 1, 3))
  mtext("B", side = 3, line = -0.5, adj = -0.04)
  
  ### PANEL C
  make_attractor_plot(columns = c(3, 3, 2), time_delay = c(0, 1, 2))
  mtext("C", side = 3, line = -0.5, adj = 0)
  
  ### PANEL D
  make_attractor_plot(columns = c(2, 1, 1), time_delay = c(0, 1, 0))
  mtext("D", side = 3, line = -0.5, adj = 0)
  
  ### PANEL D
  # make_attractor_plot(columns = c(3, 3, 2), time_delay = c(0, 1, 2))
  
  
#  par(mar = c(0, 0, 2, 0))
#  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE)
  
  if(!is.null(output))
    dev.off()
  return()
}

make_figure_7_new <- function(output = NULL)
{
  make_attractor_plot <- function(columns = c(1, 2, 3), time_delay = c(0, 0, 0))
  {
    var_labels <- c("x", "y", "z")
    labels <- ifelse(time_delay > 0, 
                     paste(var_labels[columns], "(t-", time_delay, ")", sep = ""), 
                     paste(var_labels[columns], "(t)", sep = ""))
    to_plot <- cbind(data[x-time_delay[1]*800, columns[1]], 
                     data[x-time_delay[2]*800, columns[2]], 
                     data[x-time_delay[3]*800, columns[3]])
    labels[2] <- ""
    scatterplot3d(to_plot, type = "l", mar = c(1.5, 1, 0, 0), angle = 60, 
                  tick.marks = FALSE, lwd = 1.5, color = "gray20", 
                  xlab = labels[1], ylab = labels[2], zlab = labels[3])
    return()
  }
  
  load("model_2_alt.Rdata")
  palette = c("red3", "royalblue", "green4")
  
  if(!is.null(output))
    pdf(file = output, width = 6.5, height = 3)
  #layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7), nrow = 5, byrow = TRUE), heights = c(1, 1, 1, 3, 1))
  #layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(2, 1))
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 3, byrow = FALSE))
  par(oma = c(0.5, 2, 2, 0))
  
  x <- seq(from = 121*800, to = 180*800, by = 80)
  y_labels <- c("x", "y", "z")
  ### PANEL A
  par(mar = c(1, 2, 0, 1), xaxt = "n")
  for(j in 1:NCOL(data))
  {
    plot(x, data[x, j], type = "l", lty = 1, lwd = 1.5, col = palette[j])
    mtext(y_labels[j], side = 2, line = 2.5, adj = 0.5)
    if(j == 1)
      mtext("A", side = 3, line = 0, adj = -0.05)
  }
  
  ### PANEL B
  make_attractor_plot(columns = c(2, 1, 3))
  mtext("B", side = 3, line = 0, adj = -0.05)
  
  if(!is.null(output))
    dev.off()
  return()
}

make_figures_univar <- function()
{
  load("univariate_results.Rdata")
  
  make_figure_univar(univariate_results_x, "x.EvsRho.pdf")
  make_figure_univar(univariate_results_y, "y.EvsRho.pdf")
  make_figure_univar(univariate_results_z, "z.EvsRho.pdf")
  
  return()
}

make_figure_univar <- function(univariate_results, file = "univariate.EvsRho.pdf")
{
  # colors
  curr_color <- 1
  palette(c("red", "orange", "yellow", "green", "darkgreen", "cyan", "blue", "purple", "black"))
  ylim <- c(max(0, min(univariate_results$rho)), max(univariate_results$rho))
  temp_data <- univariate_results[univariate_results$lib == lib_sizes[1],]
  
  pdf(file = file, width = 6, height = 6)
  
  plot(temp_data$E, pmax(0, temp_data$rho), type = "l", ylim = ylim, 
       xlab = "E", ylab = expression(rho), col = curr_color, lwd = 2)
  curr_color <- curr_color + 1
  for (l in lib_sizes[2:length(lib_sizes)])
  {
    temp_data <- univariate_results[univariate_results$lib == l,]
    lines(temp_data$E, pmax(0, temp_data$rho), col = curr_color, lwd = 2)
    curr_color <- curr_color + 1
  }
  legend(1, ylim[1], legend = lib_sizes, col = palette(), lwd = 2, xjust = 0, yjust = 0, title = "Library Size")
  dev.off()
  
  return()
}
