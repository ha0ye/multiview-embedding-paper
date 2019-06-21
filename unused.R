
# plot_model_attractor <- function(data_file, t = 1000:3000, x = c(1,0), y = c(2,0), z = c(3,0), 
#                                  x_label = "x", y_label = "y", z_label = "z", 
#                                  plot_file = NULL)
# {
#     load(data_file)
#     make_vec <- function(vec_definition)
#     {
#         return(data[t-vec_definition[2], vec_definition[1]])
#     }
#     t <- t + max(x[2], y[2], z[2])
#     plot3d(make_vec(x), make_vec(y), make_vec(z), type = "l", 
#            axes = FALSE, xlab = x_label, ylab = y_label, zlab = z_label)
#     axes3d(edges = c("x--", "y+-", "z--"), tick = FALSE, labels = FALSE)
#     m <- matrix(c(0.881, 0.473, 0.020, 0, 
#                   -0.062, 0.071, 0.996, 0, 
#                   0.470, -0.878, 0.092, 0, 
#                   0, 0, 0, 1), nrow = 4, byrow = TRUE)
#     par3d(userMatrix = m, windowRect = c(100, 100, 900, 800))
#     if(!is.null(plot_file))
#     {
#         rgl.postscript(plot_file, fmt = "pdf")
#     }
#     return()
# }
# 
# make_HW_plot_panel <- function(data_file, t = 500:1000, x = c(1,0), y = c(2,0), z = c(4,0),
#                                col = rgb(1.0, 0.0, 0.0),
#                                x_label = "x", y_label = "y", z_label = "z",
#                                plot_file = NULL, m = NULL, edges = c("x-+", "y--", "z--"))
# {
#     make_vec <- function(vec_definition)
#     {
#         return(data[t-vec_definition[2], vec_definition[1]])
#     }
#     t <- t + max(x[2], y[2], z[2])
# 
#     if(is.null(m))
#     {
#         m <- matrix(c(0.015, 0.598, -0.801, 0,
#                       0.930, 0.287, 0.231, 0,
#                       0.368, -0.748, -0.552, 0,
#                       0, 0, 0, 1), byrow = TRUE, nrow = 4)
#     }
# 
#     load(data_file)
#     xx <- make_vec(x)
#     yy <- make_vec(y)
#     zz <- make_vec(z)
# 
#     plot3d(NA, NA, NA,
#            xlim = range(xx), ylim = range(yy), zlim = range(zz),
#            xlab = x_label, ylab = y_label, zlab = z_label, axes = FALSE)
#     axes3d(edges = edges)
#     box3d()
#     par3d(userMatrix = m, zoom = 1.25, windowRect = c(50, 50, 850, 850))
#     lines3d(xx, yy, zz, col = col, lwd = 2)
#     if(!is.null(plot_file))
#     {
#         rgl.postscript(plot_file, fmt = "pdf")
#     }
#     return()
# }

# 
# # attractor plots of food chain model
# if(FALSE)
# {
#     plot_model_attractor("model_2_alt.Rdata")
#     plot_model_attractor("model_2_alt.Rdata", y = c(1,6), z = c(1,12))
#     plot_model_attractor("model_2_alt.Rdata", x = c(2,0), y = c(2,6), z = c(2,12))
#     plot_model_attractor("model_2_alt.Rdata", x = c(3,0), y = c(3,6), z = c(3,12))
#     plot_model_attractor("model_2_alt.Rdata", x = c(1,0), y = c(2,0), z = c(1,6))
# }
# 
# # attractor plots for Huisman-Weissing model
# if(FALSE)
# {
#     load("model_hw.Rdata")
#     t <- 500:1000
#     palette <- c(rgb(1.0, 0.0, 0.0), rgb(0.0, 0.0, 1.0))
#     point_color <- palette[1+(data[t,4] > data[t,2])]
#     
#     make_HW_plot_panel("model_hw.Rdata", t, col = point_color, 
#                     plot_file = "hw_native.pdf")
#     
#     m <- matrix(c(0.707, 0.707, 0, 0, 
#                   -0.171, 0.181, 0.968, 0, 
#                   0.684, -0.685, 0.249, 0, 
#                   0, 0, 0, 1), byrow = TRUE, nrow = 4)
#     
#     make_HW_plot_panel("model_hw.Rdata", t, col = point_color, 
#                     x = c(2, 0), y = c(2, 5), z = c(2, 10), 
#                     m = m, edges = c("x--", "y+-", "z--"), 
#                     x_label = "", y_label = "", z_label = "", 
#                     plot_file = "hw_univar_2.pdf")
#     
#     make_HW_plot_panel("model_hw.Rdata", t, col = point_color, 
#                     x = c(4, 0), y = c(4, 5), z = c(4, 10), 
#                     m = m, edges = c("x--", "y+-", "z--"), 
#                     x_label = "", y_label = "", z_label = "", 
#                     plot_file = "hw_univar_4.pdf")
# }
# 
#