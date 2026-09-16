library(ggplot2)
source('basic_function.R')

##################################################################################################
##################################################################################################
## Figure 1
##################################################################################################
##################################################################################################

make_triangle_data <- function(id, scale, dx, dy, color_group) {
  base_x <- c(0, 2, 1)
  base_y <- c(0, 0, 3)
  data.frame(x = (base_x * scale) + dx, y = (base_y * scale) + dy, 
             shape_id = id, color_group = color_group, shape_type = "Triangles")
}

make_rect_data <- function(id, scale, dx, dy, color_group) {
  base_x <- c(0, 3, 3, 0)
  base_y <- c(0, 0, 2, 2)
  data.frame(x = (base_x * scale) + dx, y = (base_y * scale) + dy, 
             shape_id = id, color_group = color_group, shape_type = "Rectangles")
}


shapes_data <- rbind(
  # Triangles
  make_triangle_data("T1", scale = 1,   dx = 1, dy = 1, color_group = "Trajectory 1"),
  make_triangle_data("T2", scale = 1.2,   dx = 4, dy = 1, color_group = "Trajectory 2"),
  make_triangle_data("T3", scale = 0.6, dx = 1, dy = 6, color_group = "Trajectory 3"),
  make_triangle_data("T4", scale = 1.5, dx = 5, dy = 5, color_group = "Trajectory 4"),
  
  # Rectangles
  make_rect_data("R1", scale = 1,   dx = 1, dy = 1, color_group = "Trajectory 5"),
  make_rect_data("R2", scale = 1.2, dx = 5, dy = 1, color_group = "Trajectory 6"),
  make_rect_data("R3", scale = 0.6, dx = 2, dy = 5, color_group = "Trajectory 7"),
  make_rect_data("R4", scale = 1.5,   dx = 4, dy = 5, color_group = "Trajectory 8")
)


shapes_data$shape_type <- factor(shapes_data$shape_type, levels = c("Triangles", "Rectangles"))


my_colors <- c(
  "Trajectory 1" = "blue", 
  "Trajectory 2" = "red", 
  "Trajectory 3" = "green", 
  "Trajectory 4" = "brown",
  "Trajectory 5" = "purple",
  "Trajectory 6" = "orange",
  "Trajectory 7" = "cyan",
  "Trajectory 8" = "magenta"
)

my_linetypes <- c(
  "Trajectory 1" = "dashed",   
  "Trajectory 2" = "dotted",   
  "Trajectory 3" = "dotdash",  
  "Trajectory 4" = "solid", 
  "Trajectory 5" = "dashed",    
  "Trajectory 6" = "dotted",
  "Trajectory 7" = "dotdash",
  "Trajectory 8" = "solid"
)


ggplot(shapes_data, aes(x = x, y = y, group = shape_id, color = color_group, linetype = color_group)) +
  geom_polygon(fill = NA, linewidth = 1.2, key_glyph = "path") + 
  scale_color_manual(values = my_colors) + 
  scale_linetype_manual(values = my_linetypes) +  # Maps the new line types
  facet_wrap(~ shape_type) +                 
  coord_fixed() +                            
  theme_void() +                             
  theme(
    legend.position = "bottom",              
    strip.text = element_blank(),            
    legend.key.width = unit(1.5, "cm"),
    panel.spacing = unit(2, "cm") 
  ) +
  # Match the guides so color and linetype merge into a single legend
  guides(
    color = guide_legend(nrow = 2, title = NULL, byrow = TRUE),
    linetype = guide_legend(nrow = 2, title = NULL, byrow = TRUE)
  )

##################################################################################################
##################################################################################################
## Figure 2
##################################################################################################
##################################################################################################

library(plotly)


theta <- seq(0, 2*pi, length.out = 50)
phi <- seq(0, pi, length.out = 50)
sphere_grid <- expand.grid(theta = theta, phi = phi)
x_base <- with(sphere_grid, sin(phi) * cos(theta))
y_base <- with(sphere_grid, sin(phi) * sin(theta))
z_base <- with(sphere_grid, cos(phi))
dim(x_base) <- dim(y_base) <- dim(z_base) <- c(50, 50)


##########################################################################################

align_shape_function <- function(m2,align_mean=c(0,0,1),alpha=10){
  f_mean_new <- intrinsic_mean_sphere(m2)
  shift_log_align <- log_rotation(f_mean_new,align_mean)
  R_translation_align <- exp_from_log_rotation_predict(shift_log_align$L, shift_log_align$theta)
  m2_align = m2 
  
  fc_list = list()
  for (t in 1:nrow(m2)) {
    fc_list[[t]] <- log_rotation(f_mean_new,m2[t,])$L
  }
  norm_f <- norm_CH(fc_list)
  scale_radius <- alpha/norm_f
  for (t in 1:nrow(m2)) {
    ## align
    log_rot <- log_rotation(f_mean_new, m2[t,])
    scaled_L <- scale_radius * log_rot$L
    scaled_theta <- scale_radius * log_rot$theta
    scaled_point_at_origin <- exp_from_log_rotation_predict(scaled_L, scaled_theta) %*% f_mean_new
    m2_align_pt <- R_translation_align %*% scaled_point_at_origin
    m2_align[t,] <- m2_align_pt / sqrt(sum(m2_align_pt^2))
  }
  return(m2_align)
}

same_shape_generation <- function(sphere_matrix, random_radius, random_translation,align_mean=c(0,0,1),alpha=10){
  m2 = m2_align = sphere_matrix
  f_mean <- intrinsic_mean_sphere(sphere_matrix)
  f_mean_new <- (f_mean + random_translation)
  f_mean_new <- f_mean_new / sqrt(sum(f_mean_new^2))
  shift_log <- log_rotation(f_mean, f_mean_new)
  R_translation <- exp_from_log_rotation_predict(shift_log$L, shift_log$theta)
  
  shift_log_align <- log_rotation(f_mean_new,align_mean)
  R_translation_align <- exp_from_log_rotation_predict(shift_log_align$L, shift_log_align$theta)
  
  for (t in 1:nrow(sphere_matrix)) {
    log_rot <- log_rotation(f_mean, sphere_matrix[t,])
    scaled_L <- random_radius * log_rot$L
    scaled_theta <- random_radius * log_rot$theta
    
    if (log_rot$theta < 1e-8) {
      scaled_point_at_origin <- f_mean
    } else {
      scaled_theta <- random_radius * log_rot$theta
      scaled_point_at_origin <- exp_from_log_rotation_predict(scaled_L, scaled_theta) %*% f_mean
    }
    m2_pt <- R_translation %*% scaled_point_at_origin
    m2[t,] <- m2_pt / sqrt(sum(m2_pt^2))
    
  }
  m2_align <- align_shape_function(m2,align_mean,alpha)
  return(list(m2,m2_align))
}


generate_random_rotation <- function() {
  theta_x <- runif(1, 0, pi)
  theta_y <- runif(1, 0, pi)
  theta_z <- runif(1, 0, pi)
  
  R_x <- matrix(c(
    1, 0,            0,
    0, cos(theta_x), -sin(theta_x),
    0, sin(theta_x), cos(theta_x)
  ), nrow = 3, byrow = TRUE)
  
  R_y <- matrix(c(
    cos(theta_y),  0, sin(theta_y),
    0,             1, 0,
    -sin(theta_y), 0, cos(theta_y)
  ), nrow = 3, byrow = TRUE)
  
  R_z <- matrix(c(
    cos(theta_z), -sin(theta_z), 0,
    sin(theta_z),  cos(theta_z), 0,
    0,             0,            1
  ), nrow = 3, byrow = TRUE)
  
  return(R_z %*% R_y %*% R_x)
}


rotate_trajectory <- function(traj_df, R) {
  traj_mat <- as.matrix(traj_df)
  rotated_mat <- traj_mat %*% R
  return(data.frame(
    x = rotated_mat[, 1],
    y = rotated_mat[, 2],
    z = rotated_mat[, 3]
  ))
}

align_and_measure <- function(Z1_df, Z2_df) {
  Z1 <- as.matrix(Z1_df)
  Z2 <- as.matrix(Z2_df)
  
  k <- nrow(Z1)
  m <- ncol(Z1)

  Z1_tilde <- Z1[, 1:(m-1), drop = FALSE]
  Z2_tilde <- Z2[, 1:(m-1), drop = FALSE]

  C <- t(Z2_tilde) %*% Z1_tilde
  
  svd_res <- svd(C)
  U <- svd_res$u
  V <- svd_res$v

  det_UV <- det(U %*% t(V))
  Omega <- diag(c(rep(1, (m - 1) - 1), sign(det_UV)))
  
  R_tilde <- U %*% Omega %*% t(V)
  
  R_star <- diag(m)
  R_star[1:(m-1), 1:(m-1)] <- R_tilde

  Z2_aligned <- Z2 %*% R_star

  Z2_aligned_df <- as.data.frame(Z2_aligned)
  colnames(Z2_aligned_df) <- colnames(Z2_df)
  
  global_trace <- sum(Z1 * Z2_aligned)
  
  mean_inner_product <- global_trace / k
  
  safe_inner_product <- min(max(mean_inner_product, -1), 1)
  
  return(list(
    R_optimal = R_star,
    Z2_aligned = Z2_aligned_df,
    dissimilarity = acos(safe_inner_product)
  ))
}

##########################################################################################

# Trajectory 1: circle

t <- seq(0, 2*pi, length.out = 200)

phi_small <- pi/8
small_circle <- data.frame(
  x = sin(phi_small) * cos(t),
  y = sin(phi_small) * sin(t),
  z = rep(cos(phi_small), 200)
)

m1 <- as.matrix(small_circle)



m2 <- same_shape_generation(m1,0.8,c(-0.1,0,0))[[1]]
m3 <- same_shape_generation(m1,1.2,c(0.3,0.4,0))[[1]]
m4 <- same_shape_generation(m1,1.6,c(0.5,0,0))[[1]]



traj11 <- data.frame(
  x = m1[,1],
  y = m1[,2],
  z = m1[,3] 
)

traj12 <- data.frame(
  x = m2[,1],
  y = m2[,2],
  z = m2[,3] 
)

traj13 <- data.frame(
  x = m3[,1],
  y = m3[,2],
  z = m3[,3] 
)

traj14 <- data.frame(
  x = m4[,1],
  y = m4[,2],
  z = m4[,3] 
)


m1_align <- same_shape_generation(m1,1,c(0,0,0))[[2]]
m2_align <- same_shape_generation(m1,0.8,c(-0.4,0,0))[[2]]
m3_align <- same_shape_generation(m1,1.2,c(0.3,0.4,0))[[2]]
m4_align <- same_shape_generation(m1,1.6,c(0.5,0,0))[[2]]



traj31 <- data.frame(
  x = m1_align[,1],
  y = m1_align[,2],
  z = m1_align[,3] 
)

traj32 <- data.frame(
  x = m2_align[,1],
  y = m2_align[,2],
  z = m2_align[,3] 
)

traj33 <- data.frame(
  x = m3_align[,1],
  y = m3_align[,2],
  z = m3_align[,3] 
)

traj34 <- data.frame(
  x = m4_align[,1],
  y = m4_align[,2],
  z = m4_align[,3] 
)




##############################

t <- seq(0, 2*pi, length.out = 400)
A <- 0.4  
B <- 0.4  

phi_curve <- A * sin(t)
theta_curve <- B * sin(2 * t)

traj21 <- data.frame(
  x = -sin(phi_curve),
  y = cos(phi_curve) * sin(theta_curve),
  z = cos(phi_curve) * cos(theta_curve)
)

m1 <- as.matrix(traj21)
m2 <- same_shape_generation(m1,0.8,c(-0.4,0,0))[[1]]
m3 <- same_shape_generation(m1,1.2,c(0.3,0.4,0))[[1]]
m4 <- same_shape_generation(m1,1.6,c(0.5,0,0))[[1]]


traj22 <- data.frame(
  x = m2[,1],
  y = m2[,2],
  z = m2[,3] 
)

traj23 <- data.frame(
  x = m3[,1],
  y = m3[,2],
  z = m3[,3] 
)

traj24 <- data.frame(
  x = m4[,1],
  y = m4[,2],
  z = m4[,3] 
)

m1_align <- same_shape_generation(m1,1,c(0,0,0),alpha = 1)[[2]]
m2_align <- same_shape_generation(m1,0.8,c(-0.4,0,0),alpha = 1)[[2]]
m3_align <- same_shape_generation(m1,1.2,c(0.3,0.4,0),alpha = 1)[[2]]
m4_align <- same_shape_generation(m1,1.6,c(0.5,0,0),alpha = 1)[[2]]



traj41 <- data.frame(
  x = m1_align[,1],
  y = m1_align[,2],
  z = m1_align[,3] 
)

traj42 <- data.frame(
  x = m2_align[,1],
  y = m2_align[,2],
  z = m2_align[,3] 
)

traj43 <- data.frame(
  x = m3_align[,1],
  y = m3_align[,2],
  z = m3_align[,3] 
)

traj44 <- data.frame(
  x = m4_align[,1],
  y = m4_align[,2],
  z = m4_align[,3] 
)


offset <- 1.5 
no_grid_axis <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE,
  showbackground = FALSE 
)

plot_ly() %>%
  add_surface(x = x_base - (3 * offset), y = y_base, z = z_base, 
            opacity = 0.2, colorscale = "Greys", showscale = FALSE, name = "Sphere 1") %>%
  add_paths(data = traj11, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'blue', dash = "dash", width = 5), name = "Trajectory 1") %>%
  add_paths(data = traj12, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'red', dash = "dot", width = 5), name = "Trajectory 2") %>%
  add_paths(data = traj13, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'green', dash = "dashdot", width = 5), name = "Trajectory 3") %>%
  add_paths(data = traj14, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'brown', dash = "solid", width = 5), name = "Trajectory 4") %>%
  add_surface(x = x_base - offset, y = y_base, z = z_base, 
              opacity = 0.2, colorscale = "Greys", showscale = FALSE, name = "Sphere 2") %>%
  add_paths(data = traj31, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'blue', dash = "dash", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj32, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'red', dash = "dot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj33, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'green', dash = "dashdot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj34, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'brown', dash = "solid", width = 5), showlegend = FALSE) %>%

  add_surface(x = x_base + offset, y = y_base, z = z_base, 
            opacity = 0.2, colorscale = "Blues", showscale = FALSE, name = "Sphere 3") %>%
  add_paths(data = traj21, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'purple', dash = "dash", width = 5), name = "Trajectory 5") %>%
  add_paths(data = traj22, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'orange', dash = "dot", width = 5), name = "Trajectory 6") %>%
  add_paths(data = traj23, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'cyan', dash = "dashdot", width = 5), name = "Trajectory 7") %>%
  add_paths(data = traj24, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'magenta', dash = "solid", width = 5), name = "Trajectory 8") %>%
  add_surface(x = x_base + (3 * offset), y = y_base, z = z_base, 
              opacity = 0.2, colorscale = "Blues", showscale = FALSE, name = "Sphere 4") %>%
  add_paths(data = traj41, x = ~(x + (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'purple', dash = "dash", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj42, x = ~(x + (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'orange', dash = "dot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj43, x = ~(x + (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'cyan', dash = "dashdot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj44, x = ~(x + (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'magenta', dash = "solid", width = 5), showlegend = FALSE) %>%
  layout(
    scene = list(
      xaxis = no_grid_axis,
      yaxis = no_grid_axis,
      zaxis = no_grid_axis,
      aspectmode = "data" 
    ),
    legend = list(
      orientation = "h",       
      xanchor = "center", 
      x = 0.5, 
      yanchor = "top", 
      y = 0.35,               
      itemwidth = 50 
    )
  )

##################################################################################################
##################################################################################################
## Figure 3
##################################################################################################
##################################################################################################


traj21 <- data.frame(
  x = -sin(phi_curve),
  y = cos(phi_curve) * sin(theta_curve),
  z = cos(phi_curve) * cos(theta_curve)
)


m1 <- as.matrix(traj21)
m2 <- same_shape_generation(m1,0.8,c(0,0,0.4))[[1]]
m3 <- same_shape_generation(m1,1.2,c(0,0.4,0.3))[[1]]
m4 <- same_shape_generation(m1,1.6,c(0,0.1,0))[[1]]


traj22 <- data.frame(
  x = m2[,1],
  y = m2[,2],
  z = m2[,3] 
)

traj23 <- data.frame(
  x = m3[,1],
  y = m3[,2],
  z = m3[,3] 
)

traj24 <- data.frame(
  x = m4[,1],
  y = m4[,2],
  z = m4[,3] 
)


set.seed(2)
R1 <- generate_random_rotation()
R2 <- generate_random_rotation()
R3 <- generate_random_rotation()
R4 <- generate_random_rotation()

traj21_rot <- rotate_trajectory(traj21, R1)
traj22_rot <- rotate_trajectory(traj22, R2)
traj23_rot <- rotate_trajectory(traj23, R3)
traj24_rot <- rotate_trajectory(traj24, R4)

traj21_align <- align_shape_function(traj21_rot,alpha = 1)
traj22_align <- align_shape_function(traj22_rot,alpha = 1)
traj23_align <- align_shape_function(traj23_rot,alpha = 1)
traj24_align <- align_shape_function(traj24_rot,alpha = 1)



traj21_final <- traj21_align 
traj22_final <- align_and_measure(traj21_align,traj22_align)[[2]]
traj23_final <- align_and_measure(traj21_align,traj23_align)[[2]]
traj24_final <- align_and_measure(traj21_align,traj24_align)[[2]]

offset <- 1.5 
no_grid_axis <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE,
  showbackground = FALSE 
)

plot_ly() %>%
  add_surface(x = x_base - (3 * offset), y = y_base, z = z_base, 
            opacity = 0.2, colorscale = "Greys", showscale = FALSE, name = "Sphere 1") %>%
  add_paths(data = traj21_rot, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'blue', dash = "dash", width = 5), name = "Trajectory 1") %>%
  add_paths(data = traj22_rot, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'red', dash = "dot", width = 5), name = "Trajectory 2") %>%
  add_paths(data = traj23_rot, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'green', dash = "dashdot", width = 5), name = "Trajectory 3") %>%
  add_paths(data = traj24_rot, x = ~(x - (3 * offset)), y = ~y, z = ~z, 
            line = list(color = 'brown', dash = "solid", width = 5), name = "Trajectory 4") %>%
  
  add_surface(x = x_base - offset, y = y_base, z = z_base, 
              opacity = 0.2, colorscale = "Blues", showscale = FALSE, name = "Sphere 2") %>%
  add_paths(data = traj21_align, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'blue', dash = "dash", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj22_align, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'red', dash = "dot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj23_align, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'green', dash = "dashdot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj24_align, x = ~(x - offset), y = ~y, z = ~z, 
            line = list(color = 'brown', dash = "solid", width = 5), showlegend = FALSE) %>%

  add_surface(x = x_base + offset, y = y_base, z = z_base, 
            opacity = 0.2, colorscale = "Purples", showscale = FALSE, name = "Sphere 3") %>%
  add_paths(data = traj21_final, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'blue', dash = "dash", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj22_final, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'red', dash = "dot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj23_final, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'green', dash = "dashdot", width = 5), showlegend = FALSE) %>%
  add_paths(data = traj24_final, x = ~(x + offset), y = ~y, z = ~z, 
            line = list(color = 'brown', dash = "solid", width = 5), showlegend = FALSE) %>%
  layout(
    scene = list(
      xaxis = no_grid_axis,
      yaxis = no_grid_axis,
      zaxis = no_grid_axis,
      aspectmode = "data" 
    ),
    legend = list(
      orientation = "h",       
      xanchor = "center", 
      x = 0.5, 
      yanchor = "top", 
      y = 0.25,               
      itemwidth = 50 
    )
  )
