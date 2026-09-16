### sphere
sphere_ts_dissimilarity_function <- function(f,g){
  dot_products <- rowSums(f * g)
  safe_dot_products <- pmin(pmax(dot_products, -1), 1)
  pointwise_distances <- acos(safe_dot_products)
  return(sum(pointwise_distances))
}

product_sphere_ts_dissimilarity_function <- function(f,g){
  inner_prod_cum <- 0
  for (k in 1:length(f)) {
    dot_products <- rowSums(f[[k]] * g[[k]])
    safe_dot_products <- pmin(pmax(dot_products, -1), 1)
    pointwise_distances <- acos(safe_dot_products)
    inner_prod_cum <- inner_prod_cum + sum(pointwise_distances)
  }
  
  return(inner_prod_cum)
}


### function
multi_function_dissimilarity_function <- function(f,g){
  measure <- sum((f-g)^2)
  return(measure)
}

product_function_dissimilarity_function <- function(f,g){
  inner_prod_cum <- 0
  for (k in 1:length(f)) {
    inner_prod_cum <- inner_prod_cum + sum((f[[k]]-g[[k]])^2)
  }
  return(inner_prod_cum)
}

function_shape_dissimilarity_function_1D <- function(f,g){
  f_mean <- mean(f)
  g_mean <- mean(g)
  
  fc <- f - f_mean
  gc <- g - g_mean

  norm_f <- sqrt(sum(fc^2))
  norm_g <- sqrt(sum(gc^2))

  tol <- 1e-12

  if (norm_f < tol && norm_g < tol) {
    inner_prod_cum <- inner_prod_cum + 1
  } else if ((norm_f < tol && norm_g >= tol) || (norm_f >= tol && norm_g < tol)) {
    inner_prod_cum <- inner_prod_cum + 0
  } else {
    f_cn <- fc / norm_f
    g_cn <- gc / norm_g
    inner_prod <- sum(f_cn*g_cn)
    inner_prod_clamped <- min(max(inner_prod, -1), 1)
  }
  return(acos(inner_prod_clamped))
}


product_function_shape_dissimilarity_function_1D <- function(f,g){
  inner_prod_cum <- 0
  for (k in 1:length(f)) {
    f_mean <- mean(f[[k]])
    g_mean <- mean(g[[k]])
    
    fc <- f[[k]] - f_mean
    gc <- g[[k]] - g_mean

    norm_f <- sqrt(sum(fc^2))
    norm_g <- sqrt(sum(gc^2))
    tol <- 1e-12
    if (norm_f < tol && norm_g < tol) {
      inner_prod_cum <- inner_prod_cum + 1
    } else if ((norm_f < tol && norm_g >= tol) || (norm_f >= tol && norm_g < tol)) {
      inner_prod_cum <- inner_prod_cum + 0
    } else {
      f_cn <- fc / norm_f
      g_cn <- gc / norm_g
      inner_prod <- sum(f_cn*g_cn)
      inner_prod_clamped <- min(max(inner_prod, -1), 1)
      inner_prod_cum <- inner_prod_cum + inner_prod_clamped
    }
  }
  return(acos(inner_prod_cum/length(f)))
}


# -------------------------------------------------------------------
# Procrustes Shape Distance Function
# -------------------------------------------------------------------
procrustes_distance <- function(Z1_raw, Z2_raw) {
  k1 <- nrow(Z1_raw)
  k2 <- nrow(Z2_raw)
  
  C1 <- diag(k1) - matrix(1, k1, k1) / k1
  C2 <- diag(k2) - matrix(1, k2, k2) / k2
  
  Z1_centered <- C1 %*% Z1_raw
  Z2_centered <- C2 %*% Z2_raw

  Z1_preshape <- Z1_centered / sqrt(sum(Z1_centered^2))
  Z2_preshape <- Z2_centered / sqrt(sum(Z2_centered^2))

  m <- ncol(Z2_preshape)
  C_mat <- t(Z2_preshape) %*% Z1_preshape
  
  svd_res <- svd(C_mat)
  U <- svd_res$u
  V <- svd_res$v

  det_UV <- det(U %*% t(V))
  Omega <- diag(c(rep(1, m - 1), sign(det_UV)))
  R_optimal <- U %*% Omega %*% t(V)

  Z2_aligned <- Z2_preshape %*% R_optimal
  shape_distance <- sqrt(sum((Z1_preshape - Z2_aligned)^2))
  
  return(shape_distance)
}


### sphere
span_basis_xy <- function(x, y) {
  u1 <- unit_norm(x)
  proj <- sum(y * u1) * u1
  v <- y - proj
  nv <- sqrt(sum(v^2))
  if (nv < 1e-12) {
    tmp <- rnorm(length(x))
    tmp <- tmp - sum(tmp * u1) * u1
    u2 <- unit_norm(tmp)
  } else {
    u2 <- v / nv
  }
  list(u1 = u1, u2 = u2)
}


log_rotation <- function(x, y, tol = 1e-12) { ## move x to y, mimic y-x
  x <- as.numeric(x); y <- as.numeric(y)
  stopifnot(length(x) == length(y))
  x <- x / sqrt(sum(x^2))
  y <- y / sqrt(sum(y^2))
  
  cxy <- sum(x * y)
  cxy <- max(min(cxy, 1), -1)
  if (1 - abs(cxy) < tol) {
    if (cxy > 0) {
      return(list(L = matrix(0, length(x), length(x)), theta = 0, u1 = x, u2 = NA))
    } else {
      j <- which.min(abs(x))
      e <- rep(0, length(x)); e[j] <- 1
      v <- e - sum(e * x) * x
      vn <- sqrt(sum(v^2))
      if (vn < tol) {
        if (length(x) < 2) stop("Need dimension >= 2 for antipodal case")
        v <- c(-x[2], x[1], if (length(x) > 2) rep(0, length(x) - 2) else NULL)
        vn <- sqrt(sum(v^2))
      }
      u2 <- v / vn
      theta <- pi
      Q <- tcrossprod(u2, x) - tcrossprod(x, u2)
      return(list(L = theta * Q, theta = theta, u1 = x, u2 = u2))
    }
  }

  u1 <- x
  v <- y - cxy * u1
  u2 <- v / sqrt(sum(v^2))
  theta <- acos(cxy)
  Q <- tcrossprod(u2, u1) - tcrossprod(u1, u2)
  
  list(L = theta * Q, theta = theta, u1 = u1, u2 = u2)
}



#### intrinsic mean
intrinsic_mean_sphere <- function(y, tol = 1e-8, max_iter = 1000) {
  n <- nrow(y)
  d <- ncol(y)

  mu <- colSums(y)
  mu <- mu / sqrt(sum(mu^2))
  
  for (iter in 1:max_iter) {
    tangent_sum <- rep(0, d)
    for (i in 1:n) {
      theta <- acos(pmin(pmax(sum(mu * y[i,]), -1), 1)) 
      if (theta > 1e-12) {
        tangent <- (theta / sin(theta)) * (y[i,] - cos(theta) * mu)
        tangent_sum <- tangent_sum + tangent
      }
    }
    tangent_mean <- tangent_sum / n
    norm_t <- sqrt(sum(tangent_mean^2))

    if (norm_t < tol) break
    mu <- cos(norm_t) * mu + sin(norm_t) * (tangent_mean / norm_t)
  }
  
  return(mu)
}

inner_product_CH <- function(list_A, list_B) {
  mean(mapply(function(mat_A, mat_B) sum(mat_A * mat_B), list_A, list_B))
}


norm_CH <- function(mat_list) {
  sqrt(inner_product_CH(mat_list, mat_list))
}

shape_dissimilarity_function_1D_ref_inv <- function(f,g){
  f_mean <- intrinsic_mean_sphere(f)
  g_mean <- intrinsic_mean_sphere(g)
  fc_list = gc_list = list()
  for (t in 1:nrow(f)) {
    fc_list[[t]] <- log_rotation(f_mean,f[t,])$L
    gc_list[[t]] <- log_rotation(g_mean,g[t,])$L
  }

  norm_f <- norm_CH(fc_list)
  norm_g <- norm_CH(gc_list)

  tol <- 1e-12

  if (norm_f < tol && norm_g < tol) {
    return(0)
  } else if ((norm_f < tol && norm_g >= tol) || (norm_f >= tol && norm_g < tol)) {
    return(pi / 2)
  } else {
    f_cn <- lapply(fc_list, function(mat) mat / norm_f)
    g_cn <- lapply(gc_list, function(mat) mat / norm_g)

    inner_prod <- inner_product_CH(f_cn, g_cn)
    inner_prod_clamped <- min(max(inner_prod, -1), 1)

    return(acos(abs(inner_prod_clamped)))
  }
  
}


product_shape_dissimilarity_function_1D <- function(f,g){
  
  inner_prod_cum <- 0
  for (k in 1:length(f)) {
    f_mean <- intrinsic_mean_sphere(f[[k]])
    g_mean <- intrinsic_mean_sphere(g[[k]])
    fc_list = gc_list = list()
    for (t in 1:nrow(f[[k]])) {
      fc_list[[t]] <- log_rotation(f_mean,f[[k]][t,])$L
      gc_list[[t]] <- log_rotation(g_mean,g[[k]][t,])$L
    }
    
    norm_f <- norm_CH(fc_list)
    norm_g <- norm_CH(gc_list)
    tol <- 1e-12

    if (norm_f < tol && norm_g < tol) {
      inner_prod_cum <- inner_prod_cum + 1
    } else if ((norm_f < tol && norm_g >= tol) || (norm_f >= tol && norm_g < tol)) {
      inner_prod_cum <- inner_prod_cum + 0
    } else {
      f_cn <- lapply(fc_list, function(mat) mat / norm_f)
      g_cn <- lapply(gc_list, function(mat) mat / norm_g)

      inner_prod <- inner_product_CH(f_cn, g_cn)

      inner_prod_clamped <- min(max(inner_prod, -1), 1)
      inner_prod_cum <- inner_prod_cum + inner_prod_clamped
    }
  }
  return(acos(inner_prod_cum/length(f)))
}

shape_dissimilarity_function_1D <- function(f,g){
  f_mean <- intrinsic_mean_sphere(f)
  g_mean <- intrinsic_mean_sphere(g)
  fc_list = gc_list = list()
  for (t in 1:nrow(f)) {
    fc_list[[t]] <- log_rotation(f_mean,f[t,])$L
    gc_list[[t]] <- log_rotation(g_mean,g[t,])$L
  }

  norm_f <- norm_CH(fc_list)
  norm_g <- norm_CH(gc_list)

  tol <- 1e-12

  if (norm_f < tol && norm_g < tol) {
    return(0)
  } else if ((norm_f < tol && norm_g >= tol) || (norm_f >= tol && norm_g < tol)) {
    return(pi / 2)
  } else {
    f_cn <- lapply(fc_list, function(mat) mat / norm_f)
    g_cn <- lapply(gc_list, function(mat) mat / norm_g)
    inner_prod <- inner_product_CH(f_cn, g_cn)
    inner_prod_clamped <- min(max(inner_prod, -1), 1)
    return(acos(inner_prod_clamped))
  }
  
}



shape_dissimilarity_function <- function(f, g, max_iter = 100, step_size = 0.1, tol = 1e-12, grad_tol = 1e-6) {
  f_mean <- intrinsic_mean_sphere(f)
  g_mean <- intrinsic_mean_sphere(g)
  
  fc_list <- list()
  gc_list <- list()
  for (t in 1:nrow(f)) {
    fc_list[[t]] <- log_rotation(f_mean, f[t,])$L
    gc_list[[t]] <- log_rotation(g_mean, g[t,])$L
  }

  norm_f <- norm_CH(fc_list)
  norm_g <- norm_CH(gc_list)
  
  if (norm_f < tol && norm_g < tol) {
    return(0)
  } else if ((norm_f < tol && norm_g >= tol) || (norm_f >= tol && norm_g < tol)) {
    return(pi / 2)
  } else {
    f_cn <- lapply(fc_list, function(mat) mat / norm_f)
    g_cn <- lapply(gc_list, function(mat) mat / norm_g)

    n <- nrow(f_cn[[1]])
    R_opt <- diag(n) 
    
    for (iter in 1:max_iter) {
      grad_E <- matrix(0, nrow = n, ncol = n)
      for (t in 1:length(f_cn)) {
        grad_E <- grad_E - 2 * (f_cn[[t]] %*% R_opt %*% g_cn[[t]])
      }

      grad_R <- R_opt %*% (0.5 * (t(R_opt) %*% grad_E - t(grad_E) %*% R_opt))

      if (max(abs(grad_R)) < grad_tol) {
        break
      }

      R_step <- R_opt + step_size * grad_R
      svd_res <- svd(R_step)

      UV_trans <- svd_res$u %*% t(svd_res$v)
      det_sign <- sign(det(UV_trans))
      
      D <- diag(1, n)
      D[n, n] <- det_sign

      R_opt <- svd_res$u %*% D %*% t(svd_res$v)
    }
    
    supremum_inner_prod <- 0
    for (t in 1:length(f_cn)) {
      transformed_g <- R_opt %*% g_cn[[t]] %*% t(R_opt)
      supremum_inner_prod <- supremum_inner_prod + sum(diag(t(f_cn[[t]]) %*% transformed_g))
    }

    inner_prod_clamped <- min(max(supremum_inner_prod, -1), 1)
    
    return(acos(inner_prod_clamped))
  }
}



transport_shape_function <- function(m2, align_mean) {
  f_mean <- intrinsic_mean_sphere(m2)
  shift_log <- log_rotation(f_mean, align_mean)
  
  R_align <- exp_from_log_rotation_predict(
    shift_log$L,
    shift_log$theta
  )
  m2_transport <- m2
  
  for (t in seq_len(nrow(m2))) {
    
    transported_point <- R_align %*% m2[t, ]
    m2_transport[t, ] <-
      transported_point / sqrt(sum(transported_point^2))
  }
  
  return(list(
    trajectory = m2_transport,
    original_mean = f_mean,
    target_mean = align_mean,
    R = R_align
  ))
}

align_shape_function <- function(m2,align_mean,alpha=10){
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
    log_rot <- log_rotation(f_mean_new, m2[t,])
    scaled_L <- scale_radius * log_rot$L
    scaled_theta <- scale_radius * log_rot$theta
    scaled_point_at_origin <- exp_from_log_rotation_predict(scaled_L, scaled_theta) %*% f_mean_new
    m2_align_pt <- R_translation_align %*% scaled_point_at_origin

    m2_align[t,] <- m2_align_pt / sqrt(sum(m2_align_pt^2))
  }
  return(m2_align)
}

# -------------------------------------------------------------------
# Spherical Procrustes Shape Distance Function
# -------------------------------------------------------------------
sphere_procrustes_distance <- function(Z1_raw, Z2_raw) {
  m <- ncol(Z1_raw)
  k <- nrow(Z2_raw)
  
  align_mean <- c(rep(0,m-1),1)
  Z1 <- align_shape_function(m2=Z1_raw,align_mean=align_mean,alpha=10)
  Z2 <- align_shape_function(m2=Z2_raw,align_mean=align_mean,alpha=10)
  
  Y1_tilde <- Z1[, 1:(m-1), drop = FALSE]
  Y2_tilde <- Z2[, 1:(m-1), drop = FALSE]

  C <- t(Y2_tilde) %*% Y1_tilde

  svd_res <- svd(C)
  U <- svd_res$u
  V <- svd_res$v

  det_UV <- det(U %*% t(V))
  Omega <- diag(c(rep(1, (m - 1) - 1), sign(det_UV)))
  
  R_tilde <- U %*% Omega %*% t(V)
  
  R_star <- diag(m)
  R_star[1:(m-1), 1:(m-1)] <- R_tilde
  Z2_aligned <- Z2 %*% R_star
  
  dot_products <- rowSums(Z1 * Z2_aligned)
  
  safe_dot_products <- pmin(pmax(dot_products, -1), 1)
 
  pointwise_distances <- acos(safe_dot_products)
  return(sum(pointwise_distances))
}


product_spherical_procrustes <- function(
    Z1_list,
    Z2_list,
    alpha = 10
) {

  if (!is.list(Z1_list) || !is.list(Z2_list)) {
    stop("Z1_list and Z2_list must be lists.")
  }
  
  L <- length(Z1_list)
  
  if (L == 0) {
    stop("The lists must contain at least one component.")
  }
  
  if (length(Z2_list) != L) {
    stop("Z1_list and Z2_list must contain the same number of components.")
  }

  R_list <- vector("list", L)
  Z1_align_list <- vector("list", L)
  Z2_align_list <- vector("list", L)
  
  component_trace <- numeric(L)

  for (l in 1:L) {
    
    Z1_raw <- Z1_list[[l]]
    Z2_raw <- Z2_list[[l]]
  
    if (!is.matrix(Z1_raw) || !is.matrix(Z2_raw)) {
      stop(
        paste0(
          "Component ", l,
          " must contain matrices."
        )
      )
    }
    
    if (!all(dim(Z1_raw) == dim(Z2_raw))) {
      stop(
        paste0(
          "Components ", l,
          " in Z1_list and Z2_list must have the same dimensions."
        )
      )
    }
    
    K <- nrow(Z1_raw)
    m <- ncol(Z1_raw)
    
    if (m < 2) {
      stop("The number of spherical dimensions m must be at least 2.")
    }
    
    align_mean <- c(rep(0, m - 1), 1)
    
    Z1 <- align_shape_function(
      m2 = Z1_raw,
      align_mean = align_mean,
      alpha = alpha
    )
    
    Z2 <- align_shape_function(
      m2 = Z2_raw,
      align_mean = align_mean,
      alpha = alpha
    )

    Z1_align_list[[l]] <- Z1
    Z2_align_list[[l]] <- Z2

    Y1_tilde <- Z1[, 1:(m - 1), drop = FALSE]
    Y2_tilde <- Z2[, 1:(m - 1), drop = FALSE]

    C <- t(Y2_tilde) %*% Y1_tilde

    svd_res <- svd(C)
    
    U <- svd_res$u
    V <- svd_res$v
    det_UV <- det(U %*% t(V))
    
    Omega <- diag(m - 1)
    
    if (det_UV < 0) {
      Omega[m - 1, m - 1] <- -1
    }
    
    R_tilde <- U %*% Omega %*% t(V)
    R_star <- diag(m)
    
    R_star[
      1:(m - 1),
      1:(m - 1)
    ] <- R_tilde
    
    R_list[[l]] <- R_star
    Z2_aligned <- Z2 %*% R_star
    row_norms <- sqrt(rowSums(Z2_aligned^2))
    
    Z2_aligned <- Z2_aligned / row_norms
    Z2_align_list[[l]] <- Z2_aligned
    dot_products <- rowSums(
      Z1 * Z2_aligned
    )
    
    # Protect against numerical error
    dot_products <- pmin(
      pmax(dot_products, -1),
      1
    )
    
    component_trace[l] <- sum(dot_products)
  }

  K <- nrow(Z1_list[[1]])
  
  mean_trace <- sum(component_trace) / (K * L)
  mean_trace <- pmin(
    pmax(mean_trace, -1),
    1
  )

  distance <- acos(mean_trace)
  
  return(
    list(
      distance = distance,
      rotations = R_list,
      Z1_aligned = Z1_align_list,
      Z2_aligned = Z2_align_list,
      component_trace = component_trace,
      mean_trace = mean_trace
    )
  )
}

D_rotProduct <- function(Z1_list, Z2_list) {
  
  result <- product_spherical_procrustes(
    Z1_list = Z1_list,
    Z2_list = Z2_list
  )
  
  return(result$distance)
}


SpheGeoDist <- function(y1,y2) {
  if (abs(length(y1) - length(y2)) > 0) {
    stop("y1 and y2 should be of the same length.")
  }
  if ( !isTRUE( all.equal(l2norm(y1),1) ) ) {
    stop("y1 is not a unit vector.")
  }
  if ( !isTRUE( all.equal(l2norm(y2),1) ) ) {
    stop("y2 is not a unit vector.")
  }
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  if (sum(y1 * y2) > 1){
    return(0)
  } else if (sum(y1*y2) < -1){
    return(pi)
  } else return(acos(sum(y1 * y2)))
}

unit_norm <- function(x) x / sqrt(sum(x^2))



geo_dist <- function(x, y) {
  # Geodesic distance on the unit sphere (radians)
  z <- sum(x * y)
  z <- max(min(z, 1), -1)
  acos(z)
}

l2norm <- function(x){
  as.numeric(sqrt(crossprod(x)))
}

SpheGeoGrad <- function(x, y, tol=1e-8) {
  dot <- sum(x * y)
  dot <- max(min(dot, 1), -1)  # clamp
  tmp <- 1 - dot^2
  if (tmp < tol) {
    return(rep(0, length(x)))
  }
  
  return(-(1/sqrt(tmp)) * x)
}

SpheGeoHess <- function(x, y, tol=1e-8) {
  dot <- sum(x * y)
  dot <- max(min(dot, 1), -1)
  tmp <- 1 - dot^2
  p <- length(x)

  if (tmp < tol) {
    return(2 * (diag(p) - outer(y, y)))
  }
  
  return(- dot * (tmp)^(-1.5) * (x %*% t(x)))
}

# set up bandwidth range
SetBwRange <- function(xin, xout, kernel_type) {
  xinSt <- unique(sort(xin))
  bw.min <- max(diff(xinSt), xinSt[2] - min(xout), max(xout) -
                  xinSt[length(xin)-1])*1.1 / (ifelse(kernel_type == "gauss", 3, 1) *
                                                 ifelse(kernel_type == "gausvar", 2.5, 1))
  bw.max <- diff(range(xin))/3
  if (bw.max < bw.min) {
    if (bw.min > bw.max*3/2) {
      #warning("Data is too sparse.")
      bw.max <- bw.min*1.01
    } else bw.max <- bw.max*3/2
  }
  return(list(min=bw.min, max = bw.max))
}

kerFctn <- function(kernel_type){
  if (kernel_type=='gauss'){
    ker <- function(x){
      dnorm(x) #exp(-x^2 / 2) / sqrt(2*pi)
    }
  } else if(kernel_type=='rect'){
    ker <- function(x){
      as.numeric((x<=1) & (x>=-1))
    }
  } else if(kernel_type=='epan'){
    ker <- function(x){
      n <- 1
      (2*n+1) / (4*n) * (1-x^(2*n)) * (abs(x)<=1)
    }
  } else if(kernel_type=='gausvar'){
    ker <- function(x) {
      dnorm(x)*(1.25-0.25*x^2)
    }
  } else if(kernel_type=='quar'){
    ker <- function(x) {
      (15/16)*(1-x^2)^2 * (abs(x)<=1)
    }
  } else {
    stop('Unavailable kernel')
  }
  return(ker)
}

LocSpheReg <- function(xin=NULL, yin=NULL, xout=NULL, optns=list()){
  
  if (is.null(xin))
    stop ("xin has no default and must be input by users.")
  if (is.null(yin))
    stop ("yin has no default and must be input by users.")
  if (is.null(xout))
    xout <- xin
  if (!is.vector(xin) | !is.numeric(xin))
    stop("xin should be a numerical vector.")
  if (!is.matrix(yin) | !is.numeric(yin))
    stop("yin should be a numerical matrix.")
  if (!is.vector(xout) | !is.numeric(xout))
    stop("xout should be a numerical vector.")
  if (length(xin)!=nrow(yin))
    stop("The length of xin should be the same as the number of rows in yin.")
  if (sum(abs(rowSums(yin^2) - rep(1,nrow(yin))) > 1e-6)){
    yin = yin / sqrt(rowSums(yin^2))
    warning("Each row of yin has been standardized to enforce sum of squares equal to 1.")
  }
  
  if (is.null(optns$bw)){
    optns$bw <- "CV" #max(sort(xin)[-1] - sort(xin)[-length(xin)]) * 1.2
  }
  if (is.character(optns$bw)) {
    if (optns$bw != "CV") {
      warning("Incorrect input for optns$bw.")
    }
  } else if (!is.numeric(optns$bw)) {
    stop("Mis-specified optns$bw.")
  }
  if (length(optns$bw) > 1)
    stop("bw should be of length 1.")
  
  if (is.null(optns$kernel))
    optns$kernel <- "gauss"
  
  if (is.numeric(optns$bw)) {
    bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$kernel)
    if (optns$bw < bwRange$min | optns$bw > bwRange$max) {
      optns$bw <- "CV"
      warning("optns$bw is too small or too large; reset to be chosen by CV.")
    }
  } 
  if (optns$bw == "CV") {
    optns$bw <- bwCV_sphe_moving_window(xin = xin, yin = yin, xout = xout, optns = optns)
  }
  yout <- LocSpheGeoReg(xin = xin, yin = yin, xout = xout, optns = optns)
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin, optns = optns)
  class(res) <- "spheReg"
  return(res)
}

LocSpheGeoReg <- function(xin, yin, xout, optns = list()) {
  k = length(xout)
  n = length(xin)
  m = ncol(yin)
  
  bw <- optns$bw
  ker <- kerFctn(optns$kernel)
  
  yout = sapply(1:k, function(j){
    mu0 = mean(ker((xout[j] - xin) / bw))
    mu1 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j]))
    mu2 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j])^2)
    s = ker((xout[j] - xin) / bw) * (mu2 - mu1 * (xin - xout[j])) /
      (mu0 * mu2 - mu1^2)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if (sum(sapply(1:n, function(i) sum(yin[i,]*y0))[ker((xout[j] - xin) / bw)>0] > 1-1e-8)){
      #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) +
      #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      # return(y0)
      ##y0 = y0 + rnorm(3) * 1e-3
      y0 = y0 + rnorm(m) * 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      # y <- y / l2norm(y)
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1e5)
    # res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument / l2norm(res$argument))
  })
  return(t(yout))
}


bwCV_sphe <- function(xin, yin, xout, optns) {
  yin <- yin[order(xin),]
  xin <- sort(xin)
  compareRange <- (xin > min(xin) + diff(range(xin))/5) & (xin < max(xin) - diff(range(xin))/5)
  
  # k-fold
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw
    folds <- numeric(length(xin))
    n <- sum(compareRange)
    numFolds <- ifelse(n > 30, 10, sum(compareRange))
    
    tmp <- c(sapply(1:ceiling(n/numFolds), function(i)
      sample(x = seq_len(numFolds), size = numFolds, replace = FALSE)))
    tmp <- tmp[1:n]
    repIdx <- which(diff(tmp) == 0)
    for (i in which(diff(tmp) == 0)) {
      s <- tmp[i]
      tmp[i] <- tmp[i-1]
      tmp[i-1] <- s
    }
    #tmp <- cut(1:n,breaks = seq(0,n,length.out = numFolds+1), labels=FALSE)
    #tmp <- tmp[sample(seq_len(n), n)]
    
    folds[compareRange] <- tmp
    
    yout <- lapply(seq_len(numFolds), function(foldidx) {
      testidx <- which(folds == foldidx)
      res <- LocSpheGeoReg(xin = xin[-testidx], yin = yin[-testidx,], xout = xin[testidx], optns = optns1)
      res # each row is a spherical vector
    })
    yout <- do.call(rbind, yout)
    yinMatch <- yin[which(compareRange)[order(tmp)],]
    mean(sapply(1:nrow(yout), function(i) SpheGeoDist(yout[i,], yinMatch[i,])^2))
  }
  bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)
  #if (!is.null(optns$bwRange)) {
  #  if (min(optns$bwRange) < bwRange$min) {
  #    message("Minimum bandwidth is too small and has been reset.")
  #  } else bwRange$min <- min(optns$bwRange)
  #  if (max(optns$bwRange) >  bwRange$min) {
  #    bwRange$max <- max(optns$bwRange)
  #  } else {
  #    message("Maximum bandwidth is too small and has been reset.")
  #  }
  #}
  res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
  res$minimum
}

bwCV_sphe_moving_window <- function(xin, yin, xout, optns) {
  # 1. Ensure data is strictly ordered by time
  ord <- order(xin)
  yin <- yin[ord, ]
  xin <- xin[ord]
  
  n <- length(xin)

  window_size <- ceiling(n * 0.8) 
  
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw

    test_indices <- (window_size + 1):n
    
    errors <- sapply(test_indices, function(i) {
      train_idx <- (i - window_size):(i - 1)
      test_idx <- i
      res <- LocSpheGeoReg(
        xin = xin[train_idx], 
        yin = yin[train_idx, , drop = FALSE], 
        xout = xin[test_idx], 
        optns = optns1
      )

      geo_dist(res, yin[test_idx, ])^2
    })

    mean(errors)
  }

  bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)

  res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
  
  return(res$minimum)
}


knn_from_dissimilarity <- function(dissimilarity_matrix, train_items, test_items, train_labels, k = 3) {

  if(is.null(names(train_labels))) {
    names(train_labels) <- train_items
  }

  D_test_train <- dissimilarity_matrix[test_items, train_items, drop = FALSE]

  predictions <- character(length(test_items))
  names(predictions) <- test_items

  for (i in 1:nrow(D_test_train)) {
    dists <- D_test_train[i, ]
    nearest_neighbors <- names(sort(dists))[1:k]
    nn_labels <- train_labels[nearest_neighbors]
    label_counts <- table(nn_labels)
    predicted_label <- names(label_counts)[which.max(label_counts)]
    predictions[i] <- predicted_label
  }
  
  return(predictions)
}

tune_knn_loocv <- function(dissimilarity_matrix, all_items, all_labels, k_values = seq(1, 21, by = 2)) {
  accuracies <- numeric(length(k_values))
  
  for (idx in seq_along(k_values)) {
    current_k <- k_values[idx]
    loocv_predictions <- character(length(all_items))
    names(loocv_predictions) <- all_items
    for (current_test_item in all_items) {
      current_train_items <- setdiff(all_items, current_test_item)
      current_train_labels <- all_labels[current_train_items]
      
      pred <- knn_from_dissimilarity(
        dissimilarity_matrix = dissimilarity_matrix,
        train_items = current_train_items,
        test_items = current_test_item,   
        train_labels = current_train_labels,
        k = current_k
      )
      
      loocv_predictions[current_test_item] <- pred
    }
    cm <- table(Predicted = loocv_predictions, Actual = all_labels)
    accuracies[idx] <- sum(diag(cm)) / sum(cm)
  }
  results_df <- data.frame(
    k = k_values,
    accuracy = accuracies
  )
  
  return(results_df)
}


classification_evaluation_function <- function(dissimilarity_matrix,
                                               all_items,
                                               train_items,
                                               train_labels,
                                               test_items,
                                               true_test_labels) {
  
  max_k <- max(21,round(sqrt(nrow(dissimilarity_matrix))))
  train_dissimilarity_matrix <- dissimilarity_matrix[train_items, train_items, drop = FALSE]
  
  tuning_results <- tune_knn_loocv(
    dissimilarity_matrix = train_dissimilarity_matrix, 
    all_items = train_items, 
    all_labels = train_labels,
    k_values = seq(1, max_k, by = 2)
  )

  k_choose <- tuning_results$k[which.max(tuning_results$accuracy)]

  predictions <- knn_from_dissimilarity(
    dissimilarity_matrix = dissimilarity_matrix,
    train_items = train_items,
    test_items = test_items,
    train_labels = train_labels,
    k = k_choose 
  )
  
  confusion_matrix_random <- table(Predicted = predictions, Actual = true_test_labels)
  accuracy_random <- sum(diag(confusion_matrix_random)) / sum(confusion_matrix_random)

  
  loocv_predictions <- character(length(all_items))
  names(loocv_predictions) <- all_items

  for (current_test_item in all_items) {
    current_train_items <- setdiff(all_items, current_test_item)
    current_train_labels <- all_labels[current_train_items]
    pred <- knn_from_dissimilarity(
      dissimilarity_matrix = dissimilarity_matrix,
      train_items = current_train_items,
      test_items = current_test_item,   
      train_labels = current_train_labels,
      k = k_choose
    )
    
    loocv_predictions[current_test_item] <- pred
  }
  

  confusion_matrix_LOOCV <- table(Predicted = loocv_predictions, Actual = all_labels)
  accuracy_LOOCV <- sum(diag(confusion_matrix_LOOCV)) / sum(confusion_matrix_LOOCV)
  
  return(list(confusion_matrix_random,accuracy_random,confusion_matrix_LOOCV,accuracy_LOOCV))
}


classification_evaluation_function <- function(
    dissimilarity_matrix,
    all_items,
    all_labels,
    train_items,
    train_labels,
    test_items,
    true_test_labels,
    max_k_cap = 21
) {
  
  # ------------------------------------------------------------
  # Helper: generate valid odd k values
  # ------------------------------------------------------------
  get_k_values <- function(n_train, max_k_cap = 21) {
    
    # During LOOCV tuning, each inner training fold contains
    # n_train - 1 observations, so k cannot exceed n_train - 1.
    max_k <- min(max_k_cap, n_train - 1)
    
    if (max_k < 1) {
      stop("Not enough observations to perform kNN tuning.")
    }
    
    # Use odd k values only
    if (max_k %% 2 == 0) {
      max_k <- max_k - 1
    }
    
    seq(1, max_k, by = 2)
  }
  
  
  # ============================================================
  # 1. TRAIN / TEST EVALUATION
  #    Tune k using only the training data
  # ============================================================
  
  train_dissimilarity_matrix <- dissimilarity_matrix[
    train_items,
    train_items,
    drop = FALSE
  ]
  
  k_values <- get_k_values(
    n_train = length(train_items),
    max_k_cap = max_k_cap
  )
  
  tuning_results <- tune_knn_loocv(
    dissimilarity_matrix = train_dissimilarity_matrix,
    all_items = train_items,
    all_labels = train_labels,
    k_values = k_values
  )
  
  # Choose k giving the highest LOOCV accuracy
  # If there is a tie, choose the smallest k
  best_accuracy <- max(tuning_results$accuracy)
  
  k_choose <- min(
    tuning_results$k[
      tuning_results$accuracy == best_accuracy
    ]
  )
  
  
  # Predict the held-out test set
  predictions <- knn_from_dissimilarity(
    dissimilarity_matrix = dissimilarity_matrix,
    train_items = train_items,
    test_items = test_items,
    train_labels = train_labels,
    k = k_choose
  )
  
  confusion_matrix_random <- table(
    Predicted = predictions,
    Actual = true_test_labels
  )
  
  accuracy_random <- sum(diag(confusion_matrix_random)) /
    sum(confusion_matrix_random)
  
  
  # ============================================================
  # 2. NESTED LOOCV ON ALL ITEMS
  #
  #    Outer loop:
  #        leave one observation out
  #
  #    Inner loop:
  #        tune k using ONLY the remaining observations
  #
  #    Therefore, the held-out observation plays no role in
  #    selecting k.
  # ============================================================
  
  loocv_predictions <- character(length(all_items))
  names(loocv_predictions) <- all_items
  
  loocv_selected_k <- integer(length(all_items))
  names(loocv_selected_k) <- all_items
  
  
  for (current_test_item in all_items) {
    
    # ----------------------------------------------------------
    # Outer LOOCV split
    # ----------------------------------------------------------
    current_train_items <- setdiff(
      all_items,
      current_test_item
    )
    
    current_train_labels <- all_labels[
      current_train_items
    ]
    
    
    # ----------------------------------------------------------
    # Dissimilarity matrix containing only outer-training data
    # ----------------------------------------------------------
    current_train_dissimilarity_matrix <-
      dissimilarity_matrix[
        current_train_items,
        current_train_items,
        drop = FALSE
      ]
    
    
    # ----------------------------------------------------------
    # Candidate k values for inner LOOCV
    # ----------------------------------------------------------
    current_k_values <- get_k_values(
      n_train = length(current_train_items),
      max_k_cap = max_k_cap
    )
    
    
    # ----------------------------------------------------------
    # Inner LOOCV:
    # tune k using outer-training observations only
    # ----------------------------------------------------------
    current_tuning_results <- tune_knn_loocv(
      dissimilarity_matrix =
        current_train_dissimilarity_matrix,
      all_items = current_train_items,
      all_labels = current_train_labels,
      k_values = current_k_values
    )
    
    
    # Choose best k
    current_best_accuracy <-
      max(current_tuning_results$accuracy)
    
    current_k <- min(
      current_tuning_results$k[
        current_tuning_results$accuracy ==
          current_best_accuracy
      ]
    )
    
    loocv_selected_k[current_test_item] <- current_k
    
    
    # ----------------------------------------------------------
    # Predict the outer held-out observation
    # ----------------------------------------------------------
    pred <- knn_from_dissimilarity(
      dissimilarity_matrix = dissimilarity_matrix,
      train_items = current_train_items,
      test_items = current_test_item,
      train_labels = current_train_labels,
      k = current_k
    )
    
    loocv_predictions[current_test_item] <- pred
  }
  
  
  # ============================================================
  # 3. NESTED-LOOCV PERFORMANCE
  # ============================================================
  
  confusion_matrix_LOOCV <- table(
    Predicted = loocv_predictions[all_items],
    Actual = all_labels[all_items]
  )
  
  accuracy_LOOCV <- sum(diag(confusion_matrix_LOOCV)) /
    sum(confusion_matrix_LOOCV)
  
  
  # ============================================================
  # 4. Return results
  # ============================================================
  
  return(
    list(
      # Train/test result
      confusion_matrix_random = confusion_matrix_random,
      accuracy_random = accuracy_random,
      
      # LOOCV result
      confusion_matrix_LOOCV = confusion_matrix_LOOCV,
      accuracy_LOOCV = accuracy_LOOCV
    )
  )
}

repeated_stratified_evaluation <- function(
    dissimilarity_matrix,
    all_labels,
    n_repeats = 100,
    train_prop = 0.75,
    seed = 123,
    max_k = NULL
) {

  all_items <- names(all_labels)
  
  if (is.null(all_items)) {
    stop("all_labels must be a named vector.")
  }

  if (!all(all_items %in% rownames(dissimilarity_matrix))) {
    stop("Some items in all_labels are missing from the dissimilarity matrix.")
  }

  class_items <- split(all_items, all_labels)

  class_sizes <- sapply(class_items, length)
  
  if (any(class_sizes < 2)) {
    stop("Each class must contain at least 2 observations.")
  }

  accuracy_vec <- numeric(n_repeats)
  k_vec <- numeric(n_repeats)
  
  confusion_matrices <- vector("list", n_repeats)
  train_sets <- vector("list", n_repeats)
  test_sets <- vector("list", n_repeats)
  predictions_list <- vector("list", n_repeats)

  set.seed(seed)
  

  for (b in seq_len(n_repeats)) {

    train_items_by_class <- lapply(
      class_items,
      function(items) {
        
        n_class <- length(items)
        
        n_train <- floor(train_prop * n_class)
        n_train <- max(1, min(n_train, n_class - 1))
        
        sample(
          items,
          size = n_train,
          replace = FALSE
        )
      }
    )
    
    train_items <- unlist(
      train_items_by_class,
      use.names = FALSE
    )
    
    test_items <- setdiff(all_items, train_items)
    
    train_labels <- all_labels[train_items]
    true_test_labels <- all_labels[test_items]

    train_sets[[b]] <- train_items
    test_sets[[b]] <- test_items

    train_dissimilarity_matrix <- dissimilarity_matrix[
      train_items,
      train_items,
      drop = FALSE
    ]
    
    if (is.null(max_k)) {
      
      current_max_k <- max(
        21,
        round(sqrt(length(train_items)))
      )
      
    } else {
      
      current_max_k <- max_k
    }

    current_max_k <- min(
      current_max_k,
      length(train_items) - 1
    )

    k_values <- seq(
      1,
      current_max_k,
      by = 2
    )
    
    
    tuning_results <- tune_knn_loocv(
      dissimilarity_matrix = train_dissimilarity_matrix,
      all_items = train_items,
      all_labels = train_labels,
      k_values = k_values
    )

    k_choose <- tuning_results$k[
      which.max(tuning_results$accuracy)
    ]
    
    k_vec[b] <- k_choose

    predictions <- knn_from_dissimilarity(
      dissimilarity_matrix = dissimilarity_matrix,
      train_items = train_items,
      test_items = test_items,
      train_labels = train_labels,
      k = k_choose
    )
    
    names(predictions) <- test_items
    
    predictions_list[[b]] <- predictions

    label_levels <- sort(unique(as.character(all_labels)))
    
    confusion_matrix <- table(
      Predicted = factor(
        predictions,
        levels = label_levels
      ),
      Actual = factor(
        true_test_labels,
        levels = label_levels
      )
    )
    
    confusion_matrices[[b]] <- confusion_matrix
    
    
    accuracy_vec[b] <- mean(
      as.character(predictions) ==
        as.character(true_test_labels)
    )
  }

  aggregate_confusion_matrix <- Reduce(
    "+",
    confusion_matrices
  )

  summary_results <- data.frame(
    Mean_Accuracy = mean(accuracy_vec),
    SD_Accuracy = sd(accuracy_vec),
    Median_Accuracy = median(accuracy_vec),
    Min_Accuracy = min(accuracy_vec),
    Max_Accuracy = max(accuracy_vec),
    Median_k = median(k_vec)
  )

  return(
    list(
      summary = summary_results,
      mean_accuracy = mean(accuracy_vec),
      sd_accuracy = sd(accuracy_vec),
      accuracies = accuracy_vec,
      selected_k = k_vec,
      aggregate_confusion_matrix = aggregate_confusion_matrix,
      confusion_matrices = confusion_matrices,
      train_sets = train_sets,
      test_sets = test_sets,
      predictions = predictions_list
    )
  )
}



### simulation
rSkewSym <- function(d,kappa) {
  m2 <- c(1:d) / sqrt(sum((c(1:d))^2))
  m1 <- rmovMF(1, kappa*m2)
  opt_trans <- log_rotation(m2,m1) ## m1-m2
  M <- opt_trans$L
  return(M)
}


AR_generate_function <- function(seed_set, sample_size, d, coef,burnin=500) {
  set.seed(seed_set)
  matrix_mean <- 0 ## should be zero
  matrix_list <- lapply(1:length(coef), function(i) rSkewSym(d, 5*d))
  for (t in (length(coef)+1):(sample_size + burnin + length(coef))) {
    matrix_save <- matrix_mean
    for (coef_num in 1:length(coef)) {
      matrix_save <- matrix_save + coef[coef_num] * (matrix_list[[t-coef_num]] - matrix_mean)
    }
    matrix_list[[t]] <- matrix_save + rSkewSym(d, 20*d)
  }
  return(matrix_list[(length(coef) + 1 + burnin):(sample_size + burnin + length(coef))])
}

# 1. Helper Functions -----------------------------------------------------

# Calculate the Euclidean norm of a vector
norm_vec <- function(x) {
  sqrt(sum(x^2))
}

# Geodesic distance (arc-cosine) between two unit vectors
geo_dist <- function(u, v) {
  val <- sum(u * v)
  val <- min(1, max(-1, val))
  acos(val)
}


get_rotation_matrix <- function(g1, g3) { ## g3-g1
  p <- length(g1)
  theta <- geo_dist(g1, g3)

  if (theta < 1e-10) {
    return(diag(p))
  }

  u1 <- g1
  proj <- g3 - sum(g3 * u1) * u1
  u2 <- proj / norm_vec(proj)
  Q <- (u2 %*% t(u1)) - (u1 %*% t(u2))

  Q2 <- Q %*% Q
  I <- diag(p)
  R <- I + sin(theta) * Q + (1 - cos(theta)) * Q2
  
  return(R)
}

exp_from_log_rotation_predict <- function(q,theta) {
  th <- theta
  L <- q
  if (is.null(th)) {
    eig <- eigen(L, symmetric = FALSE)
    V <- eig$vectors; D <- diag(exp(eig$values))
    Re(V %*% D %*% solve(V))
  } else {
    I <- diag(nrow(L))
    K <- L / max(th, 1e-16)
    I + sin(th) * K + (1 - cos(th)) * (K %*% K)
  }
}


trend_component_function <- function(u, d, inc_idx = c(1,2,3), A = 0.6) {
  mu <- c(1:d) / sqrt(sum((c(1:d)^2)))

  proj_tangent <- function(v) {
    v - sum(v * mu) * mu
  }

  trend_raw <- rep(0, d)
  trend_raw[inc_idx] <- 1
  v_dir <- proj_tangent(trend_raw)
  v_dir <- v_dir / sqrt(sum(v_dir^2))

  v <- A * (u - 0.5) * v_dir
  nv <- sqrt(sum(v^2))
  if (nv < 1e-10) return(mu)
  
  cos(nv) * mu + sin(nv) * v / nv
}

trend_generate_funtion <- function(seed_set, period_matrix, inc_ind, A){
  set.seed(seed_set)
  sample_size <- nrow(period_matrix)
  d <- ncol(period_matrix)
  mu <- c(1:d) / sqrt(sum((c(1:d)^2)))
  
  u_grid <- (1:sample_size) / sample_size
  
  trend_function_matrix <- t(sapply(u_grid, trend_component_function,
                                    d,
                                    inc_ind,
                                    A))
  
  f_mean <- intrinsic_mean_sphere(trend_function_matrix)
  coin_flip <- runif(1)
  if (coin_flip < 0.5) {
    random_a <- runif(1, min = 0, max = 2)
  } else {
    random_a <- runif(1, min = 5, max = 15)
  }
  
  new_matrix <- trend_function_matrix
  for (t in 1:nrow(trend_function_matrix)) {
    random_rataion <- random_a * log_rotation(f_mean,trend_function_matrix[t,])$L
    angle_n <- sqrt(sum((random_rataion %*% f_mean)^2))
    new_matrix[t,] <- exp_from_log_rotation_predict(random_rataion,angle_n) %*% f_mean
  }

  y_matrix <- matrix(NA,nrow = sample_size,ncol = d)
  for (t in 1:sample_size) {
    y_matrix[t,] <- get_rotation_matrix(mu,new_matrix[t,]) %*% period_matrix[t,]
  }
  
  return(y_matrix)
}

periodic_component_function <- function(sample_size,d,period){
  mu <- c(1:d) / sqrt(sum((c(1:d)^2)))
  raw_basis <- pracma::nullspace(t(mu))
  u <- raw_basis[, 1]
  v <- raw_basis[, 2]

  theta <- 0.05
  g_period <- matrix(0, nrow = period, ncol = d)
  for (t in 1:period) {
    phi <- 2 * pi * (t / period)
    g_period[t, ] <- cos(theta) * mu + sin(theta) * (cos(phi) * u + sin(phi) * v)
  }
  return(g_period)
}

period_generate_funtion <- function(seed_set,sphere_matrix){
  set.seed(seed_set)
  sample_size <- nrow(sphere_matrix)
  d <- ncol(sphere_matrix)
  mu <- c(1:d) / sqrt(sum((c(1:d)^2)))

  periodic_function_matrix <- periodic_component_function(sample_size,d,period=sample_size)
  
  
  f_mean <- intrinsic_mean_sphere(periodic_function_matrix)
  random_a <- runif(1,1,10)
  new_matrix <- periodic_function_matrix
  for (t in 1:nrow(periodic_function_matrix)) {
    random_rataion <- random_a * log_rotation(f_mean,periodic_function_matrix[t,])$L
    angle_n <- sqrt(sum((random_rataion %*% f_mean)^2))
    new_matrix[t,] <- exp_from_log_rotation_predict(random_rataion,angle_n) %*% f_mean
  }
  
  
  for (t in 1:sample_size) {
    current_ind <- t + 12 - floor((t + 12 - 1)/12) * 12
    sphere_matrix[t,] <- get_rotation_matrix(mu,new_matrix[current_ind,]) %*% sphere_matrix[t,]
  }
  
  return(sphere_matrix)
}




shaple_cluster_sim <- function(seed,function_size,sample.size=100,d_set=6,true_AR_coef=c(0.5)){
  num_groups <- 3

  group_assignments <- rep(1:num_groups, length.out = function_size)
  cluster_list <- list()
  
  for (i in 1:function_size) {
    if (group_assignments[i]==1){
      AR_list <- AR_generate_function(seed_set=(i+seed),sample_size=sample.size,d=d_set,coef=true_AR_coef)
      
      sphere_matrix <- matrix(NA,nrow = sample.size,ncol = d_set)
      for (t in 1:sample.size) {
        mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
        angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
        sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
      }
      
      cluster_list[[i]] <- trend_generate_funtion(seed_set=(i+seed), sphere_matrix, inc_ind=c(1,2,3), A=0.3)
      
    } else if (group_assignments[i]==2){
      AR_list <- AR_generate_function(seed_set=(i+seed),sample_size=sample.size,d=d_set,coef=true_AR_coef)
      
      sphere_matrix <- matrix(NA,nrow = sample.size,ncol = d_set)
      for (t in 1:sample.size) {
        mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
        angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
        sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
      }
      
      cluster_list[[i]] <- trend_generate_funtion(seed_set=(i+seed), sphere_matrix, inc_ind=c(1,2,6), A=-0.3)
      
    } else{
      AR_list <- AR_generate_function(seed_set=(i+seed),sample_size=sample.size,d=d_set,coef=true_AR_coef)
      
      sphere_matrix <- matrix(NA,nrow = sample.size,ncol = d_set)
      for (t in 1:sample.size) {
        mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
        angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
        sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
      }
      
      cluster_list[[i]] <- trend_generate_funtion(seed_set=(i+seed), sphere_matrix, inc_ind=c(1,2,5), A=-0.1)
      
    } 
  }
  
  names(cluster_list) <- group_assignments
  return(cluster_list)
}


