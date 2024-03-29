# HELPER FUNCTION FOR plot_statistic_PCA.R
# This is a simplified and modified version of the 'ggbiplot' function from
# the 'ggbiplot' package. The original function does not allow customization of
# a couple aesthetic features that need to be changed for our figures.

custom_ggbiplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
          obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
          ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
          alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
          varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE,
          arrow.size = 0.75, var.color = "darkred", overlap_fix = FALSE, ell.size = 0.75,
          point.size = 1, plot.outside = TRUE,
          ...) 
{
  library(tidyverse)
  #library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  
  nobs.factor <- sqrt(nrow(pcobj$x) - 1)
  d <- pcobj$sdev
  u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
  v <- pcobj$rotation
  
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  } else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  df.v2 <- df.v
  # if (overlap_fix) {
  #   df.v2$angle[1] <- df.v2$angle[1] + 3.2
  #   df.v2$angle[4] <- df.v2$angle[4] - 0.5
  #   df.v2$length <- sqrt(df.v2$xvar**2 + df.v2$yvar**2)
  #   df.v2$length[1] <- df.v2$length[1] * 0.945
  #   df.v2$xvar[c(1,4)] <- df.v2$length[c(1,4)] * cos(df.v2$angle[c(1,4)] * (pi/180))
  #   df.v2$yvar[c(1,4)] <- df.v2$length[c(1,4)] * sin(df.v2$angle[c(1,4)] * (pi/180))
  # }
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar),
                          arrow = arrow(length = unit(arrow.size,"picas")),
                          color = var.color, lwd=arrow.size)
  }
  
  ellipse_params <- list()
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 100), seq(pi, -pi, length = 100))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- data.frame(X1=0,X2=0,groups='')[0,]
    ell.params <- list()
    for (group in sort(unique(df.u$groups))) {
      x <- df.u %>% filter(groups==group)
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      
      ell <- rbind(ell,
                   data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = "+"), groups = x$groups[1]))
      ell.params[[group]] <- list(sigma=sigma,mu=mu,ed=ed)
      
    }
    
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups), size=ell.size)
  }
  
  # determines if a point is inside its group' ellipse
  mahalanobis_distance <- function(x, center, cov_matrix) {
    diff <- matrix(x - center, ncol = 2)
    sqrt((diff %*% solve(cov_matrix) %*% t(diff)))
  }
  is_point_out_ellipse <- function(x, y, group, ell.params) {
    params <- ell.params[[group]]
    md <- mahalanobis_distance(c(x, y), params$mu, params$sigma)[[1]]
    md > params$ed
  }
  
  df.u <- df.u %>%
    rowwise() %>%
    mutate(outside.ell = is_point_out_ellipse(xvar, yvar, groups, ell.params))
  
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      if (plot.outside) {
        df.u.out <- df.u
        df.u.in <- df.u[0,]
      } else {
        df.u.out <- df.u %>% filter(outside.ell)
        df.u.in <- df.u %>% filter(!outside.ell)
      }
      
      g <- g + geom_text(data=df.u.out,
                         aes(label = labels, color = groups), 
                         size = labels.size) +
        geom_point(data = df.u.in,
                   aes(color = groups), alpha = alpha, size=point.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  } else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha, size=point.size)
    }
    else {
      g <- g + geom_point(alpha = alpha, size=point.size)
    }
  }
  
  
  if (var.axes) {
    g <- g + geom_text(data = df.v2, aes(label = varname, 
                                        x = xvar, y = yvar, angle = angle, hjust = hjust), 
                       color = var.color, size = varname.size)
  }
  return(g)
}
