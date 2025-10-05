set.seed(111)
library(ggplot2)
# We start by defining the functions we will need later

# the regression function m
reg_m <- function(x){
  sin(1/(x/3 + 0.1))
}

block_assignment <- function(x, N) {
  cut(x, breaks = seq(0, 1, length.out = N + 1),
      include.lowest = T, right = T, labels = F)
}


# Second derivative of an order 4 polynomial used to approxiate theta_22 and sigma2
poly_2nd <- function(x, coeff) {
  coeff <- c(coeff, rep(0, max(0, 5 - length(coeff))))
  2*coeff[3] + 6*coeff[4]*x + 12*coeff[5]*x^2
}

# Function that will fit a polynomial linear model of order 4 in each block j=1, ..., N if there are
# enough data points in the block
fit_blocks <- function(x, y, N){
  liste_coeff <- vector("list", N)
  blocks_rss <- numeric(N)
  m2nd_vec <- numeric(length(x))
  
  z <- block_assignment(x, N)
  for (j in seq_len(N)){
    index <- which(z == j)
    xj <- x[index]
    yj <- y[index]
    if (length(index) < 6) {                           
      return(list(valid = F))                        
    }
    
    Xj <- cbind(1, xj, xj^2, xj^3, xj^4)
    fitj <- lm.fit(Xj, yj)
    coeffj <- coef(fitj)
    if (any(!is.finite(coeffj))) return(list(valid = F))
    
    liste_coeff[[j]] <- coeffj
    fitted_j <- drop(Xj %*% coeffj)
    resid_j  <- yj - fitted_j
    blocks_rss[j] <- sum(resid_j^2)
  }
  
   for (j in seq_len(N)){
    index <- which(z == j)
    m2nd_vec[index] <- poly_2nd(x[index], liste_coeff[[j]])
  }
  
  list(
    valid = T,
    z=z,
    betas_hat = liste_coeff,
    m2nd = m2nd_vec,
    total_rss = sum(blocks_rss)
  )
}


# Estimate theta22_hat(N) and sigma2_hat(N) 
estimates <- function(x, y, N) {
  n <- length(x)
  fit <- fit_blocks(x, y, N)
  
  if (!isTRUE(fit$valid)) return(list(valid = F))
  theta22_hat <- mean( fit$m2nd^2 )

  if ((n - 5 * N) <= 0) return(list(valid = F))
  sigma2_hat <- fit$total_rss / (n - 5 * N)
  
  list(
    valid = T,
    theta22_hat = theta22_hat,
    sigma2_hat = sigma2_hat,
    rss = fit$total_rss
  )
}

# compute the approximation of the optimal bandwidth
h_AMISE_hat <- function(n, sigma2_hat, theta22_hat, support_len = 1) {
  if (!is.finite(theta22_hat) || theta22_hat <= 0) return(NA_real_)
  n^(-1/5) * ((35 * sigma2_hat * support_len) / theta22_hat)^(1/5)
}

# Mallows' Cp function
Cp_N <- function(n, N, RSS_N, RSS_Nmax, Nmax) {
  (RSS_N / ( RSS_Nmax / (n - 5 * Nmax) )) - (n - 10 * N)
}

# compute Nmax 
Nmax_f <- function(n) {
  max( min(floor(n / 20), 5), 1 )
}



# choose the optimal value of N with respect to Mallow's Cp criterion
N_opt_Cp <- function(x, y) {
  n    <- length(x)
  Nmax <- Nmax_f(n)

  RSS_Nmax <- NA_real_
  for (Nm in Nmax:1) {
    est_try <- estimates(x, y, Nm)
    if (isTRUE(est_try$valid)) { Nmax <- Nm; RSS_Nmax <- est_try$rss; break }
  }
  if (!is.finite(RSS_Nmax)) return(NA_integer_) 
  
  best_N  <- NA_integer_
  best_Cp <- Inf
  
  for (N in 1:Nmax) {
    estN <- estimates(x, y, N)
    if (!isTRUE(estN$valid)) next
    Cp <- Cp_N(n, N, estN$rss, RSS_Nmax, Nmax)
    if (Cp < best_Cp || (is.finite(Cp) && abs(Cp - best_Cp) <= 1e-12 && N < best_N)) {
      best_Cp <- Cp
      best_N  <- N
    }
  }
  best_N
}


# vary sample size n and fix parameters alpha, beta and sigma
alpha <- 2
beta <- 2
sigma <- 1


#  ĥ for a given using the Cp-optimalN
h_hat_for_n <- function(n, alpha, beta, sigma) {
  x <- rbeta(n, alpha, beta)
  y <- reg_m(x) + rnorm(n, sd = sqrt(sigma))
  N_opt <- N_opt_Cp(x, y)
  if (is.na(N_opt)) {
    data.frame(n = n, h_hat = NA_real_, N_opt = NA_integer_)
  } else {
    est <- estimates(x, y, N_opt)
    h_hat <- h_AMISE_hat(n, est$sigma2_hat, est$theta22_hat, support_len = 1)
    data.frame(n = n, h_hat = h_hat, N_opt = N_opt)
  }
}

# values of n to try
n_vals <- c(100, 200, 400, 800, 1600)

# run and collect
res_h <- do.call(rbind, lapply(n_vals, h_hat_for_n,
                               alpha = alpha, beta = beta, sigma = sigma))

# plot 
use_repel <- requireNamespace("ggrepel", quietly = TRUE)

plt1 <- ggplot(res_h, aes(n, h_hat)) +
  geom_line(linewidth = 1.1, alpha = 0.9, na.rm = TRUE) +
  geom_point(size = 3.6, shape = 21, stroke = 1.1, fill = "white", na.rm = TRUE) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = "Sample size n",
    y = expression(hat(h)[AMISE]),
    title = expression("Computed " * hat(h)[AMISE] * " for different sample sizes " * n),
    subtitle = as.expression(bquote(
      "Using Cp-optimal block size " * N * " per dataset; " ~
        alpha == .(alpha) ~ ", " ~ beta == .(beta) ~ ", " ~ sigma^2 == .(sigma)
    ))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

plt1
ggsave("Single_Realization.png", plt1, width = 9, height = 5, dpi = 300)


df_rate <- subset(res_h, is.finite(h_hat))
mod <- lm(log(h_hat) ~ log(n), data = df_rate)
slope <- coef(mod)[2]

fit <- lm(log(h_hat) ~ log(n), data = df_rate)
a <- unname(coef(fit)[1])   # intercept on log scale
b <- unname(coef(fit)[2])   # slope (should be ~ -0.2)

# Proper log–log plot
plt2 <- ggplot(df_rate, aes(x = n, y = h_hat)) +
  geom_point(size = 2.8) +
  # draw the power-law fit h_hat = exp(a) * n^b
  geom_function(fun = function(x) exp(a) * x^b, linewidth = 1, color="blue") +
  scale_x_log10(breaks = sort(unique(df_rate$n))) +
  scale_y_log10() +
  annotate(
    "label",
    x = min(df_rate$n) * 1.05,
    y = max(df_rate$h_hat),
    label = sprintf("slope ≈ %.3f (theory: -0.2)", b),
    hjust = 0, vjust = 1, size = 3.5
  ) +
  labs(
    x = "n (log scale)",
    y = expression(hat(h)[AMISE]~"(log scale)"),
    title = "Log–log scale"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

plt2

ggsave("log_scale_plot.png", plt2, width = 9, height = 5, dpi = 300)


#Replicate the simulation 400 times and plot boxplots for each values of n
R <- 400
res_rep <- do.call(
  rbind,
  lapply(n_vals, function(nn) {
    do.call(rbind, replicate(R, h_hat_for_n(nn, alpha, beta, sigma), simplify = FALSE))
  })
)

plt3 <- ggplot(subset(res_rep, is.finite(h_hat)), aes(factor(n), h_hat)) +
  geom_boxplot(outlier.alpha = 0.25) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.6, fill = "red") +
  labs(x = "Sample size n", y = expression(hat(h)[AMISE]),
       title = expression("Observed "~hat(h)[AMISE]~" distribution across 400 replications")) +
  theme_minimal(base_size = 13)+
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )
ggsave("400_Replications.png", plt3, width = 9, height = 5, dpi = 300)

N_summary <- aggregate(N_opt ~ n, data = subset(res_rep, !is.na(N_opt)), FUN = median)

plt3_bis <-ggplot(subset(res_rep, !is.na(N_opt)), aes(x = factor(N_opt))) +
  geom_bar() +
  facet_wrap(~ n, nrow = 1) +
  labs(
    title = expression("Distribution of " * N[opt] ~ "(Cp) across 400 replications"),
    x = expression(N[opt]),
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
plt3_bis
ggsave("400_rep_N_vs_n.png", plt3_bis, width = 9, height = 5, dpi = 300)

df_mean <- aggregate(h_hat ~ n, data = subset(res_rep, is.finite(h_hat)), FUN = mean)
names(df_mean)[2] <- "mean_h"

# fit slope on log–log means and makes the plot of the means on a log-log scale
mod   <- lm(log(mean_h) ~ log(n), data = df_mean)
slope <- unname(coef(mod)[2])
plt4 <- ggplot(df_mean, aes(x = n, y = mean_h)) +
  geom_point(size = 2.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(breaks = n_vals) +
  scale_y_log10() +
  annotate(
    "label",
    x = min(df_mean$n) * 1.05,
    y = max(df_mean$mean_h),
    label = sprintf("slope ≈ %.3f (theory: -0.2)", slope),
    hjust = 0, vjust = 1, size = 3.5
  ) +
  labs(
    x = "log(n)",
    y = expression(log("mean " * hat(h)[AMISE])),
    title = "Log–log scaling using the mean "~hat(h)[AMISE]~"E over 400 replications",
    subtitle = as.expression(bquote(
      alpha == .(alpha) ~ ", " ~ beta == .(beta) ~ ", " ~ sigma^2 == .(sigma)
    ))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
plt4
ggsave("400_Replications_log_log.png", plt4, width = 9, height = 5, dpi = 300)



# build one table per n
h_vs_N <- function(n, alpha, beta, sigma) {
  x <- rbeta(n, alpha, beta)
  y <- reg_m(x) + rnorm(n, sd = sqrt(sigma))
  Nmax <- Nmax_f(n)
  
  # back off Nmax until feasible for Cp denominator
  RSS_Nmax <- NA_real_
  for (Nm in Nmax:1) {
    est <- estimates(x, y, Nm)
    if (isTRUE(est$valid)) { Nmax <- Nm; RSS_Nmax <- est$rss; break }
  }
  if (!is.finite(RSS_Nmax)) return(NULL)
  
  rows <- lapply(1:Nmax, function(N){
    estN <- estimates(x, y, N)
    if (!isTRUE(estN$valid)) return(NULL)
    data.frame(
      n = n,
      N = N,
      h_hat = h_AMISE_hat(n, estN$sigma2_hat, estN$theta22_hat),
      Cp = Cp_N(n, N, estN$rss, RSS_Nmax, Nmax)
    )
  })
  do.call(rbind, rows)
}

n_vals <- c(200, 800, 3200, 12800)
df <- do.call(rbind, lapply(n_vals, h_vs_N, alpha = alpha, beta = beta, sigma = sigma))

# optimal N per n (Cp criterion), for the dashed vertical lines
vlines <- do.call(rbind, lapply(split(df, df$n), function(d)
  data.frame(n = unique(d$n), N_opt = d$N[which.min(d$Cp)])
))

# fix a common y-scale across panels
y_rng <- range(df$h_hat[is.finite(df$h_hat)], na.rm = TRUE)

plt5 <- ggplot(df, aes(N, h_hat)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.4, shape = 21, stroke = 1.2, fill = "red", alpha = 0.9) +
  geom_vline(data = vlines_clean, aes(xintercept = N_opt), linetype = 2, linewidth = 0.6) +
  facet_wrap("n="~ n) +                     
  scale_x_continuous(breaks = seq_len(max(df$N))) +
  scale_y_continuous(limits = y_rng, expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Estimated optimal bandwidth vs block size N",
    x = "Block size (N)",
    y = expression(hat(h)[AMISE])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
plt5
ggsave("bandwidth_vs_N.png", plt5, width = 9, height = 5, dpi = 300)


# We now fix n=800, sigma=1, the optimal N, (w.r.t. Mallow's Cp criterion) and vary alpha and beta
n     <- 800
shapes <- list(c(0.5,0.5), c(5,1), c(1,3), c(2,2), c(1,5), c(5,5))

#Show the different density function for the different beta parametes pairs
# (alpha, beta) pairs
params <- data.frame(
  alpha = c(0.5, 5,   1, 2, 1, 5),
  beta  = c(0.5, 1,   3, 2, 5, 5)
)

x <- seq(1e-6, 1 - 1e-6, length.out = 1000)
labels <- with(params, paste0("Beta(", alpha, ",", beta, ")"))

df <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
  data.frame(
    x = x,
    density = dbeta(x, params$alpha[i], params$beta[i]),
    label = labels[i]
  )
}))
df$label <- factor(df$label, levels = labels)

# y-axis capped at 5
plt6 <- ggplot(df, aes(x, density, color = label)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Beta distributions",
    x = "x",
    y = "Density",
    color = "Parameters"
  ) +
  coord_cartesian(ylim = c(0, 5)) + 
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

plt6

ggsave("pdf_beta_shapes.png", plt6, width = 9, height = 5, dpi = 300)


# Fix for some rounding issues i encountered
params <- data.frame(
  alpha = vapply(shapes, `[[`, numeric(1), 1),
  beta  = vapply(shapes, `[[`, numeric(1), 2)
)
params$label <- sprintf("Beta(%s,%s)",
                        format(params$alpha, trim = TRUE, scientific = FALSE),
                        format(params$beta,  trim = TRUE, scientific = FALSE))

# single replicate
h_one <- function(alpha, beta, label) {
  x <- rbeta(n, alpha, beta)
  y <- reg_m(x) + rnorm(n, sd = sqrt(sigma))
  Nopt <- N_opt_Cp(x, y)
  if (is.na(Nopt)) return(data.frame(label = label, h_hat = NA_real_))
  est <- estimates(x, y, Nopt)
  if (!isTRUE(est$valid)) return(data.frame(label = label, h_hat = NA_real_))
  data.frame(label = label,
             h_hat = h_AMISE_hat(n, est$sigma2_hat, est$theta22_hat))
}


res_ab <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
  do.call(rbind, replicate(
    R,
    h_one(params$alpha[i], params$beta[i], params$label[i]),
    simplify = FALSE
  ))
}))

res_ab$label <- factor(res_ab$label, levels = params$label)


# Plot
plt7 <- ggplot(subset(res_ab, is.finite(h_hat)), aes(x = label, y = h_hat)) +
  geom_boxplot() +
  labs(
    title = "Estimated optimal bandwidth across Beta shapes (n = 800)",
    subtitle = sprintf("R = %d replications per shape", R),
    x = "Design density",
    y = expression(hat(h)[AMISE])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

plt7
ggsave("Boxplots_beta_shapes.png", plt7, width = 9, height = 5, dpi = 300)

