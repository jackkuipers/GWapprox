
library(igraph)
library(GWnorm)

for (p in 2*3:6){ # size of graph
n_samp_GW <- 1e6

gamma_vec <- 6*(3/2)^(0:6)

n_loop <- 100

out_data <- NULL

n_seeds <- 25

if (!file.exists(paste0("./sim_results/g_approx", p, ".rData"))) {

for (seed_number in 1:n_seeds){

set.seed(100+seed_number)
d <- 1.5 # density of 2 gives good chance of connected prime graph
edge_prob <- d*2/(p-1) # for reasonable p at least?!

counter <- 0
good_graph <- FALSE
while (!good_graph) {
  graph <- sample_gnp(p, edge_prob)
  good_graph <- check_prime_connected(graph)
  counter <- counter + 1
}

beta <- 1
adj <- as.matrix(as_adjacency_matrix(graph, type = "upper"))
prec_mat <- BDgraph::rgwish(adj = adj)
sigma <- cov2cor(solve(prec_mat))
data <- BDgraph::rmvnorm(n = 2*p, sigma = sigma) # need at least the dimension plus extras
D <- cor(data) # keep it in correlation form
# doesn't depend on scale of data at least

D_tilde <- PD_complete(graph, D)

tau <- I_Gnorm(graph, beta, D, err_flag = TRUE)$tau

for (ii in 1:length(gamma_vec)) {
  
gamma <- gamma_vec[ii]
beta <- gamma - (p+1)/2

start.time <- Sys.time()

# for D
GW_outs <- rep(NA, n_loop)
for(k in 1:n_loop){
  GW_outs[k] <- I_G_MC(graph, beta, D, n_samp = n_samp_GW)
}
max_GW <- max(GW_outs, na.rm = TRUE)
sd_estD <- sd(exp(GW_outs-max_GW), na.rm = TRUE)
IG_avD <- log(mean(exp(GW_outs-max_GW), na.rm = TRUE)) + max_GW

# for I
GW_outs <- rep(NA, n_loop)
for(k in 1:n_loop){
  GW_outs[k] <- I_G_MC(graph, beta, diag(p), n_samp = n_samp_GW)
}
max_GW <- max(GW_outs, na.rm = TRUE)
sd_estI <- sd(exp(GW_outs-max_GW), na.rm = TRUE)
IG_avI <-log(mean(exp(GW_outs-max_GW), na.rm = TRUE)) + max_GW

serr_est <- sqrt((sd_estD^2 + sd_estI^2)/n_loop)

IG_rat <- IG_avD - IG_avI
IG_rat_app <- I_G_ratio_approx(graph, beta, D)
IG_rat_app_nlo <- I_G_ratio_approx(graph, beta, D, use_nlo = TRUE)

rel_err <- IG_rat_app - IG_rat
rel_err_nlo <- IG_rat_app_nlo - IG_rat 

end.time <- Sys.time()
time.taken <- (end.time - start.time)/n_loop
out_data_local <- data.frame(seed = seed_number, p = p, tau = tau, gamma = gamma, time = time.taken, 
                             IG_rat = IG_rat, IG_rat_app = IG_rat_app, IG_rat_app_nlo = IG_rat_app_nlo, 
                             rel_err = rel_err, rel_err_nlo = rel_err_nlo,
                             serr_est = serr_est)
print(out_data_local)
out_data <- rbind(out_data, out_data_local)

}

}

save(out_data, file = paste0("./sim_results/g_approx", p, ".rData"))

}
}
