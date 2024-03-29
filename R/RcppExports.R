# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rescaled <- function(beta_temp, xi_temp, Ustar_dems) {
    .Call(`_EPAbvs_rescaled`, beta_temp, xi_temp, Ustar_dems)
}

bvsEPAcpp <- function(iterations, n_clust_betn, n_shuffle, between_covas, cluster_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, dist, subj_clust_free, beta_init, xi_init, s_init, gamma_init, w, alpha, delta, gamma_delta, eta, sigma, mu, lambda2_init, tau2_init, nu_lambda_init, nu_tau_init, a_alpha, b_alpha, a_delta, b_delta, a_gamma, b_gamma, a_eta, b_eta, ftype) {
    .Call(`_EPAbvs_bvsEPAcpp`, iterations, n_clust_betn, n_shuffle, between_covas, cluster_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, dist, subj_clust_free, beta_init, xi_init, s_init, gamma_init, w, alpha, delta, gamma_delta, eta, sigma, mu, lambda2_init, tau2_init, nu_lambda_init, nu_tau_init, a_alpha, b_alpha, a_delta, b_delta, a_gamma, b_gamma, a_eta, b_eta, ftype)
}

bvsPYcpp <- function(iterations, n_clust_betn, n_shuffle, between_covas, cluster_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, subj_clust_free, beta_init, xi_init, s_init, gamma_init, w, alpha, delta, mu, lambda2_init, tau2_init, nu_lambda_init, nu_tau_init, a_gamma, b_gamma, a_alpha, b_alpha, a_delta, b_delta, DP) {
    .Call(`_EPAbvs_bvsPYcpp`, iterations, n_clust_betn, n_shuffle, between_covas, cluster_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, subj_clust_free, beta_init, xi_init, s_init, gamma_init, w, alpha, delta, mu, lambda2_init, tau2_init, nu_lambda_init, nu_tau_init, a_gamma, b_gamma, a_alpha, b_alpha, a_delta, b_delta, DP)
}

