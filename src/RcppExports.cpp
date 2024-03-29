// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rescaled
List rescaled(arma::mat beta_temp, arma::mat xi_temp, arma::vec Ustar_dems);
RcppExport SEXP _EPAbvs_rescaled(SEXP beta_tempSEXP, SEXP xi_tempSEXP, SEXP Ustar_demsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type beta_temp(beta_tempSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xi_temp(xi_tempSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ustar_dems(Ustar_demsSEXP);
    rcpp_result_gen = Rcpp::wrap(rescaled(beta_temp, xi_temp, Ustar_dems));
    return rcpp_result_gen;
END_RCPP
}
// bvsEPAcpp
List bvsEPAcpp(int iterations, int n_clust_betn, int n_shuffle, arma::vec between_covas, arma::vec cluster_covas, arma::vec subject, arma::vec subject_dems, arma::vec Y, arma::mat Xbar, arma::mat Ustar, arma::vec Ustar_dems, arma::mat dist, arma::vec subj_clust_free, arma::mat beta_init, arma::mat xi_init, arma::mat s_init, arma::mat gamma_init, arma::vec w, arma::vec alpha, arma::vec delta, arma::vec gamma_delta, arma::vec eta, arma::uvec sigma, arma::vec mu, arma::mat lambda2_init, arma::vec tau2_init, arma::mat nu_lambda_init, arma::vec nu_tau_init, double a_alpha, double b_alpha, double a_delta, double b_delta, double a_gamma, double b_gamma, double a_eta, double b_eta, int ftype);
RcppExport SEXP _EPAbvs_bvsEPAcpp(SEXP iterationsSEXP, SEXP n_clust_betnSEXP, SEXP n_shuffleSEXP, SEXP between_covasSEXP, SEXP cluster_covasSEXP, SEXP subjectSEXP, SEXP subject_demsSEXP, SEXP YSEXP, SEXP XbarSEXP, SEXP UstarSEXP, SEXP Ustar_demsSEXP, SEXP distSEXP, SEXP subj_clust_freeSEXP, SEXP beta_initSEXP, SEXP xi_initSEXP, SEXP s_initSEXP, SEXP gamma_initSEXP, SEXP wSEXP, SEXP alphaSEXP, SEXP deltaSEXP, SEXP gamma_deltaSEXP, SEXP etaSEXP, SEXP sigmaSEXP, SEXP muSEXP, SEXP lambda2_initSEXP, SEXP tau2_initSEXP, SEXP nu_lambda_initSEXP, SEXP nu_tau_initSEXP, SEXP a_alphaSEXP, SEXP b_alphaSEXP, SEXP a_deltaSEXP, SEXP b_deltaSEXP, SEXP a_gammaSEXP, SEXP b_gammaSEXP, SEXP a_etaSEXP, SEXP b_etaSEXP, SEXP ftypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type n_clust_betn(n_clust_betnSEXP);
    Rcpp::traits::input_parameter< int >::type n_shuffle(n_shuffleSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type between_covas(between_covasSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cluster_covas(cluster_covasSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type subject(subjectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type subject_dems(subject_demsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xbar(XbarSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ustar(UstarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ustar_dems(Ustar_demsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist(distSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type subj_clust_free(subj_clust_freeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xi_init(xi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type s_init(s_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma_init(gamma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_delta(gamma_deltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda2_init(lambda2_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau2_init(tau2_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type nu_lambda_init(nu_lambda_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu_tau_init(nu_tau_initSEXP);
    Rcpp::traits::input_parameter< double >::type a_alpha(a_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type b_alpha(b_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type a_delta(a_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type b_delta(b_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type a_gamma(a_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type b_gamma(b_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a_eta(a_etaSEXP);
    Rcpp::traits::input_parameter< double >::type b_eta(b_etaSEXP);
    Rcpp::traits::input_parameter< int >::type ftype(ftypeSEXP);
    rcpp_result_gen = Rcpp::wrap(bvsEPAcpp(iterations, n_clust_betn, n_shuffle, between_covas, cluster_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, dist, subj_clust_free, beta_init, xi_init, s_init, gamma_init, w, alpha, delta, gamma_delta, eta, sigma, mu, lambda2_init, tau2_init, nu_lambda_init, nu_tau_init, a_alpha, b_alpha, a_delta, b_delta, a_gamma, b_gamma, a_eta, b_eta, ftype));
    return rcpp_result_gen;
END_RCPP
}
// bvsPYcpp
List bvsPYcpp(int iterations, int n_clust_betn, int n_shuffle, arma::vec between_covas, arma::vec cluster_covas, arma::vec subject, arma::vec subject_dems, arma::vec Y, arma::mat Xbar, arma::mat Ustar, arma::vec Ustar_dems, arma::vec subj_clust_free, arma::mat beta_init, arma::mat xi_init, arma::mat s_init, arma::mat gamma_init, arma::vec w, arma::vec alpha, arma::vec delta, arma::vec mu, arma::mat lambda2_init, arma::vec tau2_init, arma::mat nu_lambda_init, arma::vec nu_tau_init, double a_gamma, double b_gamma, double a_alpha, double b_alpha, double a_delta, double b_delta, bool DP);
RcppExport SEXP _EPAbvs_bvsPYcpp(SEXP iterationsSEXP, SEXP n_clust_betnSEXP, SEXP n_shuffleSEXP, SEXP between_covasSEXP, SEXP cluster_covasSEXP, SEXP subjectSEXP, SEXP subject_demsSEXP, SEXP YSEXP, SEXP XbarSEXP, SEXP UstarSEXP, SEXP Ustar_demsSEXP, SEXP subj_clust_freeSEXP, SEXP beta_initSEXP, SEXP xi_initSEXP, SEXP s_initSEXP, SEXP gamma_initSEXP, SEXP wSEXP, SEXP alphaSEXP, SEXP deltaSEXP, SEXP muSEXP, SEXP lambda2_initSEXP, SEXP tau2_initSEXP, SEXP nu_lambda_initSEXP, SEXP nu_tau_initSEXP, SEXP a_gammaSEXP, SEXP b_gammaSEXP, SEXP a_alphaSEXP, SEXP b_alphaSEXP, SEXP a_deltaSEXP, SEXP b_deltaSEXP, SEXP DPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type n_clust_betn(n_clust_betnSEXP);
    Rcpp::traits::input_parameter< int >::type n_shuffle(n_shuffleSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type between_covas(between_covasSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cluster_covas(cluster_covasSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type subject(subjectSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type subject_dems(subject_demsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xbar(XbarSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ustar(UstarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ustar_dems(Ustar_demsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type subj_clust_free(subj_clust_freeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xi_init(xi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type s_init(s_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma_init(gamma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda2_init(lambda2_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau2_init(tau2_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type nu_lambda_init(nu_lambda_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu_tau_init(nu_tau_initSEXP);
    Rcpp::traits::input_parameter< double >::type a_gamma(a_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type b_gamma(b_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a_alpha(a_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type b_alpha(b_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type a_delta(a_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type b_delta(b_deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type DP(DPSEXP);
    rcpp_result_gen = Rcpp::wrap(bvsPYcpp(iterations, n_clust_betn, n_shuffle, between_covas, cluster_covas, subject, subject_dems, Y, Xbar, Ustar, Ustar_dems, subj_clust_free, beta_init, xi_init, s_init, gamma_init, w, alpha, delta, mu, lambda2_init, tau2_init, nu_lambda_init, nu_tau_init, a_gamma, b_gamma, a_alpha, b_alpha, a_delta, b_delta, DP));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EPAbvs_rescaled", (DL_FUNC) &_EPAbvs_rescaled, 3},
    {"_EPAbvs_bvsEPAcpp", (DL_FUNC) &_EPAbvs_bvsEPAcpp, 37},
    {"_EPAbvs_bvsPYcpp", (DL_FUNC) &_EPAbvs_bvsPYcpp, 31},
    {NULL, NULL, 0}
};

RcppExport void R_init_EPAbvs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
