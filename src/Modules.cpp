#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4expmodn_mod) {


    class_<rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> >("model_expmodn")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_expmodn_namespace::model_expmodn, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4inv_gaussian_mod) {


    class_<rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> >("model_inv_gaussian")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_inv_gaussian_namespace::model_inv_gaussian, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4logdexp_mod) {


    class_<rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> >("model_logdexp")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_logdexp_namespace::model_logdexp, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4loggumbel_mod) {


    class_<rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> >("model_loggumbel")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_loggumbel_namespace::model_loggumbel, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4logistic_mod) {


    class_<rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> >("model_logistic")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_logistic_namespace::model_logistic, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4lognormal_mod) {


    class_<rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> >("model_lognormal")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_lognormal_namespace::model_lognormal, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4weib_mod) {


    class_<rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> >("model_weib")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_weib_namespace::model_weib, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
