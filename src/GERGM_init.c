#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GERGM_Extended_Metropolis_Hastings_Sampler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GERGM_extended_weighted_mple_objective(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GERGM_frobenius_norm(SEXP, SEXP);
extern SEXP GERGM_get_indiviual_triad_values(SEXP, SEXP, SEXP, SEXP);
extern SEXP GERGM_get_triad_weights(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GERGM_h_statistics(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GERGM_Individual_Edge_Conditional_Prediction(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GERGM_Metropolis_Hastings_Sampler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GERGM_weighted_mple_objective(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GERGM_Extended_Metropolis_Hastings_Sampler",   (DL_FUNC) &GERGM_Extended_Metropolis_Hastings_Sampler,   28},
    {"GERGM_extended_weighted_mple_objective",       (DL_FUNC) &GERGM_extended_weighted_mple_objective,       16},
    {"GERGM_frobenius_norm",                         (DL_FUNC) &GERGM_frobenius_norm,                          2},
    {"GERGM_get_indiviual_triad_values",             (DL_FUNC) &GERGM_get_indiviual_triad_values,              4},
    {"GERGM_get_triad_weights",                      (DL_FUNC) &GERGM_get_triad_weights,                       5},
    {"GERGM_h_statistics",                           (DL_FUNC) &GERGM_h_statistics,                           12},
    {"GERGM_Individual_Edge_Conditional_Prediction", (DL_FUNC) &GERGM_Individual_Edge_Conditional_Prediction, 30},
    {"GERGM_Metropolis_Hastings_Sampler",            (DL_FUNC) &GERGM_Metropolis_Hastings_Sampler,            16},
    {"GERGM_weighted_mple_objective",                (DL_FUNC) &GERGM_weighted_mple_objective,                10},
    {NULL, NULL, 0}
};

void R_init_GERGM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}