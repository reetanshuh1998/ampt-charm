// dmeson src_v14: Fortran bindings for the official EPS09s interface
// (Helenius, Eskola, Honkanen, Salgado, arXiv:1205.5359). eps09s.cpp is
// included unmodified; these wrappers replace its commented extern "C" block.
#include "eps09s.cpp"
extern "C" {
void eps09s_c_(int *eps_order, int *eps_pset, int *a, double *x, double *q,
               double *cuv, double *cdv, double *cu, double *cd,
               double *cs, double *cc, double *cb, double *cg){
  eps09s_c(*eps_order, *eps_pset, *a, *x, *q, cuv, cdv, cu, cd, cs, cc, cb, cg);
}
double ta_eps09s_(int *a, double *s){ return ta_eps09s(*a, *s); }
void test_eps09s_(){ test_eps09s(); }
}
