// eps09s_rg: offline EPS09s NLO gluon modification for the error-set systematic (README Sec. 2.3).
// stdin: "x Q2 s" per gluon leg ('# hqshw' x, VINT(52), transverse distance s in fm); stdout: r_g for sets 1..31
// at Q = sqrt(max(hqmuf*Q2, 1.69)) and T_A(s), A = 197, exactly as HQNMOD/HQSHWT (src_new/charm_mod.f) evaluate it.
// Reads eps09sNLO1..31.dat from the working directory (eps09s_errw.py runs it in deps/eps09s/).
// Build (make eps-tool): g++ -O2 -I deps/eps09s analysis/eps09s/eps09s_rg.cpp -o bin/eps09s_rg
#include "eps09s.cpp"
#include <cstdio>
#include <cmath>
int main(int argc, char **argv){
  double muf = argc > 1 ? atof(argv[1]) : 2.0;
  double x, q2, s, c[8][4];
  double cuv[4], cdv[4], cu[4], cd[4], cs[4], cc[4], cb[4], cg[4];
  while (scanf("%lf %lf %lf", &x, &q2, &s) == 3) {
    double q = sqrt(fmax(q2 * muf, 1.69)), t = ta_eps09s(197, s);
    for (int k = 1; k <= 31; k++) {
      eps09s_c(2, k, 197, x, q, cuv, cdv, cu, cd, cs, cc, cb, cg);
      printf(" %.6f", 1.0 + t * (cg[0] + t * (cg[1] + t * (cg[2] + t * cg[3]))));
    }
    printf("\n");
  }
  return 0;
}
