# reduce_dmeson_v2.awk -- node-side reduction of ana/ampt.dat for the D-meson
# framework (D-meson/docs/FRAMEWORK.md sec. 5).
# v2 (2026-10-22): one extra line per event, after the header, with second-order
#  Q-vectors of ALL charged hadrons (no prescale), 0.2 < pT < 2 GeV:
#  "Q nch A_w A_x A_y B_w B_x B_y N_w N_x N_y S_w S_x S_y" (14 fields)
#   A: -1 < eta < -0.05, B: 0.05 < eta < 1, pT-weighted (STAR TPC eta-sub plane)
#   N: 1 < eta < 2.8,    S: -2.8 < eta < -1, unweighted (PHENIX RXN-like plane)
#  Every other line is unchanged; readers that do not know the Q line skip it.
#  header : the 11 stock header fields + Nch(|eta|<0.5) + Nch(|eta|<1)
#  charm  : every charm hadron, stock line verbatim (KF px py pz m x y z t)
#  light  : "KF px py pz m w" with w = 1/prescale
#           K+-, p, pbar   |y|<1.5           w=1
#           pi+-           |y|<1.5           prescale 0.25 (w=4)
#           pi+-, K+-      1.5<=|y|<4        prescale 0.10 (w=10)
# The prescale uses awk's own rand() seeded per job (srand(seed)); it never
# touches the AMPT random sequence.
BEGIN { if (seed == "") seed = 1; srand(seed) }
function flush() { if (nh) { print hdr, nch5, nch10
    printf "Q %d", nq; for (w = 1; w <= 4; w++) printf " %.4f %.5f %.5f", qw[w], qx[w], qy[w]; printf "\n"
    for (i = 1; i <= nb; i++) print buf[i] } nb = 0 }
NF == 11 { flush(); hdr = $0; nh = 1; nch5 = 0; nch10 = 0; nq = 0; for (w = 1; w <= 4; w++) { qw[w] = 0; qx[w] = 0; qy[w] = 0 }; next }
NF >= 9 {
  kf = $1; a = (kf < 0 ? -kf : kf)
  px = $2; py = $3; pz = $4; m = $5
  p = sqrt(px*px + py*py + pz*pz); e = sqrt(p*p + m*m)
  pt = sqrt(px*px + py*py)
  if (pt > 0 && p > pz && e > pz) {
    eta = 0.5 * log((p + pz) / (p - pz)); y = 0.5 * log((e + pz) / (e - pz))
  } else { eta = 99; y = 99 }
  ch = (a == 211 || a == 321 || a == 2212 || a == 3222 || a == 3112 || a == 3312 || a == 3334 || a == 11 || a == 13)
  if (ch && eta > -0.5 && eta < 0.5) nch5++
  if (ch && eta > -1 && eta < 1) nch10++
  if (ch && pt > 0.2 && pt < 2.0) {
    w = 0
    if (eta > -1 && eta < -0.05) { w = 1; u = pt } else if (eta > 0.05 && eta < 1) { w = 2; u = pt }
    else if (eta > 1 && eta < 2.8) { w = 3; u = 1 } else if (eta > -2.8 && eta < -1) { w = 4; u = 1 }
    if (w) { c2 = (px*px - py*py) / (pt*pt); s2 = 2*px*py / (pt*pt); qw[w] += u; qx[w] += u*c2; qy[w] += u*s2; nq++ }
  }
  s = sprintf("%d", a); n = length(s)
  isc = (a > 100 && (substr(s, n-1, 1) == "4" || substr(s, n-2, 1) == "4" || (n >= 4 && substr(s, n-3, 1) == "4")))
  if (isc) { buf[++nb] = $0; next }
  ay = (y < 0 ? -y : y)
  if ((a == 321 || a == 2212) && ay < 1.5) { buf[++nb] = sprintf("%d %.4f %.4f %.4f %.4f 1", kf, px, py, pz, m); next }
  if (a == 211 && ay < 1.5) { if (rand() < 0.25) buf[++nb] = sprintf("%d %.4f %.4f %.4f %.4f 4", kf, px, py, pz, m); next }
  if ((a == 211 || a == 321) && ay < 4) { if (rand() < 0.10) buf[++nb] = sprintf("%d %.4f %.4f %.4f %.4f 10", kf, px, py, pz, m) }
}
END { flush() }
