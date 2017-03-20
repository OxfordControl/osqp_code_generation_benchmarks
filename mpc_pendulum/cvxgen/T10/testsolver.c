/* Produced by CVXGEN, 2017-03-20 02:49:33 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"
Vars vars;
Params params;
Workspace work;
Settings settings;
#define NUMTESTS 0
int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data();
  /* Solve problem instance for the record. */
  settings.verbose = 1;
  num_iters = solve();
#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;
  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();
  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;
  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif
  return 0;
}
void load_default_data(void) {
  params.x_0[0] = 0.20319161029830202;
  params.x_0[1] = 0.8325912904724193;
  params.x_0[2] = -0.8363810443482227;
  params.x_0[3] = 0.04331042079065206;
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.Q[0] = 1.8929469543476547;
  params.Q[4] = 0;
  params.Q[8] = 0;
  params.Q[12] = 0;
  params.Q[1] = 0;
  params.Q[5] = 1.896293088933438;
  params.Q[9] = 0;
  params.Q[13] = 0;
  params.Q[2] = 0;
  params.Q[6] = 0;
  params.Q[10] = 1.1255853104638363;
  params.Q[14] = 0;
  params.Q[3] = 0;
  params.Q[7] = 0;
  params.Q[11] = 0;
  params.Q[15] = 1.2072428781381868;
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.R[0] = 1.0514672033008299;
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.QN[0] = 1.4408098436506365;
  params.QN[4] = 0;
  params.QN[8] = 0;
  params.QN[12] = 0;
  params.QN[1] = 0;
  params.QN[5] = 1.0298762108785668;
  params.QN[9] = 0;
  params.QN[13] = 0;
  params.QN[2] = 0;
  params.QN[6] = 0;
  params.QN[10] = 1.456833224394711;
  params.QN[14] = 0;
  params.QN[3] = 0;
  params.QN[7] = 0;
  params.QN[11] = 0;
  params.QN[15] = 1.6491440476147607;
  params.A[0] = -0.8860508694080989;
  params.A[1] = 0.7050196079205251;
  params.A[2] = 0.3634512696654033;
  params.A[3] = -1.9040724704913385;
  params.A[4] = 0.23541635196352795;
  params.A[5] = -0.9629902123701384;
  params.A[6] = -0.3395952119597214;
  params.A[7] = -0.865899672914725;
  params.A[8] = 0.7725516732519853;
  params.A[9] = -0.23818512931704205;
  params.A[10] = -1.372529046100147;
  params.A[11] = 0.17859607212737894;
  params.A[12] = 1.1212590580454682;
  params.A[13] = -0.774545870495281;
  params.A[14] = -1.1121684642712744;
  params.A[15] = -0.44811496977740495;
  params.B[0] = 1.7455345994417217;
  params.B[1] = 1.9039816898917352;
  params.B[2] = 0.6895347036512547;
  params.B[3] = 1.6113364341535923;
  params.x_max[0] = 1.6915017425863585;
  params.x_max[1] = 0.7559880826577783;
  params.x_max[2] = 0.18443401774344848;
  params.x_max[3] = 1.3068218050470723;
  params.u_max[0] = 1.1156815247769019;
}
