This directory contains the files necessary to construct higher, SBP projection
operators

# `diagonal_sbp.m`
Constructs diagonal SBP operators with interior order 2, 4, 6, 8, and 10

# `make_opt.m`
Generates the coefficients for the 'optimized' SBP projections

# `sbp_sbp_projection.m`
Example of how to build a projection operator between two SBP multiblocks

# `sbp_sbp_test.m`
Tests to demonstrate the accuracy of the sbp_sbp_projection code

# Cached optimized coefficients

 - `optimal_10_q.mat`
 - `optimal_2_q.mat`
 - `optimal_4_q.mat`
 - `optimal_6_q.mat`
 - `optimal_8_q.mat`

# Helper functions

 - `Gaussian_quad.m`
 - `Legendre_Vandermonde.m`
 - `glue_pieces.m`
 - `make_projection.m`
 - `make_projection_aligned_P.m`
 - `make_projection_aligned_opt.m`
 - `make_projection_aligned_q.m`
 - `make_projection_g2g_h_gen.m`
 - `make_projection_g2g_hr.m`
