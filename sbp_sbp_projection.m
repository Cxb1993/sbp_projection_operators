% SBP_SBP_PROJECTION Create projection operators for two SBP blocks which have
% different numbers of grid points
% [Pa2b,Pb2a] = sbp_sbp_projection(Na,qa,Nb,qb)
%
% inputs:
%   Na:   finite difference N on side a
%   qa:   finite difference order on side a
%   Nb:   finite difference N on side b
%   qb:   finite difference order on side b
%
% outputs:
%   Pa2b: projection from side a to side b
%   Pb2a: projection from side b to side a
function [Pa2b,Pb2a] = sbp_sbp_projection(Na,qa,Nb,qb)

% Get the projection operators for each side
[Paf2g, Pag2f]  = make_projection(Na,qa);
[Pbf2g, Pbg2f]  = make_projection(Nb,qb);

% Modify the projection operators so that they go to the same order polynomials,
% which is the max(qa,qb)-1
pg   = max(qa,qb)-1;

Paf2g = kron(speye(Na),speye(pg+1,qa))*Paf2g;
Pag2f = Pag2f*kron(speye(Na),speye(qa,pg+1));

Pbf2g = kron(speye(Nb),speye(pg+1,qb))*Pbf2g;
Pbg2f = Pbg2f*kron(speye(Nb),speye(qb,pg+1));

% To create the projection operators to the same glue space
xa = linspace(-1,1,Na+1);
xb = linspace(-1,1,Nb+1);
[~, Pa2g, Pg2a, Pb2g, Pg2b] = make_projection_g2g_hr(pg, xa, xb);

% Stack the projection operators together
Pa2b = Pbg2f * Pg2b * Pa2g * Paf2g;
Pb2a = Pag2f * Pg2a * Pb2g * Pbf2g;
