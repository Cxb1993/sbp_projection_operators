% MAKE_PROJECTION_ALIGNED_P make SBP to glue projection matrices
% [Pf2g,Pg2f] = make_projection_aligned_P(N,p,m,r,q,p_SBP)
%
% inputs:
%   N:     finite difference operators size is N+1
%   p:     order of polynomials on the glue in the interior
%   m:     number of intervals in the interior
%   r:     boundary width
%   q:     coefficients of the projection operator
%   p_SBP: initerior SBP finite difference accuracy
%
% outputs:
%   Pf2g:  projection from the finite difference grid to the glue
%   Pg2f:  projection from the glue to the finite difference grid
function [Pf2g,Pg2f,M,H] = make_projection_aligned_P(N,p,m,r,q,p_SBP)
  s = r-(m/2-1);

  %% Get the SBP operator
  [~,HI] = diagonal_sbp(p_SBP, N);
  H = diag(sparse(1./diag(HI)));

  xf_b = 0:s+m-2; % fd grid
  xg_b = 0:s+m-2; % glue grid
  [~, ~, Mr_b] = glue_pieces(p, xf_b, xg_b); % get the glue mass matrix

  % Create the projection from the glue to the finite difference grid
  I = kron((s+1:N-s+1)',ones(m*(p+1),1));
  J = kron(ones(N+1-2*s,1),(1:m*(p+1))')+(p+1)*(I-1-m/2);
  qi = kron(ones(N+1-2*s,1),q(1:m*(p+1)));

  % interior
  Pg2f = sparse(I,J,qi,N+1,N*(p+1));

  % boundary
  Qb = reshape(q(m*(p+1)+1:end),r*(p+1),s)';
  Pg2f(1:s,1:r*(p+1)) = Qb;
  Qb = reshape(diag(2*mod(1:p+1,2)-1)*...
               flipud(reshape(rot90(Qb,2)',p+1,r*s)),r*(p+1),s)';
  Pg2f(end+1-s:end,end+1-r*(p+1):end) = Qb;

  % Projection the back comes from the stability relationship
  M = kron(speye(N),Mr_b);
  Pf2g = kron(speye(N),inv(Mr_b))*Pg2f'*H;
end
