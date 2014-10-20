% GLUE_PIECES    create the pieces that we need for the glue
% [U, V, M] = glue_pieces(N, VX, xg)
%
% inputs:
%   N:  order of polynomials on the glue
%   VX: finite difference node locations
%   xg: glue grid interval edges
%
% outputs:
%   U:     glue to finite difference constraints
%   V:     Vandermonde matrix for glue grid intervals defined [-1,1]
%          (finite difference to glue grid constraints)
%   M:     Mass matrix for glue grid

function [U, V, M] = glue_pieces(N, VX, xg)
  VX = shiftdim(VX)';

  K = length(VX) - 1;

  [r, ~] = Gaussian_quad(N);
  va = (1:K)'; vb = (2:K+1)';
  x = ones(N+1,1)*VX(va) + 0.5*(r+1)*(VX(vb)-VX(va));

  M = diag((2./(2*(0:N)+1)));
  Vr  = Legendre_Vandermonde(r, N);

  U = zeros(N+1, K*(N+1));

  Xmin = VX(1);
  Xmax = VX(end);

  xbr = 2*(x-Xmin)./(Xmax-Xmin) - 1;
  Vbr = Legendre_Vandermonde(xbr(:),N) * sqrt(M);

  invSqM_invVr = sqrt(diag(1./diag(M)))/Vr;
  for n=1:N+1
    Pbr = reshape(Vbr(:,n), size(x));
    tmp = invSqM_invVr * Pbr;
    U(n,:) = tmp(:);
  end

  xgbr = 2*(xg-Xmin)./(Xmax-Xmin) - 1;
  V = Legendre_Vandermonde(xgbr,N) * sqrt(M);

  M = M / 2;
end
