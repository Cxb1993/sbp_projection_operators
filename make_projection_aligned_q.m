% MAKE_PROJECTION_ALIGNED_q create a set of coefficients and null space for
% projeciont operator
% [q, nA] = make_projection_aligned_q(p,m,pb,r,p_SBP)
%
% inputs:
%   p:     order of polynomials on the glue in the interior
%   m:     number of intervals in the interior
%   pb:    order of polynomials on the glue at the boundary
%   r:     boundary width
%   p_SBP: initerior SBP finite difference accuracy
%
% outputs:
%   q:  a set of coefficients for the projection operator
%   nA: nullspace of the constraints system, i.e. q + nA*v is a solution
%       projection
function [q, nA] = make_projection_aligned_q(p,m,pb,r,p_SBP)
  s = r-(m/2-1); % note that this is s+1 in the paper

  %% Get the SBP operator
  switch(p_SBP)
  case(2)
    N = 3 + r + m;
  case(4)
    N = 9 + r + m;
  case(6)
    N = 13 + r + m;
  case(8)
    N = 17 + r + m;
  case(10)
    N = 23 + r + m;
  otherwise
    error('Invalid SBP order, must be 2,4,6,8,10')
  end

  [~,HI] = diagonal_sbp(p_SBP, N);
  H = diag(sparse(1./diag(HI)));

  xf_i = 0:m; % fd grid
  xg_i = 0:m; % glue grid
  [U_i, V_i, Mr_i] = glue_pieces(p, xf_i, xg_i);

  %%
  %% Interior
  %%

  %% glue to grid constraints:
  Ag2f_ii = U_i;
  Ag2f_ib = zeros((p+1),r*s*(p+1));
  bg2f_i  = V_i(m/2+1,:)';


  %% grid to glue constraints
  Af2g = V_i(1:m,:);
  Bf2g = Mr_i*U_i(:,(p+1)*(m/2-1)+(1:p+1))';

  Af2g_ii = kron(fliplr(Af2g'),fliplr(eye(p+1)));
  Af2g_ib = zeros((p+1)^2,r*s*(p+1));
  bf2g_i  = reshape(flipud(Bf2g),numel(Bf2g),1);

  % fliplr(reshape(1:m*(p+1),p+1,m))*Af2g
  % ([Af2g_ii,Af2g_ib]*q_test)'
  % return

  clear Af2g Bf2g

  % solve using grid to glue
  % Q = Bf2g/Af2g
  % q = reshape(fliplr(Q),m*(p+1),1);
  % q1 = Af2g_ii\bf2g_i;
  % disp([q,q1])
  % norm(q1-q)
  % return

  %% Symmetry
  % Even terms are symmetric and odd terms are skew
  Sy = [eye(m/2),-fliplr(eye(m/2))];
  Sk = [eye(m/2), fliplr(eye(m/2))];
  S_ii = kron(Sy,diag(mod(1:p+1,2)))+kron(Sk,diag(mod(0:p,2)));
  S_ib = zeros(m/2*(p+1),r*s*(p+1));
  bs_i = zeros(m/2*(p+1),1);
  clear Sy Sk

  % solve using interior
  % A = [Af2g_ii,Af2g_ib;Ag2f_ii,Ag2f_ib;S_ii,S_ib];
  % b = [bf2g_i;bg2f_i;bs_i];
  % q = A\b;
  % disp(norm(A*q - b))

  %%
  %% boundary
  %%
  xf_b = 0:s+m-2; % fd grid
  xg_b = 0:s+m-2; % glue grid
  [U_b, V_b, Mr_b] = glue_pieces(p, xf_b, xg_b);

  %% glue to grid constraints
  Ag2f_bi = zeros(s*(pb+1),m*(p+1));
  Ag2f_bb = kron(eye(s),U_b(1:pb+1,1:r*(p+1)));
  bg2f_b  = reshape(V_b(1:s,1:pb+1)',s*(pb+1),1);

  % qb  = [1:r*s*(p+1)]';
  % Qb = reshape(qb,r*(p+1),s);
  % [Ag2f_bb*qb,reshape(U_b(1:pb+1,1:r*(p+1))*Qb,size(bg2f_b))]

  %% grid to glue
  Af2g_b = V_b(1:s,1:pb+1)'*H(1:s,1:s);
  Af2g_i = V_b(s+1:s+m-1,1:pb+1)'*H(s+1:s+m-1,s+1:s+m-1);
  Bf2g   = U_b(1:pb+1,1:r*(p+1))*kron(eye(r),Mr_b(1:p+1,1:p+1));

  % stack
  bf2g_b  = reshape(Bf2g',numel(Bf2g),1);
  Af2g_bb = kron(Af2g_b,eye(r*(p+1)));

  % Q = Af2g_b \ Bf2g
  % q = Af2g_bb \ bf2g_b
  % disp(norm([reshape(Q',prod(size(Q)),1)-q]))
  T = [zeros(r-(m-1),pb+1);Af2g_i'];
  L = zeros((pb+1)*r,m);
  for i = 1:m
    L(:,i) = T(:);
    T = [zeros(1,pb+1);T(1:end-1,:)];
  end
  Af2g_bi = kron(L,eye(p+1));

  I = [zeros(1,(p+1)*(r-(m-1))),1:(m-1)*(p+1)]';
  Qi = zeros(r*(p+1),m-1);
  for i = 1:m-1
    Qi(:,i) = I;
    I = [zeros(p+1,1);I(1:end-(p+1))];
  end

  %% constraint check
  % qi = [1:m*(p+1)]';
  % qb  = [1:r*s*(p+1)]';
  % Qb = reshape(qb,r*(p+1),s);
  % [Af2g_bb*qb-reshape([Af2g_b*Qb']',size(bf2g_b)),Af2g_bi*qi-reshape([Af2g_i*Qi']',size(bf2g_b))]

  clear Af2g_b Af2g_i Bf2g qi qb Qb Qi I

  % solve using everything
  As_i   = [S_ii,S_ib];
  Ag2f_i = [Ag2f_ii,Ag2f_ib];
  Af2g_i = [Af2g_ii,Af2g_ib];
  Ag2f_b = [Ag2f_bi,Ag2f_bb];
  Af2g_b = [Af2g_bi,Af2g_bb];

  A = [As_i;Ag2f_i;Af2g_i;Ag2f_b;Af2g_b];
  b = [bs_i;bg2f_i;bf2g_i;bg2f_b;bf2g_b];

  warning('off','MATLAB:rankDeficientMatrix');
  q = A \ b;
  warning('on','MATLAB:rankDeficientMatrix');
  fprintf(['\nSBP order:                 %d'...
           '\nInterior polynomial order: %d'...
           '\nBoundary polynomial order: %d'...
           '\nConstraint check (~zero?): %e\n'],...
          p_SBP,p,pb,norm(A*q - b))
  % spy(A)
  nA = null(full(A));
  return
