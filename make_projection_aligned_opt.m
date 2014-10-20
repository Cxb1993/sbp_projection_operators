% MAKE_PROJECTION_ALIGNED_OPT make optimized SBP to glue projection matrices
% [Pf2g,Pg2f,q] = make_projection_aligned_opt(N,p,m,pb,r)
%
% inputs:
%   N:  finite difference grid size to use in the optimization
%   p:  polynomial order in the interior and sbp order - 1
%   m:  number of interior intervals
%   pb: boundary polynomial order
%   r:  boundary width
%
% outputs:
%   Pf2g: projection from finite difference to glue for size N
%   Pg2f: projection from glue to finite difference for size N
%   q:    projection coefficients for used in make_projection_aligned_P
function [Pf2g,Pg2f,q] = make_projection_aligned_opt(N,p,m,pb,r)

  % create one set of coefficients and get the nullspace for the operator
  [q, nA] = make_projection_aligned_q(p,m,pb,r,p+1);

  % If no nullspace, return found values to the user
  w = size(nA, 2);
  [Pf2g,Pg2f] = make_projection_aligned_P(N,p,m,r,q,p+1);

  if (w < 1)
    return
  end

  % Optimize the operators so that the eigenvalues are close to flat

  % Display off so that early termination message is suppressed
  opts = optimset('Algorithm', 'levenberg-marquardt','Display','off');
  x = lsqnonlin(@(x) objective(N,p,m,r,q,nA,x),zeros(w,1),[],[],opts);

  % set up the new coefficients and make the operator
  q = q+nA*x;
  [Pf2g,Pg2f] = make_projection_aligned_P(N,p,m,r,q,p+1);
end

% Objective function to optimize such that eigenvalues are as close to a
% flat line as possible.
function [r] = objective(N,p,m,r,q,nA,x)
  q0 = q + nA*x;
  [Pf2g,Pg2f] = make_projection_aligned_P(N,p,m,r,q0,p+1);
  e = sort(real(eig(full(Pg2f*Pf2g))));
  r = diff(e);
  plot(e,'*')
  axis([1 N+1 0 2])
  drawnow
end
