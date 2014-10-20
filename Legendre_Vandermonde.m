% LEGENDRE_VANDERMONDE Generalized Vandermonde from the normalized Legendre
% polynomials
% V = Legendre_Vandermonde(x,N)
%
% input:
%   x: points to evaluate the polynomials at
%   N: number of polynomials to evaluate
%
% output:
%   V: Generalized Vandermonde
%
% Based on routines given in
%   @BOOK{HesthavenWarburton2008,
%     title = {Nodal Discontinuous {G}alerkin Methods: {A}lgorithms, Analysis, and
%    Applications},
%     publisher = {Springer},
%     year = {2008},
%     author = {Hesthaven, Jan S. and Warburton, Tim},
%     volume = {54},
%     series = {Texts in Applied Mathematics},
%     doi = {10.1007/978-0-387-72067-8}
%   }

function V = Legendre_Vandermonde(x,N)
  scale = @(N) sqrt((2*N+1)/2);

  x = x(:);

  P_n_1 = ones(size(x));
  if(N >= 0)
    V(:,1) = P_n_1 .* scale(0);
  end

  P_n_0 = x;
  if(N>=1)
    V(:,2) = P_n_0 .* scale(1);
  end

  for n=2:N
    a = (2*n - 1)/n;
    c = (n - 1)*(n - 1)*(2*n)/(n*(n)*(2*n - 2));
    P_n_2 = P_n_1;
    P_n_1 = P_n_0;
    P_n_0 = a.*x.*P_n_1 - c.*P_n_2;
    V(:,n+1) = P_n_0 .* scale(n);
  end
end
