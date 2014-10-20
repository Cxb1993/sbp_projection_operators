% Create the optimized projection operators

% p:  interior sbp order
% pb: boundary accuracy
% m:  number of intervals used to project to a grid point
% r:  number of intervals used at the boundary
for p = 2:2:10
  if(p == 2)
    pb = 1;
    m  = 2;
    r  = 1;
  elseif(p == 4)
    pb = 2;
    m  = 4;
    r  = 5;
  elseif(p == 6)
    pb = 3;
    m  = 6;
    r  = 8;
  elseif(p == 8)
    pb = 4;
    m  = 8;
    % r  = 11; % 11 is the minimal number, but 12 performs better
    r = 12;
  elseif(p == 10)
    pb = 5;
    m  = 10;
    r  = 16;
  end

  % N: size of finite difference grid to use in the optimization
  N = 64;

  [Pf2g,Pg2f,q] = make_projection_aligned_opt(N,p-1,m,pb-1,r);

  name=['optimal_',num2str(p),'_q'];
  save(name,'q','pb','m','r')
end
