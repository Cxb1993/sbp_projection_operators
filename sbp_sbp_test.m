clear
Na = 100;
Nb = 230;

xa = linspace(-1,1,Na+1)';
xb = linspace(-1,1,Nb+1)';

for qa = 2:2:10
  for qb = 2:2:10
    [Pa2b,Pb2a] = sbp_sbp_projection(Na,qa,Nb,qb);

    % Check the full operator
    for n = 0:min(qa,qb)/2-1
      assert(max(abs(Pa2b*xa.^n-xb.^n)) < 1000*eps)
      assert(max(abs(Pb2a*xb.^n-xa.^n)) < 1000*eps)
    end

    % Check the interio operator
    for n = min(qa,qb)-1
      ta = Pb2a*xb.^n-xa.^n;
      tb = Pa2b*xa.^n-xb.^n;

      % Since the boundary isn't full accuracy we ensure that there are some
      % interior element which succeed by finding the boundary size and only
      % checking the interior elements
      mb = numel(find(abs(tb) > 1000*eps))/2;
      assert(mb < Nb/2)
      assert(max(abs(tb(mb+1:end-mb))) < 1000*eps);

      ma = numel(find(abs(ta) > 1000*eps))/2;
      assert(ma < Na/2)
      assert(max(abs(ta(ma+1:end-ma))) < 1000*eps);
    end
  end
end
