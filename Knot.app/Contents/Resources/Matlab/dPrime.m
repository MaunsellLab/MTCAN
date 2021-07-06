function [dP, c, beta] = dPrime(pH, pFA)
%
% Compute d-prime based on hits and false alarms 
%
% convert to Z scores and calculate d-prime and criterion
  zH = -sqrt(2) .* erfcinv(2 * pH);
  zFA = -sqrt(2) .* erfcinv(2 * pFA);
  dP = zH - zFA ;
  c = -0.5 * (zH + zFA);
  beta = exp((zFA^2 - zH^2)/2); 
end
