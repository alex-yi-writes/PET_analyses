function CR = sim1TC_ref(t, nr, K1, k2)
% sim1TC_ref:
%   simulate the reference region as a single-tissue compartment
%   ignoring vascular component, i.e. dC/dt = K1 * Cp - k2 * C.
%   here we approximate Cp ~ 1 for scaling or incorporate the shape 
%   if you have an actual plasma input (see notes).
%
%   t   : vector of mid times
%   nr  : number of frames
%   K1  : uptake constant
%   k2  : efflux constant
%
% OUTPUT:
%   CR  : model-predicted reference TAC
% 
% AUTHOUR: Yeo-Jin Yi (yeo-jin.yi.15@ucl.ac.uk)

CR = zeros(nr,1);

% for simplicity assume Cp=1.0 "constant". if you have a real 
% measured plasma input, you would incorporate that... then 
% dCR/dt = K1*1 - k2*CR. we do an euler/trapezoid integration

dt = t(1);  % from 0->t(1)
CR(1) = K1/k2 * (1 - exp(-k2*dt));  % approximate the first frame

for i=2:nr
    dt = t(i) - t(i-1);
    % one-step approximation
    CR(i) = CR(i-1) + dt*( K1 - k2*CR(i-1) );
end
end
