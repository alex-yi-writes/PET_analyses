function [CT] = simSRTM2_1_0_0(t, CR, nr, k2p, R1, BP)
% simSRTM2_1_0_0:
%   Simulate tissue TAC using the Simplified Reference Tissue Model (SRTM2),
%   where k2' is fixed and only 2 parameters (R1, BP) are being fitted.
%
% USAGE:
%   CT = simSRTM2_1_0_0(t, CR, nr, k2p, R1, BP)
%
% INPUTS:
%   t    - time midpoint array, length = nr
%   CR   - reference region TAC, length = nr
%   nr   - number of time points
%   k2p  - fixed reference efflux rate constant (k2')
%   R1   - ratio K1/K1'  (fitted externally)
%   BP   - binding potential (fitted externally)
%
% OUTPUT:
%   CT   - model-predicted target‐region TAC (length = nr)
%
% -------------------------------------------------------------------------

%% STEP 1: initialize integrated values at the first frame
% interpret t(1) as the midpoint of the first PET frame, so from time=0
% up to t(1), we do a half-trapezoid approximation for integrated reference
dt       = t(1);               % from 0 to the midpoint of the 1st frame
cri      = dt*(CR(1)/2);       % integrated reference from 0->t(1)
CT_init  = R1 * CR(1);         % approximate tissue conc. at the first midpoint
cti_last = dt*(CT_init/2);     % integrated tissue from 0->t(1), half-trapezoid

CT = CT_init;                  % store first TAC value in the output

%% STEP 2: recursively calculate CT for each subsequent frame
for i = 2:nr
    % frame i goes from t(i-1) to t(i). 
    % dt is the duration of that step:
    dt = t(i) - t(i-1);

    % incrementally integrate the reference TAC up to time t(i), again
    % via the trapezoid rule:
    cri = cri + dt * (CR(i) + CR(i-1))/2;

    % --------------------------------------------------------
    % the discrete SRTM2 update step:
    %
    %   C_T(i) = [ R1 * ( CR(i) + k2p * (integrated CR up to i) )
    %              - (k2p*R1 / (1+BP)) * ( integrated tissue + partial step )
    %            ] / [ 1 + (k2p*R1 / (1+BP)) * (dt/2) ]
    %
    % this is effectively a solution to the simplified reference‐tissue
    % model differential equation, handled via a trapezoidal integration
    % approach (same style as the previous ESRTM code).
    % --------------------------------------------------------
    numerator = R1*( CR(i) + k2p*cri ) ...
              - (k2p*R1/(1+BP))*( cti_last + dt*CT(i-1)/2 );
    denominator = 1 + (k2p*R1/(1+BP))*( dt/2 );
    CTi = numerator / denominator;

    % update integrated tissue up to this new frame boundary:
    cti_last = cti_last + dt * (CTi + CT(i-1))/2;

    % append to the output TAC array
    CT = [CT; CTi];
end
end
