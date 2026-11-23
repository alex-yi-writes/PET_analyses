function [CT] = simESRTMfixk2p_1_0_0(varargin)
%% simESRTM 1.0.0 (c) Turku PET Centre 2009
% CT = simESRTM(T,CR,nr,R1,k2,BP)
% T  = mid times (length(T)=nr)
% CR = reference region TAC (length(CR)=nr)
% nr = number of time points
% R1 = ratio K1/K1'
% k2 = k2 apparent
% BP = vector of BP in each mid time (length(BP)=nr)

% Array of time values %
t=varargin{1};

% Reference region activities %
cr=varargin{2};

% Number of values in TACs %
nr=varargin{3};
% Check for data %
if(nr<2); fprintf(1,'Not enough time points -> Abort'); return; end;

% Ratio K1/K1' %
R1=varargin{4};

% Rate constant of the model %
k2p=varargin{5};

% Binding potentials %
BP=varargin{6};


%% Calculate Tissue simulation
dt      = t(1) - 0;
cri     = dt*((cr(1)+0)/2);
CT      = R1*(cr(1) + k2p*cri)/(1+(k2p*R1/(1+BP(1)))*(dt/2));
PI      = 0;
cti_last= dt*CT(1)/2;
for i=2:nr
    dt  = t(i)-t(i-1);
    cri = cri + dt*((cr(i)+cr(i-1))/2);
    if BP(i) == BP(i-1)
        CT = [CT; (PI + R1*(cr(i) + k2p*cri) - (k2p*R1/(1+BP(i)))*(cti_last + dt*CT(i-1)/2))/(1+(k2p*R1/(1+BP(i)))*(dt/2))];
        cti_last = cti_last + dt*((CT(i)+CT(i-1))/2);
    else
        PI = PI - (R1*k2p/(1+BP(i-1)))*cti_last;
        cti_last = 0;
        CT = [CT; (PI + R1*(cr(i) + k2p*cri) - (k2p*R1/(1+BP(i)))*(cti_last + dt*CT(i-1)/2))/(1+(k2p*R1/(1+BP(i)))*(dt/2))];
        cti_last = dt*((CT(i)+CT(i-1))/2);
    end
end