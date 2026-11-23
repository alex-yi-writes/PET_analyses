function y = gamma_variate_madsen(p,t)
%gamma variate formulation of madsen1992 article

    gamma=p(1);
    alpha=p(2);
    tmax=p(3);
    delay=p(4);

    % idx = find(t > delay);
    % dt  = t(idx) - delay;
    dt    = t - delay;
    dtmax = tmax - delay; 

    % y(idx) = gamma.*(dt./dtmax).^alpha.*exp(alpha.*((tmax-t(idx))./dtmax)); 
    y = gamma.*(dt./dtmax).^alpha.*exp(alpha.*((tmax-t)./dtmax)); 

return