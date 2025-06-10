function [u,v,nu,C1p,C2p,C3p,enuDelp,ealDelp,etap,gamp] = NROPhase_Ini(Frequency,SR)
%% NRO intialization
    u = 0;
    v = 0;

    nu=2*pi*Frequency;     % initial guess for frequency (10% large than the true value)
    dt = 1/SR;
    % parameters of the non-resonant "device"
    om0=5*nu;  %om0_2=om0*om0;         % oscillator frequency
    %alpha_a=80; gama=alpha_a/2;       % damping for the "amplitude device"     
    alpha_p=0.1*24*nu;  gamp=alpha_p/2;       % damping for the "phase device"
    %factor=sqrt((om0_2-nu*nu)^2+(alpha_a*nu)^2);
    % precomputed coefficients for oscillators
    %[C1a,C2a,C3a,enuDela,ealDela,etaa]=ExSolCoefs(om0,dt,gama);  % for the "amplitude device"
    [C1p,C2p,C3p,enuDelp,ealDelp,etap]=ExSolCoefs(om0,dt,gamp);  % for the "phase device"
end

