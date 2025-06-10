function [x,y,factor,nu,C1a,C2a,C3a,enuDela,ealDela,etaa,gama] = NROAmp_Ini(Frequency,SR)
%% NRO intialization
    x = 0;
    y = 0;

    nu=2*pi*Frequency;     % initial guess for frequency (10% large than the true value)
    dt = 1/SR;
    % parameters of the non-resonant "device"
    om0=5*nu;  om0_2=om0*om0;         % oscillator frequency
    %alpha_a=80;
    alpha_a=0.75*nu; gama=alpha_a/2;       % damping for the "amplitude device"     
    factor=sqrt((om0_2-nu*nu)^2+(alpha_a*nu)^2);
    % precomputed coefficients for oscillators
    [C1a,C2a,C3a,enuDela,ealDela,etaa]=ExSolCoefs(om0,dt,gama);  % for the "amplitude device"
end

