function [cfg,model,Prediction,p,stimu] = Kalman_Ini(width,Frequency,SR)
% state initialization
stimu = 0;
 % AR model estimation
    p = 2;
    %a = arburg(raw(1:3),p);
    Prediction = zeros(1,3);
    model.A = [1.97131546347023,-0.971449115448927;1,0];%[- a(2:(p+1)); eye(p-1, p)];
    model.B = [1; 0];
    model.H = [1, 0];
    model.D = 1;

    % Setup for all the Kalman variables
    % Noise matrices
    sigmaE = 0.025;
    sigmaT = 0.1;
    sigmaV = 0.25;

    % Initialize Q with initial values
    Qe = diag([sigmaE, sigmaE]);
    Qt = diag([sigmaE, sigmaT]);
    
    cfg.duration = floor(width*SR)+1;
    cfg.frequency = Frequency;
    cfg.fs = SR;
    cfg.start = 0;
    cfg.fixed = 0;
    cfg.delay = 0;
    cfg.t0 = 0;
    cfg.idx = 0;
    cfg.lastSignal = 0;
    cfg.lastSignal1 = 0;
    cfg.lastP = 1;
    cfg.R = sigmaV*60000000;% Set a high value so that artifact will not be trusted
    cfg.Qe = Qe;
    cfg.Qt = Qt;
    cfg.Remove = 0;
    cfg.artifacts = 0;
    cfg.std = 0;
    cfg.mean = 0;
    cfg.jitterTime = 0;
    cfg.OutlierCondition = 5;%2.58;
end

