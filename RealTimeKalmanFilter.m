function [signal,cfg,model,stimu] = RealTimeKalmanFilter(origin, stimu, cfg, model)
%RealTimeKalmanFilter Using Kalman filter for each point
%   Use as:
%       [signal,cfg,model] = RealTimeKalmanFilter(origin, stimu, cfg, model)
%   Input:
%       - origin, current recorded point
%       - stimu, stimulation ON or OFF (0 or 1)
%       - cfg, parameters of estimation
%       - model, parameters of AR model
%   Output:
%       - signal, current signal without DBS artifacts
%       - cfg, updated cfg parameters
%       - model, updated AR model
%
%
%   Author   : Xuanjun Guo
%   Created  : Jan 30, 2024
%   Modified : Jul 7, 2024

% Sure stimulation duration tpeaks~fixed

%FixedThres = cfg.mean+cfg.OutlierCondition*cfg.std;
FixedThres = 300;

%% Find the artifacts
if stimu > 0
    cfg.jitterTime = cfg.jitterTime+1;
    if abs(origin) > FixedThres
        cfg.start = cfg.idx;
        cfg.fixed = cfg.start + cfg.duration;
        stimu = 0;% detected stimulation
        cfg.artifacts = 1;
        cfg.jitterTime = 0;
        %{
    elseif cfg.jitterTime > fix(18*10^-3*cfg.fs)
        cfg.start = cfg.idx;
        cfg.fixed = cfg.start + cfg.duration;
        stimu = 0;% detected stimulation
        cfg.artifacts = 1;
        cfg.jitterTime = 0;
        %}
    end
end

if cfg.idx >= cfg.start-1 && cfg.idx < cfg.fixed 
        % During artifact
        % Set noise matrices to high values so only model is trusted
        Q = cfg.Qt;
        R = cfg.R;
        u = 0;
        cfg.Remove = 1;
        cfg.artifacts = 1;
        if abs(origin-cfg.lastSignal(1,1)) < FixedThres
            cfg.artifacts = 0;
        end
else
    % During artifact-free sections
    Q = cfg.Qe;
    R = 0;
    u = 0;
    cfg.Remove = 0;
    if abs(origin-cfg.lastSignal(1,1)) > FixedThres
        cfg.artifacts = 1;
    else
        cfg.artifacts = 0;
    end
end

    % Kalman filtering
    x_time = model.A * cfg.lastSignal + model.B * u;
    P_time = model.A * cfg.lastP * model.A' + Q;

    K = P_time * model.H' * pinv(model.H * P_time * model.H' + R);
    epsilon = origin - model.H * x_time - model.D *u;

    xhat = x_time + K * epsilon;
    % Output signal
    signal = xhat(1);


%Update model
cfg.idx = cfg.idx+1;
cfg.lastP = P_time - K * model.H * P_time;
%cfg.lastSignal1 = cfg.lastSignal;
cfg.lastSignal = xhat;
end
