%% instantiate the library
clear
clc

stLabJack = 1;

if stLabJack > 0
    [stLabjack,ljudObj,ljhandle,ljchannel]=Labjack_Ini();
    if stLabjack <=0; return; end
end
%% Parameter setting (double Check)

% Stimulation pulses number
PulseNumber=300;

% Stimulation Interval/ second
ITI = 1;

%% Parameter setting (Normally don't need to Change)
tLjTtlWidth=0.5*10^-3;
stUpdateStimDur=1;


for ix=1:PulseNumber
        
        state = 1;
        ljudObj.eDO(ljhandle, ljchannel, state);
        disp(['FIO7 set to ' num2str(state)]);
        
       
        pause(tLjTtlWidth);
        
        state = 0;
        ljudObj.eDO(ljhandle, ljchannel, state);
        
        pause(ITI+rand(1))
        disp(ix);
end