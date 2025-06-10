%% instantiate the library
clear
clc

stLabJack = 1;
Level0=0;
LevelAmp = 0;% Initialization
nDac=0;  % 0 or 1

if stLabJack > 0
    [stLabjack,ljudObj,ljhandle,ljchannel]=Labjack_Ini();
    ljudObj.eDAC (ljhandle, nDac, Level0, 0, 0, 0);
    if stLabjack <=0; return; end
end
%% Parameter setting (double Check)
% stimulation duration / second
tDur=10;
% Target frequency band / Hz
Frequency = 8;
% Target Phase (must be -180 to 180)
% Cortical-alpha: 0 90 180 -90
% STN beta: 0 180
TargetPhase = 90;

%% Parameter setting (Need to be Manully set before and here is just a guess)

% Jitter Delay / second
delay = 0.0328;

% Amplitude Threshold (not noise)
THRS = 0.5;

% Artifacts detection threshold
%% Parameter setting (Normally don't need to Change)
tStimDur=0.008;% Initialized stimulation srtifacts duration / second
tLjTtlWidth=0.5*10^-3;
stUpdateStimDur=1;
% AC marker amplitude
switch TargetPhase 
    case 0
        LevelAmp = 0.05; 
    case 90
        LevelAmp = 0.1; 
    case 180 
        LevelAmp = 0.15;
    case -90 
        LevelAmp = 0.2;
    otherwise
        error("TargetPhase is not right");
end
%%
disp('Loading the library...');
lib = lsl_loadlib();

disp('Load library is OK');
%  Get possible streams

streams = lsl_resolve_all(lib);

% Display information about each stream
disp('Available LSL streams:');
nStreams=0;
for i = 1:length(streams)
    disp(['Stream ', num2str(i)]);
    disp(['  Name: ', char(streams{i}.name())]);
    disp(['  Type: ', char(streams{i}.type())]);
    disp(['  Channel Count: ', num2str(streams{i}.channel_count())]);
    disp(['  Sampling Rate: ', num2str(streams{i}.nominal_srate()), ' Hz']);
    disp(['  Unique Identifier: ', char(streams{i}.uid())]);
    disp(['  Source Id: ', char(streams{i}.source_id())]);
    disp(' ');

    disp(['  Desc: ', char()]);
    nStreams=nStreams+1;
end
disp("Number of Streams available = "+string(nStreams));
if nStreams <=0 ; return; end
%%
istream=1;
if nStreams > 1; istream=2; end
xstream=streams{istream};    
sName  =xstream.name();
sChans =xstream.channel_count();
sType  =xstream.type();

inlet = lsl_inlet(xstream);
inf   = inlet.info();
%fprintf('The stream''s XML meta-data is: \n');
%fprintf([inf.as_xml() '\n']);
%fprintf(['The cap circumference is: ' inf.desc().child('cap').child_value('size') '\n']);
%%
nChans=inf.channel_count;
SR    =inf.nominal_srate;

fprintf("Stream Name is    : %s\n", sName);
fprintf("The manufacturer  : %s\n", inf.desc().child_value('manufacturer'));
fprintf("SR                : %f\n", SR);
fprintf("Number Of Chans   : %d\n", nChans);
ch = inf.desc().child('channels').child('channel');
clear ChansInfo
fprintf("  N   Titl    Unit          Type\n");
for k = 1:inf.channel_count()
    ChansInfo.ttl(k)={ch.child_value('label')};
    ChansInfo.uni(k)={ch.child_value('unit')};
    ChansInfo.typ(k)={ch.child_value('type')};

    fprintf("%3d  %4s",k,ChansInfo.ttl{k});
    fprintf("  %8s",ChansInfo.uni{k});
    fprintf("  %8s\n",ChansInfo.typ{k});

    ch = ch.next_sibling();
end

% if nChans ~= 1; error("Number of chans must be = 1"); end  xxxxxx

%txt = input("Press any key to continue:  ","s");

%% initialization algorithm

[cfg,model,Prediction,p,stimu] = Kalman_Ini(tStimDur,Frequency,SR);
[u,v,nu,C1p,C2p,C3p,enuDelp,ealDelp,etap,gamp] = NROPhase_Ini(Frequency,SR);
[x,y,factor,~,C1a,C2a,C3a,enuDela,ealDela,etaa,gama] = NROAmp_Ini(Frequency,SR);

%% Storage data
s = zeros(1,3); % NRO eestimation
windowSize = floor(1.0*SR); %
UpdateBufferSize = 5*SR;
buffer = []; % store filter data
durationBuffer = zeros(1,UpdateBufferSize);
Phasebuffer = zeros(1,UpdateBufferSize);
Amplitudebuffer = zeros(1,UpdateBufferSize);

%% adaptive parameter(NRO, AR model,stimulation duration)
%%Stimulation duration update
updatepoint = UpdateBufferSize;
updatestep = UpdateBufferSize;
% parameters of NRO adaptation algorithm
% precomputed quantities for linear fit for frequency adaptation
tbuf=(1:UpdateBufferSize)/SR;  Sx=sum(tbuf); denom=UpdateBufferSize*sum(tbuf.^2)-Sx*Sx; 

%% Filter design
N      = 4;  % Order
Fpass1 = Frequency-0.1*Frequency;  % First Passband Frequency
Fpass2 = Frequency+0.1*Frequency;  % Second Passband Frequency
Apass  = 1;   % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, SR);
Hd = design(h, 'cheby1');
biquad = dsp.BiquadFilter('Structure', 'Direct form II','SOSMatrix',Hd.sosMatrix,'ScaleValues',Hd.ScaleValues);

%% Recording for test
RowDta = zeros(1,tDur*SR);
PhaDta = zeros(1,tDur*SR);
KlaDta = zeros(1,tDur*SR);
FltDta = zeros(1,tDur*SR);
SmpDta = zeros(1,tDur*SR);
TtlDta = zeros(1,tDur*SR);
StmDta = zeros(1,tDur*SR);
InsFrq = zeros(1,tDur*SR);
MeaDta = zeros(1,tDur*SR);
StdDta = zeros(1,tDur*SR);

flag = 0;
PreviousStimu = 0;
n = 20;
PreviousData = 0;

stCheck=0;
stStimReal=-1;

degree_ = 0;
amplitude = 0;
%% check where delay is from, just check the simple stimulation pulse, whether there wil be delay
%while true
isample=0;
nKeep=1000;
ikeep=0;  % 1000;
phaseOld=0;
iStmPha = 0;
stStimTesting=1; stArtefact=-1; nArts=5; iArts=0; 


xTargetPhase = TargetPhase-(90+fix((delay*Frequency)*360));  %   -90+90; -> = 0 deg
if xTargetPhase > 180
    xTargetPhase = xTargetPhase - 360;
else
    if xTargetPhase < -180
        xTargetPhase = xTargetPhase + 360;
    end
end

while isample < (tDur+2)*SR
    % Start output marker

    ljudObj.eDAC (ljhandle, nDac, LevelAmp, 0, 0, 0);
    
    % get chunk from the inlet
    [chunk,stamps] = inlet.pull_chunk();    
    nSamples=length(stamps);
    if nSamples > 0
        %tmx=(1.0*isample)/SR;
        %stx=sprintf(" %6.2f, nSamples= %d",tmx,nSamples);
        %disp(stx);
        if stArtefact > 0
            for iArts=1:5
                chunk(1,iArts)=chunk(1,iArts)+1000; 
            end
        end
        for is=1:nSamples
            raw = chunk(1,is);

%            if stArtefact > 0 && iArts < nArts
%                raw=raw+1000; 
%                iArts=iArts+1;
%            end
            if iArts >= nArts; stArtefact=-1; iArts=0; end

            isample = isample+1;
            SmpDta(isample) = nSamples;
  
        if isample < p+1
            %% memory data for AR model
            Prediction = shiftArray(Prediction, raw,p+1);
            s = shiftArray(s,raw,3);
            buffer = shiftArray(buffer,raw,windowSize);
            cfg.std = std(buffer);
            cfg.mean = mean(buffer);
        elseif isample == p+1
            %% AR model parameter
            Prediction = shiftArray(Prediction, raw,p+1); 

            if raw ~= Prediction(end-1)
                a = arburg(Prediction(1:p+1),p);
                model.A = [- a(2:(p+1)); eye(p-1, p)];
            end

            s = shiftArray(s,raw,3);
            buffer = shiftArray(buffer,raw,windowSize);
            cfg.std = std(buffer);
            cfg.mean = mean(buffer);
            MeaDta(isample) = cfg.mean;
            StdDta(isample) = cfg.std;
        else

            %% Kalman filter to remove artifacts
            %stCheck=0;            
            [signal,cfg,model,stCheck] = RealTimeKalmanFilter(raw, stCheck, cfg, model);
            durationBuffer = shiftArray(durationBuffer,cfg.artifacts,UpdateBufferSize);
            %% AR model parameter
            Prediction = shiftArray(Prediction,signal,p+1);

            if signal ~= Prediction(end-1) && (isnan(signal) == 0)
                a = arburg(Prediction(1:p+1),p);
                model.A = [- a(2:(p+1)); eye(p-1, p)];
            end

            %% Shift Filter 
            buffer = shiftArray(buffer,signal,windowSize);
            cfg.std = std(buffer);
            cfg.mean = mean(buffer);
            MeaDta(isample) = cfg.mean;
            StdDta(isample) = cfg.std;

            if length(buffer) == windowSize
                %% Filter
                ModBuffer = double(buffer);
                Filtered = step(biquad,ModBuffer);
                bp = Filtered(end);
                %% Phase Setimation
                s = shiftArray(s,bp,3);
                [u,v] = OneStep(u,v,gamp,etap,enuDelp,ealDelp,C1p,C2p,C3p,s(1),s(2),s(3));
                z = v/nu;
                phase_ = atan2(-z,u); Phasebuffer = shiftArray(Phasebuffer,phase_,UpdateBufferSize);
                degree_ = rad2deg(phase_);
                
                [x,y] = OneStep(x,y,gama,etaa,enuDela,ealDela,C1a,C2a,C3a,s(1),s(2),s(3));
                z = y/nu;  amplitude_=factor*sqrt(z*z+x*x); % new amplitude value
                Amplitudebuffer = shiftArray(Amplitudebuffer,amplitude_,UpdateBufferSize);
                
                KlaDta(isample) = signal;
                FltDta(isample) = bp;
                RowDta(isample) = raw;
                TtlDta(isample) = flag; flag = 0;
                StmDta(isample) = stCheck;
                PhaDta(isample) = phase_;
                InsFrq(isample) = (phase_-phaseOld)*SR/(2*pi);
                phaseOld=phase_;
            end                
        end                     
            if isample>windowSize*1/10*n  %% Find Target Phase
                if iStmPha < 1 % only have amplitude threshold at the first stimulation
                   if (abs(degree_-xTargetPhase)<(360*Frequency/SR/2)  && amplitude_>THRS ) && (cfg.Remove < 1) && (stCheck < 1 ) && (PreviousStimu < 1)
                       % cfg.Remove <1 :Not in artifacts
                       % stCheck < 1: Already detected artifacts
                       % PreviousStimu < 1: Do not give multiple pulses
                       if stLabJack > 0
                            state = 1;
                            flag = 1;
                            ljudObj.eDO(ljhandle, ljchannel, state);
                            pause(tLjTtlWidth);
                            state = 0;
                            ljudObj.eDO(ljhandle, ljchannel, state);
                            PreviousStimu = 1;
                            stCheck=1;% wait for give stimulation sign
                            tSamplesPass=1.0*is/SR;
                            iStmPha=iStmPha+1;
                            stx=sprintf(" istim= %d, tSamplesPass= %6.3f",iStmPha,tSamplesPass);
                            disp(stx)
                       else
                            if stStimTesting > 0
                                state = 1;
                                flag = 1;
                                stArtefact=1;
                                pause(tLjTtlWidth);
                                PreviousStimu = 1;
                                stCheck=1;% wait for give stimulation sign
                            end
                       end
                   else
                        PreviousStimu = 0;                                              
                   end
                else
                    if (abs(degree_-xTargetPhase)<(360*Frequency/SR/2)) && (cfg.Remove < 1) && (stCheck < 1 ) && (PreviousStimu < 1)
                       % cfg.Remove <1 :Not in artifacts
                       % stCheck < 1: Already detected artifacts
                       % PreviousStimu < 1: Do not give multiple pulses
                       if stLabJack > 0
                            state = 1;
                            flag = 1;
                            ljudObj.eDO(ljhandle, ljchannel, state);
                            pause(tLjTtlWidth);
                            state = 0;
                            ljudObj.eDO(ljhandle, ljchannel, state);
                            PreviousStimu = 1;
                            stCheck=1;% wait for give stimulation sign
                            tSamplesPass=1.0*is/SR;
                            iStmPha=iStmPha+1;
                            stx=sprintf(" istim= %d, tSamplesPass= %6.3f",iStmPha,tSamplesPass);
                            disp(stx)
                       else
                            if stStimTesting > 0
                                state = 1;
                                flag = 1;
                                stArtefact=1;
                                pause(tLjTtlWidth);
                                PreviousStimu = 1;
                                stCheck=1;% wait for give stimulation sign
                            end
                       end
                   else
                        PreviousStimu = 0;                                              
                   end
                end
        
             end
        
             end
        
        end % end of nSamples
        %% Update stimulation duration estimation && AR model
        %%
        if (isample>updatepoint)
            if amplitude_>THRS
                % AR model
                a = arburg(double(buffer(end-p:end)),p);
                model.A = [- a(2:(p+1)); eye(p-1, p)];
            end            
%{
            % artifact duration
            dx = diff(durationBuffer);
            rising  = find(dx > 0);
            falling = find(dx < 0);
            if length(rising) >= 1 && length(falling) >= 1
                if rising(1) > falling(1)
                falling(1) = [];
                end
                if rising(end) > falling(end)
                    rising(end) = [];
                end
                durationList = (falling - rising);
                if stUpdateStimDur > 0
                    cfg.duration = round(max(durationList))+1; 
                end
                updatepoint=updatepoint+updatestep;  % ??????
            else 
                updatepoint=updatepoint+updatestep;  % ??????
            end
%}
        end    

end

% Stop ouput marker
ljudObj.eDAC (ljhandle, nDac, Level0 , 0, 0, 0);
%%
StmIndex = find(TtlDta);
nStms=length(StmIndex);
figttl=sprintf("Freq=%4.1f, Phase=%4.1f, iStmPhaPhases= %d, nStms= %d",Frequency,(TargetPhase),iStmPha,nStms);
THRSNew = prctile(abs(FltDta),50);
figttl1=sprintf("Noise Threshold = %4.1f",THRSNew );

InsFrqOffline = diff(unwrap(PhaDta))*SR/(2*pi);
InsFrqOffline(end+1)=0;

CleanMean = mean(KlaDta);
CleanStd = std(KlaDta);
Clean = KlaDta(abs(KlaDta-CleanMean)<2.58*CleanStd);

%ylm=max(Clean);
%ylm=ylm+0.1*ylm;
ylm = 300;

stPlotDta=1;
if stPlotDta > 0
fig = figure;
%set(fig,'Position', [2824 303 904 902]);
%set(fig,'Position', [2249 105 663 764]); % Beta Lab
%set(fig,'Position', [323 199 702 699]); % Alpha Lab
set(fig,'Position', [18.1429 124.4286 701.7143 699.4286]);%Xuanjun Laptop

ndx=length(RowDta);
tmx=(1:ndx)/SR;

ax1 = subplot(5,1,1);hold on;
plot(tmx,RowDta,'-k');grid on;
plot(tmx,KlaDta,'-g');grid on;
legend('RawDta', 'KlaDta'); legend('boxoff');
ylim([-ylm ylm]);
xlim([tmx(1) tmx(end)]);
title(figttl)


ax2 = subplot(5,1,2);hold on;
plot(tmx,KlaDta,'-k');grid on;
plot(tmx,FltDta,'-g');grid on;
plot(tmx,THRSNew*ones(1,length(tmx)),'-r');grid on;
ylim([-ylm ylm]);
legend('KlaDta', 'FltDta','THRS'); legend('boxoff');
xlim([tmx(1) tmx(end)]);
title(figttl1);

ax3 = subplot(5,1,3);hold on;
yyaxis left
plot(tmx,FltDta,'-b');grid on;
ylim([-ylm ylm]);
ax3.YColor = 'b';

yyaxis right;
plot(tmx,TtlDta,'-g');grid on;
ylim([0 1.1]);
xlim([tmx(1) tmx(end)]);
legend('FltDta','Stim'); legend('boxoff');
xlim([tmx(1) tmx(end)]);
ax3.YColor = 'g';

ax4 = subplot(5,1,4);hold on;
yyaxis left
plot(tmx,PhaDta,'-b');grid on;
ax4.YColor = 'b';
yyaxis right;
plot(tmx,InsFrq,'-r');grid on;
plot(tmx,InsFrqOffline,'-g');grid on;
legend('PhaDta','InsFrq','InsFrqOff'); legend('boxoff');
xlim([tmx(1) tmx(end)]);
ax4.YColor = 'r';
ylim([Frequency-3 Frequency+3]);

ax5 = subplot(5,1,5);hold on;
plot(tmx,SmpDta,'-k'); grid on;
legend('Samples'); legend('boxoff');
title('LSL samples')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
xlim([tmx(1) tmx(end)]);
end


%% Aligned to TTLpulse

if nStms <=0; warning("No Stimulation is found"); return; end

clear RawMat KlaMat FltMat PhaMat TtlMat;
period = fix(SR/Frequency);
nst=fix(1.2*period);
nen=fix(1.2*period);
dtm=1.0/SR;
tmx_=(-nst*dtm):dtm:(nen*dtm);
% Artifacts Index
peaks = findpeaks(RowDta);
ampl = max(peaks);
thr = 0.5*ampl;
Index_ = [];
idx = 1;
Stimidx = 1;
gap = floor(0.25*SR/Frequency);
while idx < length(RowDta)
    if RowDta(idx) > thr
        Index_ = [Index_ idx];
        Stimidx = Stimidx+1;
        idx = idx + gap;
    end
    idx = idx + 1;
end

if length(Index_)>length(StmIndex)
    error("Stimulation artifacts is abnormally small, please check the connection");
end
DelayList = (Index_ - StmIndex(1:length(Index_)))./SR;
%OutlierList = (RowDta(Index_+cfg.duration-3) - MeaDta(Index_+cfg.duration-3))./StdDta(Index_+cfg.duration-3);
%Outlier = min(abs(OutlierList));

for k = 1:nStms
    start = StmIndex(k) - nst;
    final = StmIndex(k) + nen;
    if start > 0 && final < length(RowDta)
        RawMat(k,:) = RowDta(start:final);
        KlaMat(k,:) = KlaDta(start:final);
        FltMat(k,:) = FltDta(start:final);
        PhaMat(k,:) = PhaDta(start:final);
        TtlMat(k,:) = TtlDta(start:final);
    end
end

fig=figure;
%set(fig,'Position', [3735 297 636 913]); % Home
%set(fig,'Position', [2923 189 730 848]); % Work
%set(fig,'Position', [2918 119 666 752]); % Beta Lab
%set(fig,'Position', [323 199 702 699]); % Alpha Lab
set(fig,'Position', [721 121 701.7143 699.4286]);%Xuanjun Laptop
%{
ax1=subplot(3,1,1);
yyaxis left;grid on;
plot(tmx_,FltMat,'-r');hold on;
ax1.YColor = 'r';
%{
yyaxis right;grid on;
plot(tmx_,PhaMat,'-b');
plot(tmx_,TtlMat,'-k');
legend('FltDta', 'PhaDta','TTLDta'); legend('boxoff');
ax1.YColor = 'b';
%}
legend('FltDta'); legend('boxoff');
xlim([tmx_(1) tmx_(end)]);
xlabel('Time/ms')
title(figttl);
%}

ax2=subplot(2,1,1); hold on;
yyaxis left;grid on;
plot(tmx_,TtlMat,'-g');
ax2.YColor = 'g';
%{
yyaxis right;grid on;
plot(tmx_,TtlMat,'-k');
ax2.YColor = 'k';
%}
legend('TTLDta'); legend('boxoff');
xlim([tmx_(1) tmx_(end)]);
xlabel('Time/ms')
figttl2=sprintf("Align to TTL pulse( Mean Delay = %f)",mean(DelayList));%, Artifact Duration = %f , Outlier = %4.1f)",mean(DelayList))%,cfg.duration/SR,Outlier);
title(figttl2)


ax3=subplot(2,1,2); hold on;
yyaxis left;grid on;
plot(tmx_,RawMat,'-k');
ax3.YColor = 'k';
ylim([-ylm ylm]);
%{
yyaxis right;grid on;
plot(tmx_,TtlMat,'-b');
ax3.YColor = 'b';
%}
legend('RowDta'); legend('boxoff');
xlim([tmx_(1) tmx_(end)]);
xlabel('Time/ms')

%% Aligned to artifacts
%{


clear RawMat KlaMat FltMat PhaMat TtlMat;
period = fix(SR/Frequency);
nst=fix(1.2*period);
nen=fix(1.2*period);
dtm=1.0/SR;
tmx_=(-nst*dtm):dtm:(nen*dtm);
%{
dx = diff(StmDta);
falling = find(dx < 0);
Index_ = falling;
%}
peaks = findpeaks(RowDta);
ampl = max(peaks);
thr = 0.7*ampl;

Index_ = [];
idx = 1;
Stimidx = 1;
gap = floor(0.75*SR/Frequency);
while idx < length(RowDta)
    if RowDta(idx) > thr
        Index_ = [Index_ idx];
        Stimidx = Stimidx+1;
        idx = idx + gap;
    end
    idx = idx + 1;
end

for k = 1:length(Index_)
    start = Index_(k) - nst;
    final = Index_(k) + nen;
    if start > 0 && final < length(RowDta)
        RawMat(k,:) = RowDta(start:final);
        KlaMat(k,:) = KlaDta(start:final);
        FltMat(k,:) = FltDta(start:final);
        PhaMat(k,:) = PhaDta(start:final);
        TtlMat(k,:) = TtlDta(start:final);
    end
end

fig=figure;
%set(fig,'Position', [3735 297 636 913]); % Home
%set(fig,'Position', [2923 189 730 848]); % Work
%set(fig,'Position', [2918 119 666 752]); % Beta Lab
set(fig,'Position', [323 199 702 699]); % Alpha Lab

ax1=subplot(3,1,1);
yyaxis left;grid on;
plot(tmx_,RawMat,'-k');ylim([-ylm ylm]);
ax1.YColor = 'k';
%{
yyaxis right;grid on;
plot(tmx_,RawMat,'-k');ylim([-ylm ylm]);
legend('FltDta', 'Stimulation'); legend('boxoff');
ax1.YColor = 'k';
%}
legend('Stimulation'); legend('boxoff');
xlim([tmx_(1) tmx_(end)]);
xlabel('Time/ms')
title(figttl);

ax2=subplot(3,1,2); hold on;
yyaxis left;grid on;
plot(tmx_,KlaMat,'-r');
ax2.YColor = 'r';
%{
yyaxis right;grid on;
plot(tmx_,RawMat,'-k');ylim([-ylm ylm]);
ax2.YColor = 'k';
%}
legend('KlaDta'); legend('boxoff');
xlim([tmx_(1) tmx_(end)]);
xlabel('Time/ms')

ax3=subplot(3,1,3); hold on;
yyaxis left;grid on;
plot(tmx_,FltMat,'-b');
ax3.YColor = 'b';
%{
yyaxis right;grid on;
plot(tmx_,RawMat,'-k');ylim([-ylm ylm]);
ax3.YColor = 'k';
%}
legend('FltDta'); legend('boxoff');
xlim([tmx_(1) tmx_(end)]);
xlabel('Time/ms')
title('Align to Artifacts')
%%
%}

%% Save the file
%a = [0.5 0.03];
a = [THRSNew mean(DelayList)];
fid = fopen('Alpha_config.txt','wt');
fprintf(fid,'%f\n',a);
fclose(fid);