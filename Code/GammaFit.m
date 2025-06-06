function [results] = GammaFit(x, TRinSeconds, EvLengthInSeconds)

% Script to fit single gamma functions to experimental data
% Based on scripts from Baxter Rogers
% Modified by Chris Asplund, October 2005
% Turned into a function, 3 November 2005
% Combined with CanonicalGammaInfo by Alvin Wong, October 2018

% Event-related averaging of time series
%
% $Id: eravg.m,v 1.10 2005/06/16 20:38:36 rogersbp Exp $

% $Log: eravg.m,v $
% Revision 1.10  2005/06/16 20:38:36  rogersbp
% Now with fewer parameters!
%

global TR 
global events
global fullhrf

events = x;

%for i = 1:6,
%events = area6(i, 1:20);
TR = TRinSeconds;

%for i = 1:20,
%    res(i) = mean([d1(i) d2(i) d3(i) d4(i) d5(i) d6(i)]);
%end

rct = 1;  % result count, used later

options = optimset('MaxFunEvals',100000);
options = optimset('MaxIter',10000);
    
%%%   delay  amp   ons  disp
x0 = [ 3.5    0.8     0    1]; 
lb = [  0    0.01   0    0.05];
ub = [  20    5     5    4.75];
    
[result resnorm resid] = lsqnonlin(@eravg_fiterr,x0,lb,ub,options);
      
results(rct).delay = result(1);
results(rct).amplitude = result(2);
results(rct).onset = result(3);
results(rct).dispersion = result(4);
results(rct).resnorm = resnorm;
results(rct).resid = resid;
    
hold on;

TRT = 1/1000; % Pretend TR is 1/1000 of actual for upsampling purposes
eventlength = EvLengthInSeconds/1*1000;

%Builds hrf function (subtraction of one gammma function from another)
p = [results.delay 16 results.dispersion 1 Inf results.onset];
hrf = spm_hrf(TRT,p);
hrf = hrf/(max(hrf)-min(hrf))*results.amplitude;
hrf = hrf;
hrf = hrf(1:eventlength)';

%Finds maximum, baseline, minimum, halfmax (for FWHM)
[MaxValue, MaxIndex] = max(hrf);
HalfMax = (MaxValue-0)/2 + 0;
OnsetMax = (MaxValue-0)/10 + 0;
[MinValue, MinIndex] = min(hrf);

for i = MaxIndex:-1:1
    HMPoint1 = i;
    if hrf(i) < HalfMax
        break
    end
end

for i = MaxIndex:eventlength
    HMPoint2 = i;
    if hrf(i) < HalfMax
        break
    end
end

for i = MaxIndex:-1:1
    RisePoint = i;
    if hrf(i) < OnsetMax
        break
    end
end

FWHM = HMPoint2 - HMPoint1;

t = 0:TRT:eventlength*TRT-TRT;
plot(t,hrf(1,:));

results.hrf = hrf;
results.fullhrf = fullhrf(end,:);
results.MaxValue = MaxValue;
results.MaxIndexInSeconds = MaxIndex * TRT;
results.MinValue = MinValue;
results.MinIndexInSeconds = MinIndex * TRT;
results.Baseline = 0;
results.MAmplitude = MaxValue - 0;
results.Range = MaxValue - MinValue;
results.FWHMInSeconds = FWHM * TRT;
results.OnsetTimeInSeconds = RisePoint * TRT;
results.HMPoint1InSeconds = HMPoint1 * TRT;
results.HMPoint2InSeconds = HMPoint2 * TRT;

Era = nanmean(events, 1);
PeakEra = max(Era);
TTPEra = find(abs(Era-PeakEra) < 0.001) * 1;

results.PeakEventAvg = PeakEra;
results.TimeToPeakEvAvg = TTPEra;

end