function [ signalData ] = SFilt( signal, time )
% SFilt
%   Savitzky-Golay smoothing filter and RMS extraction

signalFiltered = sgolayfilt(signal,7,23);
signalNoise = signal - signalFiltered;
signalRMS = sqrt(1/max(time) * trapz(time,signalFiltered.^2));
signalRMSNoise = sqrt(1/max(time) * trapz(time,signalNoise.^2));

signalData.signal = signalFiltered;
signalData.noise = signalNoise;
signalData.amplitudeRMS = signalRMS;
signalData.RMSNoise = signalRMSNoise;
end