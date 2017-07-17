function [ P, f ] = fourierTransform( q, t )
% Fourier transform of signal
%
% [ P, f ] = fourierTransform( q, t )
%
% Takes vectors time t (in s) and quantity q in input, returns vectors
% absolute value P (normalized to 1) and frequency f (in Hz). If no t is entered
% the samplig period is assumed to be 1 s.
%

if nargin > 1
    % calculate sampling period
    dt = mean(gradient(t));
else
    dt = 1;
end

% remove f = 0 Hz component
q = q - mean(q);

%% old fft, gives mirror-image result
% calculate fft
% Q = fft(q);
% P = Q.*conj(Q);
% P = P/max(P);
% f = 1/dt/length(Q)*(0:(length(Q)-1));

%% new fft
% calculate fft
L = 2^(nextpow2(length(q)) + 3);
Q = fft(q,L);
P2 = abs(Q)/L;
P1 = P2(1:L/2);
P = P1/max(P1);
f = 1/dt*(0:(L/2 - 1))/L;

end

