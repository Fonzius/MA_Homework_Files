function [S, magS, angleS, f, df] = myFFT(s, fs, varargin)
 
%This function computes fft of a signal
%If a pad value is provided, the number of computed samples is equal to
%padding value. Otherwise a power of two value is compute from actual
%signal length in order to find the best compromise between computational
%speed and frequency resolution.
% 2018 Riccardo De Lucia

%% INPUT ARGUMENTS
%s = the signal to be transformed
%fs = the signal sampling frequency
%paddingK = the desired number of samples for fft

%% OUTPUT ARGUMENTS
%S = the complex fft of signal s
%magS = the signal amplitude spectrum
%angleS = the signal phase spectrum
%f = the frequency bins array
%df = the fft frequency resolution

%% COMPUTATION
if(nargin <= 2)
   %find best power of two number of samples for fft
   paddingK = 2^(ceil(log2(length(s))));
else
   paddingK = varargin{1};
end

S = fft(s, paddingK); 

magS = abs(S);
angleS = angle(S);

df = fs/paddingK;

%compute the frequency bins array
f = (0:df:fs);
f = f(1:length(magS));
  
end
