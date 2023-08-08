function [] = plotFFT_linearFreqScale(magS, angleS, f, df, fs, maxFreq, h)
  %this function plots the magnitude and phase spectrum of the provided
  %signal in frequency linear scale
  % 2018 Riccardo De Lucia

  %% INPUT PARAMETERS
  %magS = the signal magnitude spectrum
  %angleS = the signal phase spectrum
  %f = the frequency bins array
  %df = the fft frequency resolution
  %fs = the sampling frequency
  %maxFreq = the maximum frequency to be plotted (maximum value = fs/2)
  %h = the figure handle for the plot
   
  %% COMPUTATION 
  
  %set the maximum frequency to be plotted
  if(maxFreq>0 && maxFreq<fs/2)
      freq = maxFreq;
  else
      freq = fs/2;
  end
  
  %trim accordingly the arrays
  maxBin = ceil(freq/df);
  magS = magS(1:maxBin);
  angleS = angleS(1:maxBin);
  f = f(1:maxBin);
  
  %set the current figure
  figure(h);
  
  %get axes handle for future reference to specific subplots
  ax1 = subplot(2, 1, 1);
  ax2 = subplot(2, 1, 2);
  
  %allowing plots overlapping
  hold on
  %select subplot by setting current axes
  axes(ax1)
  plot(f(1:end), magS(1:end));
  title('Magnitude spectrum');
  ylabel('Magnitude');
  xlabel('Frequency');
  leg = sprintf('df=%.2f Hz', df);
  legend(leg);
  grid on;
  axis([1 f(end) min(magS) max(magS)+0.1*max(magS)]);
  hold off
  
  %allowing plots overlapping
  hold on
  %select subplot by setting current axes
  axes(ax2);
  plot(f(1:end), angleS(1:end));
  title('Phase spectrum');
  ylabel('Phase');
  xlabel('Frequency');
  leg = sprintf('df=%.2f Hz', df);
  legend(leg);
  grid on;
  axis([1 f(end) min(angleS) max(angleS)]);
  hold off
