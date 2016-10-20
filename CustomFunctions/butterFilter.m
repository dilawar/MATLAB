function onlyButterworth = butterFilter(myData, samplingRate, lowPassCutoffFreq, highPassCutoffFreq, filterOrder)

nyquistsamplingRate = samplingRate/2;

%Butterworth
%low pass
lowWn =lowPassCutoffFreq/nyquistsamplingRate;
[b,a]=butter(filterOrder, lowWn, 'low'); %weirdly sets up a bandstop.. NOTE: we have adjusted this in Wn
y=filtfilt(b,a,myData);
%high pass
highWn =highPassCutoffFreq/nyquistsamplingRate;
[b,a]=butter(filterOrder, highWn, 'high'); %weirdly sets up a bandstop.. NOTE: we have adjusted this in Wn
onlyButterworth=filtfilt(b,a,y);

end