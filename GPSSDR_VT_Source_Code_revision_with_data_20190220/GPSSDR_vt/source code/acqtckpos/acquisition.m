function  Acquired = acquisition(file,signal,acq)
%Purpose:
%   Perform signal acquisition using parallel code phase search algorithm
%   Fine frequency is found using long FFT.
%Inputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal      - parameters related to signals,a structure
%Outputs:
%	Acquired	- acquisition results, i.e. space vehicle (SV),signal-to-noise 
%                   ratio (SNR), Doppler, code delay, and fine carrier
%                   frequency
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% Initialization 
Acquired.sv = []; 
Acquired.SNR = []; 
Acquired.Doppler = []; 
Acquired.codedelay = []; 
Acquired.fineFreq = [];
sampleindex = 1: signal.Sample;

%% read data
fseek(file.fid,file.skip*signal.Sample*file.dataPrecision*file.dataType,'bof');
if file.dataPrecision == 2  % 16bits
    rawsignal = fread(file.fid,signal.Sample*file.dataType*4,'int16')'; % 4ms data for acquisition
    sin_rawsignal = rawsignal(1:2:length(rawsignal));
    cos_rawsignal = rawsignal(2:2:length(rawsignal));
    rawsignal = sin_rawsignal - mean(sin_rawsignal) + 1i*(cos_rawsignal-mean(cos_rawsignal));
else
    rawsignal = fread(file.fid,signal.Sample*file.dataType*4,'int8')'; % 4ms data for acquisition
end

%% acquisition process 
for freqband = 1 : acq.freqNum
    dopplershift = acq.freqMin + acq.freqStep*(freqband-1);
    carrier(freqband,:) = exp(1i*2*pi*(signal.IF + dopplershift) * sampleindex ./ signal.Fs); 
end

for svindex = 1:32
    svindex
    ocode = generateCAcode(svindex);
    ocode = [ocode ocode];
    scode = ocode(ceil(sampleindex.*(signal.codeFreqBasis/signal.Fs)));
    correlation = zeros(acq.freqNum,signal.Sample);
    for idx = 1 : 4  % 4ms data
        for freqband = 1 : acq.freqNum        
        	replica = scode;
            temp1 = rawsignal(1+(idx-1)*signal.Sample:idx*signal.Sample).* carrier(freqband,:);
            temp2 = conj(fft(temp1));
            temp3 = fft(replica);
            correlation(freqband,:) = correlation(freqband,:) + abs(ifft(temp3.*temp2)).^2;
        end
    end
    [peak, fbin] = max(max(abs(correlation')));
    [peak, codePhase] = max(max(abs(correlation)));
    Doppler = acq.freqMin + acq.freqStep * (fbin-1);
    
    codechipshift = ceil(signal.Fs/signal.codeFreqBasis);
    SNR = 10 * log10((peak^2)/(sum(correlation(fbin,[1:codePhase-codechipshift codePhase+codechipshift:end]).^2)...
                /length(correlation(fbin,[1:codePhase-codechipshift codePhase+codechipshift:end]))));

    if SNR >= 18  % acquisition thredhold
        Acquired.sv        = [Acquired.sv        svindex];
        Acquired.SNR       = [Acquired.SNR       SNR];
        Acquired.Doppler   = [Acquired.Doppler   Doppler];
        Acquired.codedelay = [Acquired.codedelay codePhase-1];
    end
end


%% fine frequency calculation 
acq.L = 10;
fseek(file.fid,file.skip*signal.Sample*file.dataPrecision*file.dataType,'bof');
if file.dataPrecision == 2
    rawsignal = fread(file.fid,signal.Sample*file.dataType*(acq.L+1),'int16')';   
    sin_rawsignal = rawsignal(1:2:length(rawsignal));
    cos_rawsignal = rawsignal(2:2:length(rawsignal));
    longrawsignal = sin_rawsignal - mean(sin_rawsignal) + 1i*(cos_rawsignal-mean(cos_rawsignal));
else
    longrawsignal = fread(file.fid,signal.Sample*file.dataType*(acq.L+1),'int8')';  
end

for svindex = 1 : length(Acquired.sv)     
    caCode = generateCAcode(Acquired.sv(svindex));    
    codeValueIndex = floor((1/signal.Fs*(1:acq.L*signal.Sample))/(1/signal.codeFreqBasis));    
    longCaCode = caCode((rem(codeValueIndex, signal.codelength) + 1));
    CarrSignal = longrawsignal((signal.Sample-Acquired.codedelay(svindex)):...
                    (signal.Sample-Acquired.codedelay(svindex))+acq.L*signal.Sample - 1).* longCaCode;
     
    fftlength = length(CarrSignal) * 20;
    fftSignal = abs(fft(CarrSignal, fftlength));
    halffftlength = ceil((fftlength)/2);
    [~, FreqPeakIndex] = max(fftSignal(1:halffftlength));
    fineDoppler = FreqPeakIndex * (signal.Fs/fftlength);    
    
    if file.dataType == 2
        if Doppler>= 0
            fineDoppler = FreqPeakIndex * (signal.Fs/fftlength);
        else
            fineDoppler = -FreqPeakIndex * (signal.Fs/fftlength);
        end
    end    
    
    Acquired.fineFreq = [Acquired.fineFreq  fineDoppler];
        
    fprintf('SV[%2d] SNR = %2.2f Code phase = %5d  Doppler frequency = %5f\n', ...
            Acquired.sv(svindex),Acquired.SNR(svindex),Acquired.codedelay(svindex), ...
            Acquired.fineFreq(svindex)-signal.IF);
end

end % end function