function [TckResultCT, navSolutionsCT] = trackingCT_POS(file,signal,track,cmn, Acquired,cnslxyz,eph,sbf)
%Purpose:
%   Conventional tracking and positioning using EKF and WLS
%Inputs: 
%	file        - parameters related to the data file to be processed,a structure 
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris
%	sbf         - parameters used for pseudorange estimation
%
%Outputs:
%	TckResultCT         - conventional tracking results
%	navSolutionsCT   	- navigation solutions in conventional tracking mode
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
sv_clk              = zeros(1,32);
clkBias_kf         	= 0; 
clkBias_wls         	= 0; 
clkDrift            = 0;  
estusr_kf           = cnslxyz;%
estusr_wls           = cnslxyz;%
estVel              = zeros(1,3); 

Spacing = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing];

sv      = Acquired.sv;
f0      = signal.codeFreqBasis;
fL1     = signal.Fc;
fs      = signal.Fs ;
pdi     = track.pdi; 
t       = signal.ms;
svlength    = length(sv);
datalength  = track.msToProcessCT;

% Kalman Filter Parameter
num_state   = 8;

% error state vector
error_state = zeros(num_state,1);
total_state = [cnslxyz,zeros(1,5)]';%zeros(num_state,1);

% system transition matrix
Dynamic_Model = diag(zeros(1,num_state));
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model*pdi*t;

% error covariance matrix
state_cov = 1e5*diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]);

% process (System) noise noise covariance matrix
process_noise(1:3,1:3) = diag(ones(1,3)*2e-1);
process_noise(4:6,4:6) = diag(ones(1,3)*1e-1);
process_noise(7,7) = 1e-1;
process_noise(8,8) = 1e-2;

% measurement noise covariance matrix
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*3e-1;
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-1;

% parameters for measurement noise variance update
flag_corrCovEst2 = 1;
counterUptR = 0;
counter_r = 0;
thresUptR = 200/pdi;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(track.DLLBW,track.DLLDamp,track.DLLGain);
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);

% initialize tracking parameters using acquisition results
for svindex = 1:length(sv)
    prn                     = sv(svindex);
    codetemp                = generateCAcode(prn);
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
    AcqCodeDelay(svindex)   = Acquired.codedelay(svindex);
    file_ptr(svindex)       = signal.Sample - AcqCodeDelay(svindex) -1 + file.skip *fs*t;  
    carrFreq(svindex)       = Acquired.fineFreq(svindex);
    AcqFreq(svindex)        = Acquired.fineFreq(svindex);
    
    oldcodedelay_pos(svindex) = 0;
    oldabsoluteSample_pos(svindex) = 0;
end

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

eph_idx     = ones(1,svlength);

corrUpdateSec   = 0.1;
corrUpt         = corrUpdateSec/(pdi*t);
counter_corr    = corrUpt-1 * ones(svlength,1); 
        
% Tracking parameters
carrNco      = zeros(1,svlength);
oldCarrNco  = zeros(1,svlength);
oldCarrError       = zeros(1,svlength);
codeNco         = zeros(1,svlength);
code_outputLast     = zeros(1,svlength);
DLLdiscriLast       = zeros(1,svlength);
remChip             = zeros(1,svlength); 
codeFreq            = ones(1,svlength)*f0; 
remCarrPhase        = zeros(1,svlength);
carrError           = zeros(1,svlength);
codeError           = zeros(1,svlength);
delayValue          = zeros(svlength,datalength/pdi); 
 
    
h = waitbar(0,['Conventional Tracking, Length: ',num2str(datalength),' ms,', '  Please wait...']);
%%
for msIndex = 1: datalength/pdi 
    waitbar(msIndex/(datalength/pdi),h)
    for svindex = 1 :svlength
        prn = sv(svindex);
        
        % read raw data file
        numSample = round((signal.codelength*pdi-remChip(svindex))/(codeFreq(svindex)/signal.Fs));  
               
        delayValue(svindex,msIndex) = numSample - signal.Sample*pdi;
        
        fseek(file.fid, file_ptr(svindex)*file.dataPrecision*file.dataType,'bof');  
        if file.dataPrecision == 2
            rawsignal = fread(file.fid,numSample*file.dataType,'int16')';
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal = (sin_rawsignal - mean(sin_rawsignal)) + 1i.*(cos_rawsignal-mean(cos_rawsignal));
        else
            rawsignal = fread(file.fid,numSample*file.dataType,'int8')'; %
        end

        file_ptr(svindex)   = file_ptr(svindex) + numSample;  
        
        t_CodeEarly 	= (0 + Spacing(1) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * ...
                                codeFreq(svindex)/signal.Fs + Spacing(1) + remChip(svindex));
        t_CodePrompt 	= (0 + Spacing(2) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * ...
                                codeFreq(svindex)/signal.Fs + Spacing(2) + remChip(svindex));
        t_CodeLate   	= (0 + Spacing(3) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * ...
                                codeFreq(svindex)/signal.Fs + Spacing(3) + remChip(svindex));
                            
        CodeEarly 	= Code(svindex,(ceil(t_CodeEarly) + 1));
        CodePrompt	= Code(svindex,(ceil(t_CodePrompt) + 1));
        CodeLate	= Code(svindex,(ceil(t_CodeLate) + 1));
        remChip(svindex) = t_CodePrompt(numSample) + codeFreq(svindex)/signal.Fs - signal.codelength*pdi;
        
        CarrTime = (0:numSample)./signal.Fs;
        Wave = 2*pi*((carrFreq(svindex)).*CarrTime) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample+1), 2*pi);
        carrsig = exp(1i.* Wave(1:numSample));
        InphaseSignal    = imag(rawsignal .* carrsig);
        QuadratureSignal = real(rawsignal .* carrsig);
        
        E_i = sum(CodeEarly    .*InphaseSignal);  
        E_q = sum(CodeEarly    .*QuadratureSignal);
        P_i	= sum(CodePrompt   .*InphaseSignal);  
        P_q = sum(CodePrompt   .*QuadratureSignal);
        L_i	= sum(CodeLate     .*InphaseSignal);  
        L_q = sum(CodeLate     .*QuadratureSignal);
        
        % Calculate CN0
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                meanZk  = mean(Zk(svindex,:));
                varZk   = var(Zk(svindex,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_CT(snrIndex(svindex),svindex) =  abs(10*log10(1/(1*t*pdi) * NA2/(2*varIQ)));
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        % Implement code loop filter and generate NCO command
        E = sqrt(E_i^2+E_q^2);
        L = sqrt(L_i^2+L_q^2);
        codeError(svindex) = 0.5*(E-L)/(E+L);  % DLL discriminator
        codeNco(svindex) = code_outputLast(svindex) + (tau2code/tau1code)*(codeError(svindex)...
                                    - DLLdiscriLast(svindex)) + codeError(svindex)* ((pdi*t)/tau1code);
        DLLdiscriLast(svindex) = codeError(svindex);
        code_outputLast(svindex) = codeNco(svindex);
        codeFreq(svindex) = signal.codeFreqBasis - codeNco(svindex);
        
        % PLL discriminator
        carrError(svindex) = atan(P_q/P_i)/(2*pi);  % PLL discriminator
        carrNco(svindex) = oldCarrNco(svindex) + (tau2carr/tau1carr)*(carrError(svindex) ...
                                        - oldCarrError(svindex)) + carrError(svindex) * ((pdi*t)/tau1carr);
        oldCarrNco(svindex) = carrNco(svindex); 
        oldCarrError(svindex) = carrError(svindex);
        carrFreq(svindex)  = AcqFreq(svindex) + carrNco(svindex);  % Modify carrier freq 
        
        % Data Recording
        TckResultCT(prn).P_i(msIndex)             = P_i;
        TckResultCT(prn).P_q(msIndex)             = P_q;
        TckResultCT(prn).E_i(msIndex)             = E_i;
        TckResultCT(prn).E_q(msIndex)             = E_q;
        TckResultCT(prn).L_i(msIndex)             = L_i;
        TckResultCT(prn).L_q(msIndex)             = L_q;
        TckResultCT(prn).carrError(msIndex)       = carrError(svindex);
        TckResultCT(prn).codeError(msIndex)       = codeError(svindex);
        TckResultCT(prn).codeFreq(msIndex)        = codeFreq(svindex);
        TckResultCT(prn).carrFreq(msIndex)        = carrFreq(svindex);
        TckResultCT(prn).numSample(msIndex)       = numSample;
        TckResultCT(prn).remChip(msIndex)         = remChip(svindex);        
        TckResultCT(prn).remCarrPhase(msIndex)    = remCarrPhase(svindex);
        TckResultCT(prn).absoluteSample(msIndex)  = ftell(file.fid); 
        TckResultCT(prn).absoluteSampleCodedelay(msIndex)  = mod(TckResultCT(prn).absoluteSample(msIndex)/...
                                                                (file.dataPrecision*file.dataType),fs*t);
        TckResultCT(prn).codedelay(msIndex)       = mod(TckResultCT(prn).absoluteSample(msIndex)/ ...
                                                        (file.dataPrecision*file.dataType),fs*t);
        TckResultCT(prn).delayValue(msIndex)      = delayValue(svindex,msIndex);
        
    end % end for svindex in Tracking
    
	%% Interpolate
    for svindex=1:svlength
    	prn = sv(svindex);
        x1 = [oldabsoluteSample_pos(svindex)/(file.dataPrecision*file.dataType) ...
                TckResultCT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType)];
      	y2 = [oldcodedelay_pos(svindex) TckResultCT(prn).codedelay(msIndex)];
      	x_s = file.skip*fs*t + msIndex*pdi*fs*t ;%+ (sbf.nav1(prn))*fs*t 
        codedelay_pos(msIndex, svindex) = interp1(x1,y2,x_s);
    end
    oldcodedelay_pos = codedelay_pos(msIndex, :);
    oldabsoluteSample_pos(svindex) = TckResultCT(prn).absoluteSample(msIndex);
    
   
    %% Pseudorange calculation and correction     
    [pseudorange, ~] = pr_est(Acquired, sbf, signal,TckResultCT, msIndex);
 
    usr_clk = clkBias_kf;
    estusr = estusr_kf;
    
    for svindex = 1 : svlength
        prn = sv(svindex);
        
        %
        time = eph(prn).TOW1(eph_idx(svindex)) - sbf.sfb1(prn) * 1/50 ...
              	- sbf.nav1(prn) * 1/1000 + msIndex*pdi*1/1000;
            
        % difference for all satellites in unit of samples        
        fs_new = (f0/codeFreq(svindex))*fs;
        diff_of_dat_pos(svindex) = sbf.sfb1(prn)*20*fs*t ...     
                                        + sbf.nav1(prn)*fs*t ...
                                        + (codedelay_pos(msIndex,svindex));                                    
        % time of transmition                          
        tot_est_pos(svindex) = eph(prn).TOW1(eph_idx(svindex)) ...   % TOW1: TOW of current subframe 1
                                    - diff_of_dat_pos(svindex)/fs ...% 
                                    + pdi*msIndex*t ...
                                    - (1/cmn.cSpeed)*(sv_clk(prn)); 

        [svxyz(svindex,:), sv_vel(svindex,:), sv_clk(prn), sv_clk_vel(prn), grpdel] = ...
                                        svPosVel(prn,eph,tot_est_pos(svindex),eph_idx(svindex));
        
        % C/A-code pseudorange corrected for satellite clock and Tgd 
        prvec(svindex)      = pseudorange(svindex) - usr_clk + sv_clk(prn) - grpdel*cmn.cSpeed ;  
        
        % Adjust satellite position coordinates for earth rotation correction
        svxyzr(svindex,:)   = erotcorr(svxyz(svindex,:),prvec(svindex));   
        svenu               = xyz2enu(svxyzr(svindex,:),estusr(1:3));
        
        % tropospheric and ionospheric delay correction	
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) ==  corrUpt
            svenu           = xyz2enu(svxyzr(svindex,:),estusr(1:3));%
            el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
            az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
            az(svindex)     = az_rad(svindex)*180/pi;
            el(svindex)     = el_rad(svindex)*180/pi;
            temp            = xyz2llh(estusr(1:3));
            user_ll         = [temp(1:2).*180/pi temp(3)];
            ionodel(svindex)        = ionocorr(tot_est_pos(svindex),svxyzr(svindex,:), estusr(1:3));
            tropodel_unb3(svindex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svindex))); 
            counter_corr(svindex)   = 0;     
        end
        
        prvec(svindex) = prvec(svindex) - ionodel(svindex) - tropodel_unb3(svindex); 
    end 

    
    %% Position cal using LS method

    % Position cal using KF
    for svindex = 1 : svlength
        r                       = sqrt(sum((svxyzr(svindex,:) - estusr_kf(1:3)).^2));
        a(svindex,:)            = (svxyzr(svindex,:)-estusr_kf(1:3))/r;
        H(svindex,:)            = [-a(svindex,:) 0 0 0 1 0];
        H(svindex+svlength,:)   = [0 0 0 -a(svindex,:) 0 1];
        
        d_p(svindex,1)          = (TckResultCT(Acquired.sv(svindex)).carrFreq(msIndex) ...
                                        - signal.IF)*cmn.cSpeed/fL1 + sum( sv_vel(svindex,:) .* a(svindex,:));
 
        % measured pseudorange rate
        prr_measured(svindex)	= (TckResultCT(Acquired.sv(svindex)).carrFreq(msIndex) - signal.IF)*cmn.cSpeed/fL1 ;
        
        %predicted pseudorange rate 
        prr_predicted(svindex)	= (estVel-sv_vel(svindex,:))*a(svindex,:)';
       
        d_p(svindex,1) = prr_predicted(svindex) - prr_measured(svindex) - clkDrift + sv_clk_vel(prn);% 
        pr_delta(svindex,1) = prvec(svindex) - norm(svxyzr(svindex,:) - estusr_kf(1:3));% ; 
    end
    
    Z = [pr_delta; d_p]; % measurements
    	
	error_state = zeros(num_state,1);
    error_state = Transistion_Matrix * error_state; % predict
    state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise;
    
    kalman_gain = state_cov * transpose(H) * inv(H * state_cov * transpose(H) +  mesurement_noise);
    
    counterUptR = counterUptR + 1; % counter for update measurement noise variance, R
    recordR(counterUptR,:) = (Z - H*error_state);% measurement innovation
    
    error_state = error_state + kalman_gain*(Z-H*error_state);
    state_cov = (eye(num_state) - kalman_gain * H) * state_cov ;
 
    total_state  = total_state + error_state;
    estusr_kf   = total_state(1:3)';
    estVel      = total_state(4:6)';
    clkBias_kf  = total_state(7);
    clkDrift    = total_state(8);
       
    %record results
 	llh     = xyz2llh(estusr_kf);  
    L_b     = llh(1);
    lamda_b = llh(2);
    C_e_n = [ -sin(lamda_b)           cos(lamda_b)         	 0;...
              -sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)	 cos(L_b);...
              -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)	-sin(L_b);];  
    usrenuvel(msIndex,:) =  C_e_n * estVel';
    
    usrenu(msIndex,:)                   = xyz2enu(estusr_kf(1:3),cnslxyz);
    usrllh(msIndex,:)                   = xyz2llh(estusr_kf(1:3));
    usrllh(msIndex,1:2)                 = usrllh(msIndex,1:2)*180/pi;
    
    navSolutionsKF.usrTime(msIndex,:)        = time;  
    navSolutionsKF.usrPos(msIndex,:)      = estusr_kf;
    navSolutionsKF.usrVel(msIndex,:)      = estVel;
    navSolutionsKF.usrPosENU(msIndex,:)   = usrenu(msIndex,:);
    navSolutionsKF.usrVelENU(msIndex,:)   = usrenuvel(msIndex,:);
    navSolutionsKF.usrPosLLH(msIndex,:)   = usrllh(msIndex,:);
    navSolutionsKF.error_state(msIndex,:)       = error_state;
    navSolutionsKF.clkDrift(msIndex)       = clkDrift;
    navSolutionsKF.clkBias(msIndex)       = clkBias_kf;
    navSolutionsKF.newZ(msIndex,:)        = Z;
    navSolutionsKF.meas_inno(msIndex,:)   = (Z - H * error_state);
    navSolutionsKF.predicted_z(msIndex,:) = H * error_state;
    navSolutionsKF.kalman_gain(:,:,msIndex)	 = kalman_gain;
    navSolutionsKF.state_cov(msIndex,:)      = diag(state_cov);
    navSolutionsKF.satEA(msIndex,:)      = el;
    navSolutionsKF.satAZ(msIndex,:)      = az;
   
    % predict postion and clkBias at next epoch   
	total_state  = Transistion_Matrix * total_state;
    estusr_kf = total_state(1:3)';
    clkBias_kf = total_state(7)';
 
    % adaptively update measurement noise covariance R
    if flag_corrCovEst2 == 1 && counterUptR == thresUptR
        tmpR = diag(var(recordR));
        mesurement_noise(1:svlength,1:svlength)=tmpR(1:svlength,1:svlength) ;% ;  
        mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength)= ...
                                tmpR(svlength+1:2*svlength,svlength+1:2*svlength) ;% 
        
        for idx = 1 : svlength
            if mesurement_noise(idx,idx) >= 12000
                mesurement_noise(idx,idx) = 12000;
            elseif mesurement_noise(idx,idx) <= 0.01 	
                mesurement_noise(idx,idx) = 0.01;	
            end
            if mesurement_noise(idx+svlength,idx+svlength) >= 400
                mesurement_noise(idx+svlength,idx+svlength) = 400;
            elseif mesurement_noise(idx+svlength,idx+svlength) <= 0.01 
                mesurement_noise(idx+svlength,idx+svlength) = 0.01;	
            end
        end
        counterUptR = 0;
        counter_r = counter_r + 1;
        navSolutionsKF.R(counter_r,:) = diag(mesurement_noise);
    end

    % 
    TOW_USR_CT(msIndex) = time;
    fprintf('CT_KF: index = %4d TOW = %f E = %f N = %f U = %f  B = %f D = %f\n',msIndex, ...
        TOW_USR_CT(msIndex),usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3), clkBias_kf,clkDrift); 
 
    navSolutionsCT = navSolutionsKF;
    
end % end for msIndex
    
close(h);

save(['navSolCT_',file.fileName], 'navSolutionsCT','eph','TOW_USR_CT');
save(['tckRstCT_',file.fileName], 'TckResultCT','CN0_CT'); 

