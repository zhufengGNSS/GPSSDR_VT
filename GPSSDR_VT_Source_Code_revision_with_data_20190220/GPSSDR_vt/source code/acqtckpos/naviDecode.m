function [ephemeris ALLTckResult for_prest] = naviDecode(Acquired, ALLTckResult)
%Purpose
%   decode navigation data from inphase prompt (P_i) channel.
%   *MUST be used after tracking process to get ALLTckResult
%   *Navigation Data length MUST be longer than 30 seconds
%Inputs:
%	Acquired        - acquisition results
%	ALLTckResult  	- conventional tracking results
%Outputs:
%	ephemeris       - ephemeris obtained from the navigation bit stream
%	ALLTckResult	- revised conventional tracking results 
%	for_prest       - parameters used for pseudorange estimation including
%                       the first navigation data (nav1) point and the
%                       begining of subframe 1 (sfb1)
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
% Setting the paramethers of C/N0 estimation (PRM)
M = 20; K =5; T = 1e-3;

% Setting the parameter of navigation data decode
preamble = [-1 1 1 1 -1 1 -1 -1];

flag_sf     = zeros(length(Acquired.sv),5);
flag_eph    = zeros(length(Acquired.sv),1);

sv          = Acquired.sv;
[ephemeris] = ini_eph(Acquired);

for svindex = 1 : length(Acquired.sv)
    case1_index = 0;
    flag_sfb1   = 0;
    prn         = sv(svindex);
    ephemeris(prn).updateflag = 0;
    fprintf('\n');
    RawNavigationData                           = ALLTckResult(prn).P_i(:);
    NaviDatams(find(RawNavigationData>=0))      = 1;
    NaviDatams(find(RawNavigationData<0))       = -1;
    NaviDataXORms(find(RawNavigationData>=0))   = 0;
    NaviDataXORms(find(RawNavigationData<0))    = 1;
    msdatalength                                = length(NaviDatams);
    
    for startms = 2:msdatalength
        % check = abs(sum(NaviData(startms:startms+19)))/20
        if NaviDatams(startms) ~= NaviDatams(startms-1)
            break;
        end
    end
    
    % startms = mod(startms,20);, if startms == 0, startms = 1;, end;
    % numpadd = 21 - startms;
    ALLTckResult(prn).nav1  = startms-1;
    for_prest.nav1(prn)     = startms-1;
    
    % Lock Detectors (using C/N0 to determine whether the signal is successful tracked or not)
    WBP = zeros(1,K);
    index_k = 0;
    NBP_I = 0;
    NBP_Q = 0;
    idx = 0;
    idx2 = 0;
    for (index = startms : msdatalength)
        idx = idx + 1;
        index_k = mod(ceil(idx / M),K); if (index_k)==0,index_k = K;,end
        WBP(index_k) = WBP(index_k) + (ALLTckResult(prn).P_i(index)^2+ALLTckResult(prn).P_q(index)^2);
        %         fprintf('[%3d] [%3d] %2d %f %f\n',idx,index,NaviData(index),ALLTckResult(prn).P_i(index),ALLTckResult(prn).P_q(index))
        NBP_I = NBP_I + ALLTckResult(prn).P_i(index);
        NBP_Q = NBP_Q + ALLTckResult(prn).P_q(index);
        %         fprintf('[%3d] [%d] %2d %f %f %f \n',idx,index_k,NaviData(index),WBP(index_k),NBP_I,NBP_Q)
        if mod(idx,M) == 0
            NBP(index_k) = (NBP_I^2 + NBP_Q^2);
            NP(index_k) = NBP(index_k) / WBP(index_k);
            NBP_I = 0; NBP_Q =0;
        end
        if mod(idx,(M*K)) == 0
            idx2 = idx2 + 1;
            u = 1/K * sum(NP);
            ALLTckResult(prn).CN0(idx2) = 10*log10(1/T * (u-1)/(M-u));
            WBP = zeros(1,K);
        end
    end
    fprintf('[%2d] CN0 = %f\n',svindex,mean(ALLTckResult(prn).CN0));
    
    % Bit Synchronization (to check which ms is the start of a navigation bit)
    % change NaviData from ms to 50bps
    idx = 0; 
    idx2 = 0;
    for index = startms:msdatalength
        idx     = idx + 1;
        temp    = floor(idx/(M*K))+1;
        if (msdatalength-index) > 100 % 100 because C/N0
            if ALLTckResult(prn).CN0(floor(idx/(M*K))+1) >= 1 % 30dB-Hz by Parkinson
                if mod(idx,20) == 0
                    idx2                = idx2 + 1;
                    NaviDataXOR(idx2)   = NaviDataXORms(index);
                    NaviData(idx2)      = NaviDatams(index);
                end
            else
                break;
            end
        else
            break;
        end
    end
    fprintf('sv %2d successful tracked bit = %5d startms = %2d\n',prn,idx2,startms);
    
    % Matching the Preamble of the subframe (data demodulation)
    navidatalength  = length(NaviDataXOR);
    flag            = 0; % if all parity check ok, flag = 1;
    for index = 8 : navidatalength
        if ALLTckResult(prn).CN0(floor((index-7)/K)+1) >= 30 && flag==0
            if length(index : navidatalength) > 360 % 300 next subframe
                if abs(sum((NaviData(index-7:index).*preamble)))>7.99 && abs(sum((NaviData(index-7+300:index+300).*preamble)))>7.99 % find preamble
                    temp        = NaviData(index-7:end);
                    end_HOW     = sum(temp(60:-1:59)) ;
                    end_HOW2    = sum(temp(360:-1:359));                  
                    if end_HOW ~= 0 && end_HOW2~=0
                        % parity check
                        if flag == 0
                            [pass NaviDataXOR] = paritychk_James(NaviDataXOR,index-7);
                            if ( pass == 1)
                                flag = 1;
                                fprintf('parity successful svindex = %2d index =%4d\n',svindex,index)
                            end
                        end % end flag = 0
                        
                        if flag == 1
                            num_sf = floor(length(NaviDataXOR(index-7:end))/300);
                            for idx = 1:num_sf
                                subframe(idx,:)     = NaviDataXOR(index-7+300*(idx-1):index-7+299+300*(idx-1));
                                ephemeris(prn).TOW  = [ephemeris(prn).TOW (bin2dec_GPSSDR(subframe(idx,47:-1:31))-1)* 6]; % -1 because of z-count2TOW transformation
                                subframeid          = bin2dec_GPSSDR(subframe(idx,52:-1:50));
                                fprintf('[%2d] - %4d TOW = %6d subframe = %d end_HOW = %2d\n',svindex, index, ephemeris(prn).TOW(end), subframeid,end_HOW);
                                
                                switch subframeid
                                    case 1
                                        case1_index                         = case1_index + 1;
                                        ALLTckResult(prn).sfb1(case1_index) = (index-7)+(idx-1)*300;
                                        if flag_sfb1 == 0
                                            for_prest.sfb1(prn) = (index-7)+(idx-1)*300;
                                            fprintf('================ [%2d] - %4d ================\n',prn, for_prest.sfb1(prn));
                                            flag_sfb1 =1;
                                        end
                                        
                                        ephemeris(prn).sfb1         = [ephemeris(prn).sfb1 (index-7)+(idx-1)*300];
                                        ephemeris(prn).weeknum      = [ephemeris(prn).weeknum bin2dec_GPSSDR(subframe(idx,70:-1:61)) + 1024];
                                        ephemeris(prn).TOW1         = [ephemeris(prn).TOW1 (bin2dec_GPSSDR(subframe(idx,47:-1:31)) -1 )* 6];
                                        
                                        ephemeris(prn).N = [ephemeris(prn).N bin2dec_GPSSDR(subframe(idx,76:-1:73))];
                                        if (ephemeris(prn).N(end) <= 6), accuracy = 2^(1+ephemeris(prn).N(end)/2);,
                                        elseif (ephemeris(prn).N(end) > 6), accuracy = 2^(ephemeris(prn).N(end)-2);,
                                        elseif (ephemeris(prn).N(end)==15),disp('satellite in risk (N=15)!\n');,end
                                        
                                        ephemeris(prn).health = [ephemeris(prn).health bin2dec_GPSSDR(subframe(idx,82:-1:78))];
                                        if (ephemeris(prn).health(end)~=0),fprintf('SV in bad health %d\n',ephemeris(prn).health(end));,end;
                                        
                                        ephemeris(prn).IODC = [ephemeris(prn).IODC bin2dec_GPSSDR(subframe(idx,218:-1:211))];                                         
                                        ephemeris(prn).TGD = [ephemeris(prn).TGD comp2dec(subframe(idx,204:-1:197),-31)]; % comp for complement
                                        ephemeris(prn).toc = [ephemeris(prn).toc bin2dec_GPSSDR(subframe(idx,234:-1:219)) * 2^4];
                                        if (ephemeris(prn).toc(end)) > 604784, fprintf('Fault toc not in effective range. toc = %d\n',ephemeris(prn).toc(end));,end;
                                        
                                        ephemeris(prn).af2 = [ephemeris(prn).af2 comp2dec(subframe(idx,248:-1:241),-55)]; 
                                        ephemeris(prn).af1 = [ephemeris(prn).af1 comp2dec(subframe(idx,264:-1:249),-43)];
                                        ephemeris(prn).af0 = [ephemeris(prn).af0 comp2dec(subframe(idx,292:-1:271),-31)];
                                        flag_sf(svindex,1) = 1;
                                        
                                    case 2
                                        ephemeris(prn).IODE2 = [ephemeris(prn).IODE2 bin2dec_GPSSDR(subframe(idx,68:-1:61))];
                                        %                             IODE2 = bin2dec_GPSSDR(subframe(idx,68:-1:61))
                                        ephemeris(prn).Crs = [ephemeris(prn).Crs comp2dec(subframe(idx,84:-1:69),-5)]; 
                                        ephemeris(prn).deltan = [ephemeris(prn).deltan comp2dec(subframe(idx,106:-1:91),-43) * pi];
                                        ephemeris(prn).M0 = [ephemeris(prn).M0 comp2dec([subframe(idx,144:-1:121) subframe(idx,114:-1:107)],-31) * pi];
                                        ephemeris(prn).Cuc = [ephemeris(prn).Cuc comp2dec(subframe(idx,166:-1:151),-29)];
                                        ephemeris(prn).ecc = [ephemeris(prn).ecc bin2dec_GPSSDR([subframe(idx,204:-1:181) subframe(idx,174:-1:167)]) * 2^(-33)];
                                        ephemeris(prn).Cus = [ephemeris(prn).Cus comp2dec(subframe(idx,226:-1:211),-29)];
                                        %                             ephemeris(prn).Cus
                                        %                             subframe(211:226)
                                        %                             comp2dec(subframe(226:-1:211),-29)
                                        ephemeris(prn).sqrta = [ephemeris(prn).sqrta bin2dec_GPSSDR([subframe(idx,264:-1:241) subframe(idx,234:-1:227)]) * 2^(-19)];
                                        ephemeris(prn).toe = [ephemeris(prn).toe bin2dec_GPSSDR(subframe(idx,286:-1:271)) * 2^4];
                                        if (ephemeris(prn).toe(end) > 604784), fprintf('Fault toe not in effective range. toe = %d\n',ephemeris(prn).toe(end));,end;
                                        flag_sf(svindex,2) = 1;
                                        
                                    case 3
                                        ephemeris(prn).Cic = [ephemeris(prn).Cic comp2dec(subframe(idx,76:-1:61),-29)];
                                        ephemeris(prn).omegae = [ephemeris(prn).omegae comp2dec([subframe(idx,114:-1:91) subframe(idx,84:-1:77)],-31) * pi];
                                        ephemeris(prn).Cis = [ephemeris(prn).Cis comp2dec(subframe(idx,136:-1:121),-29)];
                                        ephemeris(prn).i0 = [ephemeris(prn).i0 comp2dec([subframe(idx,174:-1:151) subframe(idx,144:-1:137)],-31) * pi];
                                        ephemeris(prn).Crc = [ephemeris(prn).Crc comp2dec(subframe(idx,196:-1:181),-5)];
                                        ephemeris(prn).w = [ephemeris(prn).w comp2dec([subframe(idx,234:-1:211) subframe(idx,204:-1:197)],-31) * pi];
                                        ephemeris(prn).omegadot = [ephemeris(prn).omegadot comp2dec(subframe(idx,264:-1:241),-43) * pi];
                                        ephemeris(prn).IODE3 = [ephemeris(prn).IODE3 bin2dec_GPSSDR(subframe(idx,278:-1:271))];
                                        %                             IODE3 = bin2dec_GPSSDR(subframe(idx,278:-1:271))
                                        ephemeris(prn).idot = [ephemeris(prn).idot comp2dec(subframe(idx,292:-1:279),-43) * pi];
                                        flag_sf(svindex,3) = 1;
                                        
                                    case 4
                                        flag_sf(svindex,4) = 1;
                                    case 5
                                        flag_sf(svindex,5) = 1;
                                        
                                end % end switch
                                
                                
                                
                                if flag_sf(svindex,1)==1 & flag_sf(svindex,2)==1 & flag_sf(svindex,3)==1 & flag_sf(svindex,4)==1 & flag_sf(svindex,5)==1 ...
                                        & ephemeris(prn).health(end)==0
                                    fprintf('SV%2d all Ephemeris update!\n',prn);
                                    % flag_eph(svindex) = 1;
                                    
                                    ephemeris(prn).updateflag = 1;
                                    ephemeris(prn).updatetime = [ephemeris(prn).updatetime (index + (idx)*300) * 20 + (startms-1)]; % unit:  ms
                                    ephemeris(prn).updatetime_tow = [ephemeris(prn).updatetime_tow ephemeris(prn).TOW(end) + 6] ; % unit: second (+6 because of it's next frame)
                                    % if flag_sf(svindex,1)==1
                                    % end
                                    flag_sf(svindex,1)=0;
                                    flag_sf(svindex,2)=0;
                                    flag_sf(svindex,3)=0;
                                    flag_sf(svindex,4)=0;
                                    flag_sf(svindex,5)=0;
                                end
                                
                                
                            end % end for num_sf
                        end % end flag = 1
                    end % end_HOW
                end % find preamble
            end % end if length
        end % end if C/N0 > 30
    end % end index of NaviDataXOR
end % end for svindex
 