function [pass NaviDataXOR] = paritychk_James(NaviDataXOR,idx_sf1)
%Purpose
%   Vector tracking and positioning
%Inputs: 
%	NaviDataXOR    	- navigation data 
%	idx_sf1         - 
%Outputs:
%	pass            - 
%	NaviDataXOR   	- 
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
pass = 1;
% from xor to multiplication (0->1 1->-1)
NaviDataXOR(find(NaviDataXOR==1)) = -1; 
NaviDataXOR(find(NaviDataXOR==0)) =  1; 

datalength = floor(length(NaviDataXOR(idx_sf1:end))/30) * 30;

%     The Parity Matrix
%     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
H = [ 1  1  1  0  1  1  0  0  0  1  1  1  1  1  0  0  1  1  0  1  0  0  1  0;...
      0  1  1  1  0  1  1  0  0  0  1  1  1  1  1  0  0  1  1  0  1  0  0  1;...
      1  0  1  1  1  0  1  1  0  0  0  1  1  1  1  1  0  0  1  1  0  1  0  0;...
      0  1  0  1  1  1  0  1  1  0  0  0  1  1  1  1  1  0  0  1  1  0  1  0;...
      1  0  1  0  1  1  1  0  1  1  0  0  0  1  1  1  1  1  0  0  1  1  0  1;...
      0  0  1  0  1  1  0  1  1  1  1  0  1  0  1  0  0  0  1  0  0  1  1  1];

for idx = idx_sf1 : 30 : datalength
    D30star = NaviDataXOR(idx-1);
    D29star = NaviDataXOR(idx-2);
    NaviDataXOR(idx:idx+23) = D30star * NaviDataXOR(idx:idx+23);
    d = NaviDataXOR(idx:idx+23);
    Df = [D29star D30star D29star D30star D30star D29star];
    for idx2 = 1 : 6
        temp = H(idx2,:) .* d;
        p(idx2) = prod([Df(idx2) temp(find(temp~=0))]);
    end    
    if (p ~= NaviDataXOR(idx+24:idx+29))
        pass = 0;
        disp('\n\nParity Check Failed!!\n\n');
%         fprintf('idx_sf1 = %4d idx = %4d\n',idx_sf1,idx);
%         disp(p)
%         disp(NaviDataXOR(idx+24:idx+29))    
    end 
end
% pass
NaviDataXOR = -NaviDataXOR;
NaviDataXOR = (NaviDataXOR+1)./2;

end % function end