function h = hmat(svmat,usrpos)
%Purpose:
%   Find direction cosine matrix,H
%Inputs:
%   svmat    - satellite position in ECEF coordinate
%   usrpos   - user position in ECEF coordinate
%Outputs:
%       h    - direction cosine matrix
%-------------------------------------------------------------------------- 

%% 
N = max(size(svmat));
if N < 4,
    error('insufficient number of satellites')
else
    tmppos = usrpos;
    [m,n] = size(tmppos);
    if m > n, tmppos = tmppos'; end,
    h = ones(N,4);
    % linear expansion in user position, the precision depends on the user coordinates
    for i = 1:N,
        tmpvec = tmppos - svmat(i,:);
        h(i,1:3) = tmpvec./norm(tmpvec);
    end,
end,
