function [estusr, dop] = olspos(prvec,svxyzmat,initpos,tol)
%OLSPOS Compute position from satellite positions and pseudoranges
%	via ordinary least squares.
%
%estusr = OLSPOS(prvec,svxyzmat,initpos,tol)
%
%INPUTS
%prvec = vector of 'measured' pseudoranges for satellites
%specified in svxyzmat
%svxyzmat(i,1:3) = position of satellite i in user defined
%                   cartesian coordinates.
%initpos = optional argument. Initial 'estimate' of user state:
%           three-dimensional position and clock offset
%           (in user defined coordinates). Used to speed up
%           iterative solution. Initial clock offset is
%           optional with default value = 0.
%tol = optional argument. Tolerance value used to determine
%           convergence of iterative solution. Default value = 1e-3
%
%	OUTPUTS
%estusr(1:3) = estimated user x, y, and z coordinates
%estusr(4) = estimated user clock offset
%	Note: all four elements of estusr are in the same units
%           as those used in prvec
 
%M. & S. Braasch 11-96
%Copyright (c) 1996 by GPSoft
%All Rights Reserved.
%
if nargin<4,tol=1e-3;end
if nargin<3,initpos=[0 0 0 0];end
if nargin<2,error('insufficient number of input arguments'),end
[m,n]=size(initpos);
if m>n, estusr=initpos';else,estusr=initpos;end
if max(size(estusr))<3 
    error('must define at least 3 dimensions in INITPOS')
end
if max(size(estusr))<4,estusr=[estusr 0];end
numvis=max(size(svxyzmat));
beta=[1e9 1e9 1e9 1e9];
maxiter=10;
iter=0;
% while ((iter<maxiter)&&(norm(beta)>tol)),
while ( (norm(beta)>tol)) 
    for N = 1:numvis 
        pr0 = norm(svxyzmat(N,:)-estusr(1:3));
        y(N,1) = prvec(N) - pr0 - estusr(4);
    end 
    H = hmat(svxyzmat,estusr(1:3));
    beta = H\y;
    estusr=estusr+beta';
    iter=iter+1;
end

% Calculate DOP
dop = zeros(1,4);
Q = inv(H'*H);
dop(1)  = sqrt(trace(Q));                       % GDOP
dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
dop(4)  = sqrt(Q(3,3));                         % VDOP


