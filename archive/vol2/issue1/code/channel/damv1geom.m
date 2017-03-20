function [x,y]=damv1geom(bs,s)
%DAMV1GEOM	Gives geometry data for the damv1geom PDE model.
%
%   NE=DAMV1GEOM gives the number of boundary segments
%
%   D=DAMV1GEOM(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left-hand regions.
%   Row 4 contains the number of the right-hand regions.
%
%   [X,Y]=DAMV1GEOM(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

nbs=17;

if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 % start parameter value
  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 % end parameter value
  1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 % left hand region
  0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 % right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error('Non-existent boundary segment number')
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 & n==1,
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) | n~=size(s,2),
  error('bs must be scalar or of same size as s');
end

if ~isempty(s),

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=(0.59999999999999987-(0.64999999999999991))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(0.64999999999999991);
y(ii)=(0.34999999999999987-(0.14999999999999991))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(0.14999999999999991);
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=(0.54999999999999982-(0.59999999999999987))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(0.59999999999999987);
y(ii)=(0.49999999999999989-(0.34999999999999987))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(0.34999999999999987);
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=(0.44999999999999984-(0.54999999999999982))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(0.54999999999999982);
y(ii)=(0.59999999999999998-(0.49999999999999989))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(0.49999999999999989);
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=(0.44999999999999984-(0.24999999999999981))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(0.24999999999999981);
y(ii)=(-0.34999999999999998-(-0.4499999999999999))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(-0.4499999999999999);
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=(0.54999999999999982-(0.44999999999999984))*(s(ii)-d(1,5))/(d(2,5)-d(1,5))+(0.44999999999999984);
y(ii)=(-0.25-(-0.34999999999999998))*(s(ii)-d(1,5))/(d(2,5)-d(1,5))+(-0.34999999999999998);
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=(0.59999999999999987-(0.54999999999999982))*(s(ii)-d(1,6))/(d(2,6)-d(1,6))+(0.54999999999999982);
y(ii)=(-0.10000000000000009-(-0.25))*(s(ii)-d(1,6))/(d(2,6)-d(1,6))+(-0.25);
end

% boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=(0.64999999999999991-(0.59999999999999987))*(s(ii)-d(1,7))/(d(2,7)-d(1,7))+(0.59999999999999987);
y(ii)=(0.099999999999999964-(-0.10000000000000009))*(s(ii)-d(1,7))/(d(2,7)-d(1,7))+(-0.10000000000000009);
end

% boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=(-1.3999999999999999-(-1.3999999999999999))*(s(ii)-d(1,8))/(d(2,8)-d(1,8))+(-1.3999999999999999);
y(ii)=(-0.45000000000000007-(0.69999999999999984))*(s(ii)-d(1,8))/(d(2,8)-d(1,8))+(0.69999999999999984);
end

% boundary segment 9
ii=find(bs==9);
if length(ii)
x(ii)=(1.25-(1.25))*(s(ii)-d(1,9))/(d(2,9)-d(1,9))+(1.25);
y(ii)=(0.69999999999999984-(-0.45000000000000007))*(s(ii)-d(1,9))/(d(2,9)-d(1,9))+(-0.45000000000000007);
end

% boundary segment 10
ii=find(bs==10);
if length(ii)
x(ii)=(0.34999999999999976-(0.24999999999999986))*(s(ii)-d(1,10))/(d(2,10)-d(1,10))+(0.24999999999999986);
y(ii)=(0.64999999999999991-(0.69999999999999984))*(s(ii)-d(1,10))/(d(2,10)-d(1,10))+(0.69999999999999984);
end

% boundary segment 11
ii=find(bs==11);
if length(ii)
x(ii)=(0.44999999999999984-(0.34999999999999976))*(s(ii)-d(1,11))/(d(2,11)-d(1,11))+(0.34999999999999976);
y(ii)=(0.59999999999999998-(0.64999999999999991))*(s(ii)-d(1,11))/(d(2,11)-d(1,11))+(0.64999999999999991);
end

% boundary segment 12
ii=find(bs==12);
if length(ii)
x(ii)=(0.64999999999999991-(0.6499999999999998))*(s(ii)-d(1,12))/(d(2,12)-d(1,12))+(0.6499999999999998);
y(ii)=(0.099999999999999964-(-0.4499999999999999))*(s(ii)-d(1,12))/(d(2,12)-d(1,12))+(-0.4499999999999999);
end

% boundary segment 13
ii=find(bs==13);
if length(ii)
x(ii)=(0.64999999999999991-(0.64999999999999991))*(s(ii)-d(1,13))/(d(2,13)-d(1,13))+(0.64999999999999991);
y(ii)=(0.69999999999999984-(0.14999999999999991))*(s(ii)-d(1,13))/(d(2,13)-d(1,13))+(0.14999999999999991);
end

% boundary segment 14
ii=find(bs==14);
if length(ii)
x(ii)=(0.24999999999999986-(-1.3999999999999999))*(s(ii)-d(1,14))/(d(2,14)-d(1,14))+(-1.3999999999999999);
y(ii)=(0.69999999999999984-(0.69999999999999984))*(s(ii)-d(1,14))/(d(2,14)-d(1,14))+(0.69999999999999984);
end

% boundary segment 15
ii=find(bs==15);
if length(ii)
x(ii)=(1.25-(0.6499999999999998))*(s(ii)-d(1,15))/(d(2,15)-d(1,15))+(0.6499999999999998);
y(ii)=(0.69999999999999984-(0.69999999999999984))*(s(ii)-d(1,15))/(d(2,15)-d(1,15))+(0.69999999999999984);
end

% boundary segment 16
ii=find(bs==16);
if length(ii)
x(ii)=(0.24999999999999981-(-1.3999999999999999))*(s(ii)-d(1,16))/(d(2,16)-d(1,16))+(-1.3999999999999999);
y(ii)=(-0.4499999999999999-(-0.45000000000000007))*(s(ii)-d(1,16))/(d(2,16)-d(1,16))+(-0.45000000000000007);
end

% boundary segment 17
ii=find(bs==17);
if length(ii)
x(ii)=(1.25-(0.64999999999999991))*(s(ii)-d(1,17))/(d(2,17)-d(1,17))+(0.64999999999999991);
y(ii)=(-0.45000000000000007-(-0.4499999999999999))*(s(ii)-d(1,17))/(d(2,17)-d(1,17))+(-0.4499999999999999);
end

end
