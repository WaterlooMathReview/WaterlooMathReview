function [x,y]=damv2geom(bs,s)
%DAMV2GEOM	Gives geometry data for the damv2geom PDE model.
%
%   NE=DAMV2GEOM gives the number of boundary segments
%
%   D=DAMV2GEOM(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left-hand regions.
%   Row 4 contains the number of the right-hand regions.
%
%   [X,Y]=DAMV2GEOM(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

nbs=22;

if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 % start parameter value
  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 % end parameter value
  1 1 1 1 1 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 % left hand region
  0 0 0 0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 0 0 % right hand region
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
x(ii)=(-0.19999999999999996-(-0.24999999999999978))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(-0.24999999999999978);
y(ii)=(-0.59999999999999998-(-0.69999999999999996))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(-0.69999999999999996);
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=(-0.099999999999999881-(-0.19999999999999996))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(-0.19999999999999996);
y(ii)=(-0.4499999999999999-(-0.59999999999999998))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(-0.59999999999999998);
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=(2.2204460492503131e-016-(-0.099999999999999881))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(-0.099999999999999881);
y(ii)=(-0.34999999999999998-(-0.4499999999999999))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(-0.4499999999999999);
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=(0.14999999999999991-(2.2204460492503131e-016))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(2.2204460492503131e-016);
y(ii)=(-0.25-(-0.34999999999999998))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(-0.34999999999999998);
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=(0.5-(0.34999999999999987))*(s(ii)-d(1,5))/(d(2,5)-d(1,5))+(0.34999999999999987);
y(ii)=(-0.1000000000000001-(-0.15000000000000002))*(s(ii)-d(1,5))/(d(2,5)-d(1,5))+(-0.15000000000000002);
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=(-0.20000000000000021-(-0.25))*(s(ii)-d(1,6))/(d(2,6)-d(1,6))+(-0.25);
y(ii)=(0.59999999999999953-(0.69999999999999951))*(s(ii)-d(1,6))/(d(2,6)-d(1,6))+(0.69999999999999951);
end

% boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=(-0.1000000000000001-(-0.20000000000000021))*(s(ii)-d(1,7))/(d(2,7)-d(1,7))+(-0.20000000000000021);
y(ii)=(0.44999999999999962-(0.59999999999999953))*(s(ii)-d(1,7))/(d(2,7)-d(1,7))+(0.59999999999999953);
end

% boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=(0-(-0.1000000000000001))*(s(ii)-d(1,8))/(d(2,8)-d(1,8))+(-0.1000000000000001);
y(ii)=(0.34999999999999953-(0.44999999999999962))*(s(ii)-d(1,8))/(d(2,8)-d(1,8))+(0.44999999999999962);
end

% boundary segment 9
ii=find(bs==9);
if length(ii)
x(ii)=(0.14999999999999969-(0))*(s(ii)-d(1,9))/(d(2,9)-d(1,9))+(0);
y(ii)=(0.24999999999999956-(0.34999999999999953))*(s(ii)-d(1,9))/(d(2,9)-d(1,9))+(0.34999999999999953);
end

% boundary segment 10
ii=find(bs==10);
if length(ii)
x(ii)=(0.49999999999999978-(0.34999999999999964))*(s(ii)-d(1,10))/(d(2,10)-d(1,10))+(0.34999999999999964);
y(ii)=(0.099999999999999534-(0.14999999999999958))*(s(ii)-d(1,10))/(d(2,10)-d(1,10))+(0.14999999999999958);
end

% boundary segment 11
ii=find(bs==11);
if length(ii)
x(ii)=(-1.45-(-1.45))*(s(ii)-d(1,11))/(d(2,11)-d(1,11))+(-1.45);
y(ii)=(-0.69999999999999996-(0.69999999999999996))*(s(ii)-d(1,11))/(d(2,11)-d(1,11))+(0.69999999999999996);
end

% boundary segment 12
ii=find(bs==12);
if length(ii)
x(ii)=(1.2-(1.2))*(s(ii)-d(1,12))/(d(2,12)-d(1,12))+(1.2);
y(ii)=(0.69999999999999996-(-0.69999999999999996))*(s(ii)-d(1,12))/(d(2,12)-d(1,12))+(-0.69999999999999996);
end

% boundary segment 13
ii=find(bs==13);
if length(ii)
x(ii)=(0.25-(0.14999999999999991))*(s(ii)-d(1,13))/(d(2,13)-d(1,13))+(0.14999999999999991);
y(ii)=(-0.20000000000000007-(-0.25))*(s(ii)-d(1,13))/(d(2,13)-d(1,13))+(-0.25);
end

% boundary segment 14
ii=find(bs==14);
if length(ii)
x(ii)=(0.34999999999999987-(0.25))*(s(ii)-d(1,14))/(d(2,14)-d(1,14))+(0.25);
y(ii)=(-0.15000000000000002-(-0.20000000000000007))*(s(ii)-d(1,14))/(d(2,14)-d(1,14))+(-0.20000000000000007);
end

% boundary segment 15
ii=find(bs==15);
if length(ii)
x(ii)=(0.24999999999999978-(0.14999999999999969))*(s(ii)-d(1,15))/(d(2,15)-d(1,15))+(0.14999999999999969);
y(ii)=(0.19999999999999962-(0.24999999999999956))*(s(ii)-d(1,15))/(d(2,15)-d(1,15))+(0.24999999999999956);
end

% boundary segment 16
ii=find(bs==16);
if length(ii)
x(ii)=(0.34999999999999964-(0.24999999999999978))*(s(ii)-d(1,16))/(d(2,16)-d(1,16))+(0.24999999999999978);
y(ii)=(0.14999999999999958-(0.19999999999999962))*(s(ii)-d(1,16))/(d(2,16)-d(1,16))+(0.19999999999999962);
end

% boundary segment 17
ii=find(bs==17);
if length(ii)
x(ii)=(0.49999999999999978-(0.49999999999999978))*(s(ii)-d(1,17))/(d(2,17)-d(1,17))+(0.49999999999999978);
y(ii)=(-0.1000000000000001-(-0.69999999999999996))*(s(ii)-d(1,17))/(d(2,17)-d(1,17))+(-0.69999999999999996);
end

% boundary segment 18
ii=find(bs==18);
if length(ii)
x(ii)=(0.5-(0.5))*(s(ii)-d(1,18))/(d(2,18)-d(1,18))+(0.5);
y(ii)=(0.69999999999999951-(0.099999999999999534))*(s(ii)-d(1,18))/(d(2,18)-d(1,18))+(0.099999999999999534);
end

% boundary segment 19
ii=find(bs==19);
if length(ii)
x(ii)=(-0.25-(-1.45))*(s(ii)-d(1,19))/(d(2,19)-d(1,19))+(-1.45);
y(ii)=(0.69999999999999951-(0.69999999999999996))*(s(ii)-d(1,19))/(d(2,19)-d(1,19))+(0.69999999999999996);
end

% boundary segment 20
ii=find(bs==20);
if length(ii)
x(ii)=(1.2-(0.49999999999999978))*(s(ii)-d(1,20))/(d(2,20)-d(1,20))+(0.49999999999999978);
y(ii)=(0.69999999999999996-(0.69999999999999951))*(s(ii)-d(1,20))/(d(2,20)-d(1,20))+(0.69999999999999951);
end

% boundary segment 21
ii=find(bs==21);
if length(ii)
x(ii)=(-0.24999999999999978-(-1.45))*(s(ii)-d(1,21))/(d(2,21)-d(1,21))+(-1.45);
y(ii)=(-0.69999999999999996-(-0.69999999999999996))*(s(ii)-d(1,21))/(d(2,21)-d(1,21))+(-0.69999999999999996);
end

% boundary segment 22
ii=find(bs==22);
if length(ii)
x(ii)=(1.2-(0.5))*(s(ii)-d(1,22))/(d(2,22)-d(1,22))+(0.5);
y(ii)=(-0.69999999999999996-(-0.69999999999999996))*(s(ii)-d(1,22))/(d(2,22)-d(1,22))+(-0.69999999999999996);
end

end
