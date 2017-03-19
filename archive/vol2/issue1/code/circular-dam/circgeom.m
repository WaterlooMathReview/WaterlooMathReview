function [x,y]=circgeom(bs,s)
%CIRCGEOM	Gives geometry data for the circgeom PDE model.
%
%   NE=CIRCGEOM gives the number of boundary segments
%
%   D=CIRCGEOM(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left-hand regions.
%   Row 4 contains the number of the right-hand regions.
%
%   [X,Y]=CIRCGEOM(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

nbs=8;

if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 % start parameter value
  1 1 1 1 1 1 1 1 % end parameter value
  0 0 0 0 2 2 2 2 % left hand region
  1 1 1 1 1 1 1 1 % right hand region
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
x(ii)=(1-(1))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(1);
y(ii)=(-1-(1))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(1);
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=(-1-(1))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(1);
y(ii)=(-1-(-1))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(-1);
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=(-1-(-1))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(-1);
y(ii)=(1-(-1))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(-1);
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=(1-(-1))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(-1);
y(ii)=(1-(1))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(1);
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=0.39999999999999991*cos(1.5707963267948966*s(ii)+(-3.1415926535897931))+(0);
y(ii)=0.39999999999999991*sin(1.5707963267948966*s(ii)+(-3.1415926535897931))+(0);
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=0.39999999999999991*cos(1.5707963267948966*s(ii)+(-1.5707963267948966))+(0);
y(ii)=0.39999999999999991*sin(1.5707963267948966*s(ii)+(-1.5707963267948966))+(0);
end

% boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=0.39999999999999991*cos(1.5707963267948966*s(ii)+(0))+(0);
y(ii)=0.39999999999999991*sin(1.5707963267948966*s(ii)+(0))+(0);
end

% boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=0.39999999999999991*cos(1.5707963267948966*s(ii)+(-4.7123889803846897))+(0);
y(ii)=0.39999999999999991*sin(1.5707963267948966*s(ii)+(-4.7123889803846897))+(0);
end

end
