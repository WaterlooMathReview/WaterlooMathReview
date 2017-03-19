%% Computes the three outward surface normals 
%
% Input the triangle and output the surface normals for that vetor
% 
%%
function [n1,n2,n3] = trianglenorm(a,b,c)

%the corners of the triangle
x1=a(1,1);
x2=b(1,1);
x3=c(1,1);

y1=a(2,1);
y2=b(2,1);
y3=c(2,1);

%compute the normals
[n1x,n1y]=normal(x1,y1,x2,y2);
[n2x,n2y]=normal(x2,y2,x3,y3);
[n3x,n3y]=normal(x3,y3,x1,y1);

%return the normals
n1=[n1x,n1y];
n2=[n2x,n2y];
n3=[n3x,n3y];
end
