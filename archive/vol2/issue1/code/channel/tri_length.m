%% Calculates the length of a segment
%
% Inputs are the x_1,y_1,x_2,y_2 co-ordinates
%
%
%%

function [l1 l2 l3] = tri_length(a,b,c)

%the corners of the triangle
x1=a(1,1);
x2=b(1,1);
x3=c(1,1);

y1=a(2,1);
y2=b(2,1);
y3=c(2,1);

%compute the normals
l1=sqrt(abs(x1-x2)^2+abs(y1-y2)^2);
l2=sqrt(abs(x2-x3)^2+abs(y2-y3)^2);
l3=sqrt(abs(x3-x1)^2+abs(y3-y1)^2);
end