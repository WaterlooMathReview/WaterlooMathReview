%% Computes the outwards norm of an edge
%
%  Recieves as input the four co-ordinates of the line
%  and returns vector normal to the surface
%
%%
function [nx,ny] = normal(x1,y1,x2,y2)

x=x2-x1;
y=y2-y1;
N=sqrt(x^2+y^2);
nx=y/N;
ny=-x/N;
end



