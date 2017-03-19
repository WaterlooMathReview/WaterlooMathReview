function jy_flux = jy(h, u, v)

g=9.81;
jy_flux=[h*v,h*u*v,(1/2)*g*h*h + h*v*v];
end 