function jx_flux = jx(h, u, v)

g=9.81;
jx_flux = [h*u, (1/2)*g*h*h + h*u*u, h*u*v];
end

