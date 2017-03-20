%% Calculates the flux through a surface
%
% This function computes the flux through the three surfaces 
% add more here 
%
%
% current_tri_info is a 4 vector that is the centre of the triangle
% current_tri_info(1) is the label
% current_tri_info(2,3,4) are h,u,v respectively
% in main this is given by TriData(:,i) where i is the tri considered
%
%
% adj_tri_info is a 4x3 matrix that is the same as above except with labels
% 1,2,3 to distinguish the three sides
% 
%
%%

function interface = riemann_solver(current_tri_info,T, adj_tri_info)

%gravity
g=9.81;
%convenient definition
c = sqrt(g*current_tri_info(2));
%index of the triangle 
tri_index =current_tri_info(1);

%calculate the flux through each of the three surfaces
for i = 1:3
    
    %grab the normal vectors through the surface
    n_x = T(2*i, tri_index);
	n_y = T(2*i+1, tri_index);
    
    %compute the Roe average state, (Q^{+}+Q^{-})/2
	h_avg = (current_tri_info(2) + adj_tri_info(2,i))/2;
    %display(h_avg);
    u_avg = (current_tri_info(3) + adj_tri_info(3,i))/2;
    %display(u_avg);
    v_avg = (current_tri_info(4) + adj_tri_info(4,i))/2;
    %display(v_avg);
	
    %we had that the left and right states are
    h_minus = current_tri_info(2);
    u_minus= current_tri_info(3);
    v_minus = current_tri_info(4);
    h_plus= adj_tri_info(2,i);
    u_plus= adj_tri_info(3,i);
    v_plus= adj_tri_info(4,i);
    q_plus=[h_plus, u_plus, v_plus];
    q_minus=[h_minus,u_minus,v_minus];
    
    %Compute the Roe Average matrix 
    %compute the right eigenvector matrix, eq 8 of A&C
	R = [0, 1, 1; n_y, u_avg - c*n_x, u_avg + c*n_x; -n_x, v_avg - c*n_y, v_avg + c*n_y];
    %compute the left eigenvector matrix, eq 9 of A&C
	L = [ -(u_avg*n_y - v_avg*n_x), n_y, -n_x; (u_avg*n_y + v_avg*n_x)/(2*c) + 1/2, -n_x/(2*c), -n_y/(2*c); -(u_avg*n_x + v_avg*n_y)/(2*c) + 1/2, n_x/(2*c), n_y/(2*c)];
    %compute the eigenvalue matrix, eq 7of A&C
	Lambda = [abs(u_avg*n_x + v_avg*n_y), 0, 0; 0, abs(u_avg*n_x + v_avg*n_y - c), 0; 0, 0, abs(u_avg*n_x + v_avg*n_y + c)];
	%F = J^xn_x + J^yn_y
	fplus = jx(h_plus,u_plus,v_plus)*n_x + jy(h_plus,u_plus,v_plus)*n_y;
    
	fminus = jx(h_minus,u_minus,v_minus).*n_x + jy(h_minus,u_minus,v_minus).*n_y;
    
    %
	interface(i,:) = (1/2)*(fplus' + fminus') - 0.5*R*Lambda*L*(q_plus'-q_minus');
%%end

end