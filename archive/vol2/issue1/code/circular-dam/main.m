clear all;

%set-up the geometry, see man file for more specific information
% p is the point matrix and contains the verticies of each triangle
% e is the edge matrix
% t is the triangle matrix
[p,e,t]=initmesh('circgeom');

%refines the grid. Uncomment when ready
[p,e,t]=refinemesh('circgeom',p,e,t);
%[p,e,t]=refinemesh('circgeom',p,e,t);

%% if unable to run pdetools uncomment
%{

p=dlmread('p');
e=dlmread('e');
t=dlmread('t');

%}
%% Constants
%number of verticies
NumPoints=length(p);
%number of triangles by counting number of columns in t
NumTri=size(t,2);
%number of boundary points
NumBound=length(e);

%gravity
g=9.8;

%pretty sure this does nothing
%epsilon=0.1;
%set(axes,'ZLim',[0.5 1.5],'YLim',[-1,1],'XLim',[-1,1])

%the time step!
Deltat = 0.001;
totaltime=100;
%% Create the data matrix, TriData that contains the velocity information
% TriData(1) is the triangle label
% TriData(2) is h
% TriData(3) is uh
% Tridata(4) is vh


TriData=zeros(4,NumTri);
%needed for time step evolution
TriDataNext=zeros(4,NumTri);

% Create a matrix, TriInfo that contains info about triangle (ccw orient)
% the first index is the triangle index (can use this number to call
% labels)
% 2,3 are the x and y components of the n1 vector
% 4,5 are the x and y components of the n2 vector
% 6,7 are the x and y components of the n3 vector
% 8 contains boundary info, 0 for non, 1 for boundary
% 9,10,11 are the lengths of the n1,n2,n3 sides respectively 
% 12 is the area of the triangle
TriInfo=zeros(12,NumTri);

%compute some important info
for i=1:NumTri
    %compute the norms
    [a b c]=trianglenorm(p(:,t(1,i)),p(:,t(2,i)),p(:,t(3,i)));
    TriInfo(1,i)=i;
    TriInfo(2:3,i)=a(1:2);
    TriInfo(4:5,i)=b(1:2);
    TriInfo(6:7,i)=c(1:2);
    
    %compute the side lengths of each triangle as describe above
    [a b c]=tri_length(p(:,t(1,i)),p(:,t(2,i)),p(:,t(3,i)));
    TriInfo(9,i)=a;
    TriInfo(10,i)=b;
    TriInfo(11,i)=c;
    
    %now compute the area of the triangle using Heron's formula
    s=(a+b+c)/2;
    TriInfo(12,i)=sqrt(s*(s-a)*(s-b)*(s-c));
end

%% Find the Boundary and order labels properly
%
% This is entirely problem specific! The 5th entry of the edge matrix
% determines what edgesegment it is. IF we want to have interior initial
% conditions that are non-trivial we will draw them in pdetool. However
% these are not actually boundary conditions. Thus we must neglect them.
% The 5th entry tells us that!
%  
% For this problem the edges are 1,2,3,4 and the circle is 5,6,7,8. We DOT
% NOT care about the edges 5,6,7,8 thus we must not count them. Hence we
% add another condition that ensures that e(5,j)<= 4
%
%
%
%
% compute if it is a boundary point
% this works by simply looping through the edge matrix and checking if
% the three possible sides are on the boundary


%define sorted MATLAB boundary matrix 
EdgeSort=sort(e(1:2,:));
for j=1:NumTri
    for i=1:NumBound
        %it is an actual boundary!
        if( e(5,i)<5 )
            %check first side
            if( (t(1,j)==EdgeSort(1,i) && t(2,j)==EdgeSort(2,i)) || (t(2,j)==EdgeSort(1,i) && t(1,j)==EdgeSort(2,i)) )
                TriInfo(8,j)=1;
            end
            %check the second side
            if( (t(3,j)==EdgeSort(1,i) && t(2,j)==EdgeSort(2,i)) || (t(2,j)==EdgeSort(1,i) && t(3,j)==EdgeSort(2,i)) )
                TriInfo(8,j)=1;
            end
            %check the third side
            if( (t(1,j)==EdgeSort(1,i) && t(3,j)==EdgeSort(2,i)) || (t(3,j)==EdgeSort(1,i) && t(1,j)==EdgeSort(2,i)) )
                TriInfo(8,j)=1;
            end
        end
    end
end

% Finds the adjacent cells. The three values are the labels of the
% neighbouring cells. See m-file for description of search. 
%
% For example if EdgeMatrix(:,TriInfo(1,140)) = (255,144,229) tells us that
% the labels of the neighbouring triangles are 255,144,229 in no particular order. 
% If it returns zero than that triangle is on the boundary 
% 
EdgeMatrix=edgefind(NumTri,t);


%Now order those so that EdgeMatrix has the edges going in the ccw
%direction, i.e. EdgeMatrix(1,i) corresponds to the normal n1 of the
%triangle i, EdgeMatrix(2,i) to n2, EdgeMatrix(3,i) to n3
for i=1:NumTri
    %not a boundary
    if(TriInfo(8,i)~=1)
        PropOrder=order_triangles_nb(i,EdgeMatrix,t);
        for k=1:3
            EdgeMatrix(k,i)=PropOrder(k);
        end
    %sorts the boundary, the boundy will ALWAYS have a 0 for that edge
    else
        PropOrder=order_triangles_b(i,EdgeMatrix,t);
        for k=1:3
            EdgeMatrix(k,i)=PropOrder(k);
        end  

    end
end




%Now that booking keeping is done we can begin actual numerical simulation

%initial conditions: for this problem an initial circular dam
for i=1:NumTri
    TriData(1,i)=i;
    %circle radius
    rad1=0.4;
    z=sqrt(p(1,t(1:3,i)).^2+p(2,t(1:3,i)).^2);
    if( z(1) <= rad1 && z(2)<= rad1 && z(3)<=rad1)
        TriData(2,i)=1.9;
    else
        TriData(2,i)=1.0;
    end
    
    %initial velocities are zero
    TriData(3,i)=0;
    TriData(4,i)=0;
end


%needed for making animated gifs (uncomment as necessary)
%vidObj = VideoWriter('circular_dam_refined.avi');
%vidObj.Quality=75;
%open(vidObj);


%just check that the initial data looks good, comment out
pdeplot(p,e,t,'xydata',TriData(2,:),'zdata',TriData(2,:),'mesh','on')

for timestep=1:totaltime
%do the actual Riemann solving for each triangle

for i=1:NumTri
    %if the triangle is not on the boundary
    if(TriInfo(8,i)~=1)
        %we need the three fluxes through the other cells 

        %create a matrix that contains information about the three adajcent cells
        SideFluxes=TriData(:,EdgeMatrix(:,i));
	%compute the flux through each of the triangle edges using the above information
        Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
        %do the time-step and compute the next flux
        TriDataNext(1,i)=i;
        TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
    %it must be on the boundary

    %there is a more clever way to do this 
    else
     %the boundary is n1
     if(EdgeMatrix(1,i)==0)
	%the normal points in x
        if(abs(TriInfo(2,i))==1)
	    %set-up the flux matrix 
            SideFluxes=[[TriData(1,i);TriData(2,i);-TriData(3,i);TriData(4,i)],TriData(:,EdgeMatrix(2:3,i))];
	    %compute the flux matrix
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
	%the normal point in y
        else
            SideFluxes=[[TriData(1,i);TriData(2,i);TriData(3,i);-TriData(4,i)],TriData(:,EdgeMatrix(2:3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        end
     %the boundary is n2
     end
     if (EdgeMatrix(2,i)==0)
         if(abs(TriInfo(2,i))==1)
            SideFluxes=[TriData(:,EdgeMatrix(1,i)),[TriData(1,i);TriData(2,i);-TriData(3,i);TriData(4,i)],TriData(:,EdgeMatrix(3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        else
            SideFluxes=[TriData(:,EdgeMatrix(1,i)),[TriData(1,i);TriData(2,i);TriData(3,i);-TriData(4,i)],TriData(:,EdgeMatrix(3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
         end
     end
     %the boundary is n3
     if (EdgeMatrix(3,i)==0)
        if(abs(TriInfo(2,i))==1)
            SideFluxes=[TriData(:,EdgeMatrix(1:2,i)),[TriData(1,i);TriData(2,i);-TriData(3,i);TriData(4,i)]];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        else
             SideFluxes=[TriData(:,EdgeMatrix(1:2,i)),[TriData(1,i);TriData(2,i);TriData(3,i);-TriData(4,i)]];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        end
     end
    end
end
TriData(2:4,:)=TriDataNext(2:4,:);

pdeplot(p,e,t,'xydata',TriData(2,:),'zdata',TriData(2,:),'mesh','on','colorbar','off');
axis([-1 1 -1 1 0.9 2])
currFrame=getframe;
%writeVideo(vidObj,currFrame);
end 

%close(vidObj);

%{
pdeplot(p,e,t,'xydata',TriData(2,:),'zdata',TriData(2,:),'mesh','on')
M(k)=getframe;
%}
