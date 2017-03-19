clear all;

%set-up the geometr, see man file for info
[p,e,t]=initmesh('damv2geom');

%improve the grid. Uncomment when ready
[p,e,t]=refinemesh('damv2geom',p,e,t);
%[p,e,t]=refinemesh('damv2geom',p,e,t);

%%% if not available uncomment the following
%{
p=dlmread('p');
e=dlmread('e');
t=dlmread('t');
%}
%% Constants
%number of points
NumPoints=length(p);
%number of triangles
NumTri=size(t,2);
%number of boundary points
NumBound=length(e);

%gravity
g=9.8;
epsilon=0.1;
%set(axes,'ZLim',[0.5 1.5],'YLim',[-1,1],'XLim',[-1,1])

%the time step!
Deltat = 0.001;

%% Create the data matrix, TriData that contains all the physical info 
% the first index is the point label
% TriData(2) is h
% TriData(3) is uh
% Tridata(4) is vh

TriData=zeros(4,NumTri);
TriDataNext=zeros(4,NumTri);
BottomData=zeros(4,NumTri);
BottomDataNext=zeros(4,NumTri);

%% Create a matrix, TriInfo that contains info about triangle (ccw orient)
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
    
    %compute the side lengths of each cell as describe above
    [a b c]=tri_length(p(:,t(1,i)),p(:,t(2,i)),p(:,t(3,i)));
    TriInfo(9,i)=a;
    TriInfo(10,i)=b;
    TriInfo(11,i)=c;
    
    %now compute the area of the triangle using Heron's formula
    s=(a+b+c)/2;
    TriInfo(12,i)=sqrt(s*(s-a)*(s-b)*(s-c));
end

%% Find the Boundary
%
% This is entirely problem specific! The 5th entry of the edge matrix
% determines what edgesegment it is. IF we want to have interior initial
% conditions that are non-trivial we will draw them in pdetool. However
% these are not actually boundary conditions. Thus we must neglect them.
% The 5th entry tells us that!
%  
% For this problem the edges are 1-16 but they are all exterior. HOWEVER we
% do have that the following boundary conditions
% 19,20,21,22 are the river
% 11 is the ocean
% 12 is the channel boundary
% 1-10,13-18 is the gate 
%
% compute if it is a boundary point
% this works by simply looping through the edge matrix and checking if
% the three possible sides are on the boundary

%define sorted MATLAB boundary matrix
EdgeSort=sort(e(1:2,:));
for j=1:NumTri
    for i=1:NumBound
            %check first side
            if( (t(1,j)==EdgeSort(1,i) && t(2,j)==EdgeSort(2,i)) || (t(2,j)==EdgeSort(1,i) && t(1,j)==EdgeSort(2,i)) )
                TriInfo(8,j)=e(5,i);
            end
            %check the second side
            if( (t(3,j)==EdgeSort(1,i) && t(2,j)==EdgeSort(2,i)) || (t(2,j)==EdgeSort(1,i) && t(3,j)==EdgeSort(2,i)) )
                TriInfo(8,j)=e(5,i);
            end
            %check the third side
            if( (t(1,j)==EdgeSort(1,i) && t(3,j)==EdgeSort(2,i)) || (t(3,j)==EdgeSort(1,i) && t(1,j)==EdgeSort(2,i)) )
                TriInfo(8,j)=e(5,i);
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
    if(TriInfo(8,i)==0)
        PropOrder=order_triangles_nb(i,EdgeMatrix,t);
        for k=1:3
            EdgeMatrix(k,i)=PropOrder(k);
        end
    %sorts the boundary, the boundary will ALWAYS have a 0 for that edge
    else
        PropOrder=order_triangles_b(i,EdgeMatrix,t);
        for k=1:3
            EdgeMatrix(k,i)=PropOrder(k);
        end  

    end
end




%% Now we are ready to begin

%initial conditions: everything is at rest initially
for i=1:NumTri
    TriData(1,i)=i;
    TriData(2,i)=1;
    TriData(3,i)=0;
    TriData(4,i)=0;
end


%
%vidObj = VideoWriter('water_heigher.avi');
%vidObj.Quality=75;
%open(vidObj);
%just check that the initial data looks good, comment out
pdeplot(p,e,t,'xydata',TriData(2,:),'zdata',TriData(2,:),'mesh','on')

for timestep=1:400
%do the actual Riemann solving
for i=1:NumTri
    if(timestep<=15)
        ocean_u=5;
        ocean_v=0;
        ocean_h=1.4;
        ocean_vector=[(sqrt(ocean_h)-(TriData(3,i)-ocean_u)/sqrt(g))^2;ocean_u;ocean_v];
    else
        ocean_u=0;
        ocean_v=0;
        ocean_h=1;
        ocean_vector=[(sqrt(ocean_h)-(TriData(3,i)-ocean_u)/sqrt(g))^2;ocean_u;ocean_v];
    end
    %if it is not a boundary point
    if(TriInfo(8,i)==0)
        %we need the three fluxes through the other cells 
        %do something
        %create a matrix that contains the fluxes of the three other sides\
        SideFluxes=TriData(:,EdgeMatrix(:,i));
        Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
        %compute the next flux
        TriDataNext(1,i)=i;
        TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
    %it must be on the boundary
    else
     %the boundary is n1
     if(EdgeMatrix(1,i)==0)
        %the surface is on the river bank, i.e. TriInfo(8,i) = 20,21,16,17
        if((TriInfo(8,i)==22) || (TriInfo(8,i)==21) || (TriInfo(8,i)==20) || (TriInfo(8,i)==19))
            SideFluxes=[[TriData(1,i);TriData(2,i);TriData(3,i);-TriData(4,i)],TriData(:,EdgeMatrix(2:3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        %the surface is on the gate    
        elseif ( TriInfo(8,i)==(1||2||3||4||5||6||7||8||9||10||13||14||15||16||17||18))
            %the angle (CHANGES FOR each boundary)
            theta=atan2(TriInfo(3,i),TriInfo(2,i));
            SideFluxes=[[TriData(1,i);TriData(2,i);-TriData(3,i)*cos(2*theta)-TriData(4,i)*sin(2*theta);TriData(4,i)*cos(2*theta)+TriData(3,i)*sin(2*theta)],TriData(:,EdgeMatrix(2:3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        %the point is on the open sea boundary    
        elseif (TriInfo(8,i)==11)
            %ocean boundary conditions
            SideFluxes=[[TriData(1,i);ocean_vector],TriData(:,EdgeMatrix(2:3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);            
        %it is the non-ocean boundary
        else

            SideFluxes=[[TriData(1,i);TriData(2,i);TriData(3,i);TriData(4,i)],TriData(:,EdgeMatrix(2:3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);            
           
        end
     %the boundary is n2
     end
     if (EdgeMatrix(2,i)==0)
        %the surface is on the river bank, i.e. TriInfo(8,i) = 14,15,16,17
        if((TriInfo(8,i)==22) || (TriInfo(8,i)==21) || (TriInfo(8,i)==20) || (TriInfo(8,i)==19))
            SideFluxes=[TriData(:,EdgeMatrix(1,i)),[TriData(1,i);TriData(2,i);TriData(3,i);-TriData(4,i)],TriData(:,EdgeMatrix(3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        %the surface is on the gate    
        elseif ( TriInfo(8,i)==(1||2||3||4||5||6||7||8||9||10||13||14||15||16||17||18))
            %the angle (CHANGES FOR each boundary)
            theta=atan2(TriInfo(5,i),TriInfo(4,i));
            SideFluxes=[TriData(:,EdgeMatrix(1,i)),[TriData(1,i);TriData(2,i);-TriData(3,i)*cos(2*theta)-TriData(4,i)*sin(2*theta);TriData(4,i)*cos(2*theta)+TriData(3,i)*sin(2*theta)],TriData(:,EdgeMatrix(3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);            
        %the point is on the open sea boundary    
        elseif (TriInfo(8,i)==11)
            %ocean boundary conditions
            SideFluxes=[TriData(:,EdgeMatrix(1,i)),[TriData(1,i);ocean_vector],TriData(:,EdgeMatrix(3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);            
        else
            
          
            SideFluxes=[TriData(:,EdgeMatrix(1,i)),[TriData(1,i);TriData(2,i);TriData(3,i);TriData(4,i)],TriData(:,EdgeMatrix(3,i))];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);            
          
            
        end
     end
     %the boundary is n3
     if (EdgeMatrix(3,i)==0)
        %the surface is on the river bank, i.e. TriInfo(8,i) = 14,15,16,17
        if((TriInfo(8,i)==22) || (TriInfo(8,i)==21) || (TriInfo(8,i)==20) || (TriInfo(8,i)==19))
            SideFluxes=[TriData(:,EdgeMatrix(1:2,i)),[TriData(1,i);TriData(2,i);TriData(3,i);-TriData(4,i)]];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);
        %the surface is on the gate    
        elseif ( TriInfo(8,i)==(1||2||3||4||5||6||7||8||9||10||13||14||15||16||17||18))
            %the angle (CHANGES FOR each boundary)
            theta=atan2(TriInfo(7,i),TriInfo(6,i));
            SideFluxes=[TriData(:,EdgeMatrix(1:2,i)),[TriData(1,i);TriData(2,i);-TriData(3,i)*cos(2*theta)-TriData(4,i)*sin(2*theta);TriData(4,i)*cos(2*theta)+TriData(3,i)*sin(2*theta)]];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);            
        %the point is on the open sea boundary    
        elseif (TriInfo(8,i)==11)
            %ocean boundary conditions
            SideFluxes=[TriData(:,EdgeMatrix(1:2,i)),[TriData(1,i);ocean_vector]];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);            
        else
            
            SideFluxes=[TriData(:,EdgeMatrix(1:2,i)),[TriData(1,i);TriData(2,i);TriData(3,i);TriData(4,i)]];
            Flux=riemann_solver(TriData(:,i),TriInfo,SideFluxes);
            TriDataNext(1,i)=i;
            TriDataNext(2:4,i)=TriData(2:4,i)'-Deltat*(Flux(1,:)*TriInfo(9,i)+Flux(2,:)*TriInfo(10,i)+Flux(3,:)*TriInfo(11,i))/TriInfo(12,i);  
        end
     end
    end
end
TriData(2:4,:)=TriDataNext(2:4,:);
pdeplot(p,e,t,'xydata',TriData(2,:),'zdata',TriData(2,:),'mesh','on','colorbar','off');
axis([-1 1 -1 1 0.0 2])
view([140 34]);
%view([-180 0]);
currFrame=getframe;
%writeVideo(vidObj,currFrame);
end 

%close(vidObj);

%{
pdeplot(p,e,t,'xydata',TriData(2,:),'zdata',TriData(2,:),'mesh','on')
M(k)=getframe;
%}

%}
