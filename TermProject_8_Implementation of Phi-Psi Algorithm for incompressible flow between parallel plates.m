%Flow through parallel plates

% setting up the parameters
clear all
clc
Lf = 2;    %length of the free stream
Lp = 2;    %length of the plate
II = 50;      %number of nodes along the x axis 
IO = 25;      %number of nodes in free stream 
JJ = 25;      %number of nodes along the Y direction
h=1;
Re = 1;      % reynolds number
E_shi = 0.000001;    % error in SHI
E_omg = 0.000001;  %error1 in OMEGA 
Ly= 1;     %distance between the plate and the center line 
Lx = Lp+Lf; %total length of the flow simulation


% calculating the values of Delta_X, Beta, F
Del_x = (Lx)/II; % grid space in x direc
Del_y = Ly/JJ;       % grid space in Y direc 
B  = (Del_x/Del_y);% Aspect ratio, Beta
zeta  =  (cos((pi*Del_x/Lx))+(B^2)*sin((pi*Del_y)/Ly))/(((B^2) +1)^2);
F     =  2*(1-sqrt(1-(zeta^2)))/zeta ; %relaxation parameter 

Del_t =  0.5/((1+B)/Del_x+4*(1+B^2)/(Re*(Del_x^2))); % delta t (time)


% setting up initial and boundary conditions
%u = 1; v = 0; omg = 0; shi = y ; % flow at -inf is assumed to be at II=1


u = zeros(II+1,JJ+1);

v= ones(II+1,JJ+1);

vor = zeros(II+1,JJ+1);
for i= 1:II+1
    for j= 1:JJ+1
    shi = (j-1)*Del_y; 
    end
end

% iitializing the parameters at the boundary
% initializing at the wall
for i=IO+1:II+1
u(i,JJ+1)= 0;
v(i,JJ+1)= 0;
shi(i,JJ+1)= 1;
vor(i,JJ+1) = 3;
end

%initialinzing at the outlet

for j= 1:JJ+1
    y = (j-1)*Del_y;
    u(II+1,j)= 1.5*y-0.5*(y^2);
    v(II+1,j)= 0;
    shi(II+1,j)= 1.5*y-0.5*(y^3);
    vor(II+1,j)= 3*y;
end
%solving for vorticity at interior points
count=0;
loop1 = 0;
loop2 = 0;
while count<1
loop1 = loop1+1;    
for i=2:II
    for j= 2:JJ
            %count=0;
        vor_max = 1; %vor(IO+1,JJ+1);
        %while count<1
        if u(i,j)>=0
            Delx_uw = (u(i,j)*vor(i,j)-u(i-1,j)*vor(i-1,j))/Del_x;
        else
            Delx_uw = (u(i+1,j)*vor(i+1,j)-u(i,j)*vor(i,j))/Del_x;
        end
        
        if v(i,j)>=0
            Dely_vw = (v(i,j)*vor(i,j)-v(i,j-1)*vor(i,j-1))/Del_y;
        else
            Dely_vw = (v(i,j+1)*vor(i,j+1)-v(i,j)*vor(i,j))/Del_y;
        end
        Delsq_vor = (vor(i+1,j)+vor(i-1,j)+(B^2)*vor(i,j+1)-2*((B^2)+1)*vor(i,j))/(Del_x^2);
        
        Del_vor = Del_t*(-Delx_uw-Dely_vw+2*(Delsq_vor)/Re);
        
        vorn(i,j)= vor(i,j)+ Del_vor;
        
        if abs((vorn(i,j)-vor(i,j))/vor_max) < E_omg
            count = 1;
        else 
            count =0;
        end
            
            vor(i,j)= vorn(i,j);
        end
   end

 

 % solving for stream function in the interior
 count2 = 0;
 while count2<1
     
    count2 = 
    loop2 = loop2+1;
     for i=2:II
          for j=2:JJ
            %count2 = 0;
            shi_max = 1;
            %while count2<1
            DSTR = shi(i+1,j)+shi(i-1,j)+(B^2)*shi(i,j+1)+(B^2)*shi(i,j-1)-2*((B^2)+1)*shi(i,j)+vor(i,j)*(Del_x^2);
     
            shin(i,j)= shi(i,j)+F*DSTR/(2*((B^2)+1));
     
            if abs((shin(i,j)-shi(i,j))/shi_max)<= E_shi
         
                count2= 1;
           % else
            %    count2=0;
            end
    
            shi(i,j)=shin(i,j);
          end
 
      end
 end
 % calculating the velocities in the interior 
 
 for i=2:II
 for j=2:JJ
 u(i,j)= B*(shi(i,j+1)-shi(i,j-1))/(2*Del_x);
 v(i,j)= (shi(i-1,j)-shi(i+1,j))/(2*Del_x);
 end
 end
  
 % clculating center line and stagnation velocity u 
 
 for i= 2:II
 u(i,1)= shi(i,2)*B/Del_x;
 end
 for i= 2:IO+1
 u(i,JJ+1)= (1-shi(i,JJ))*B/Del_x;
 end
 %calculating vorticity at the walls
for i= IO+2:II 
 %count3 = 0;
    %while count <1
    vor_temp = (1-shi(i,JJ)*(B^2)/(Del_x^2));
    Dvor=vor(i,JJ+1)-vor_temp;
    if (abs(Dvor/vor_max)<E_omg)
        count=1;
    else
        count=0;
    end
    
end
 % crak nikoklson scheme for leading edge vorticity
 vor(IO+1,JJ+1)=0.5*(vor(IO+2,JJ+1)+vor(IO+2,JJ-1)+vor(IO,JJ-1))-vor(IO+2,JJ)-vor(IO+1,JJ-1)-vor(IO,JJ)+2*vor(IO+1,JJ);
end

%printing of the values
 
     for j=1:JJ+1
         y=(j-1)*Del_y;
         plot(u(60,j),y(j))
     end

     
     
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
     
 










        
        
        
         
            
            
            
            
            
            
   





















