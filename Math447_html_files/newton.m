clear all;
close all;

% Project 2

%--------------------------------------------------------------------------
%                             Define values
%--------------------------------------------------------------------------

parts=3; %use A,B,C if 1, rho, phi, theta if 2, do part 5 if 3, do part 5 plot mode if 4


c=299792.458;

% x = 0;
% y = 0;
% z = 6370;
% d = 0;

x = 0;
y = 0;
z = 6370;
d = 0.0001;

if (parts==1)|(parts==2)|(parts==3)

if parts==1
A1=15600 ;B1=7540 ;C1=20140 ;t1=0.07074 ;
A2=18760 ;B2=2750 ;C2=18610 ;t2=0.07220 ;
A3=17610 ;B3=14630 ;C3=13480 ;t3=0.07690 ;
A4=19170 ;B4=610 ;C4=18390 ;t4=0.07242 ;
end

if parts==2;
    
    %----------------------------------------------------------------------
    %                         Part 4
    %----------------------------------------------------------------------
    
    rho=26570;
    phi1=abs(rand())*pi/2
    phi2=abs(rand())*pi/2
    phi3=abs(rand())*pi/2
    phi4=abs(rand())*pi/2
    theta1=abs(rand())*2*pi
    theta2=abs(rand())*2*pi
    theta3=abs(rand())*2*pi
    theta4=abs(rand())*2*pi
    
    A1=rho*cos(phi1)*cos(theta1);
    A2=rho*cos(phi2)*cos(theta2);
    A3=rho*cos(phi3)*cos(theta3);
    A4=rho*cos(phi4)*cos(theta4);
    
    B1=rho*cos(phi1)*sin(theta1);
    B2=rho*cos(phi2)*sin(theta2);
    B3=rho*cos(phi3)*sin(theta3);
    B4=rho*cos(phi4)*sin(theta4);
    
    C1=rho*sin(theta1);
    C2=rho*sin(theta2);
    C3=rho*sin(theta3);
    C4=rho*sin(theta4);
    
    R1=sqrt(A1^2+B1^2+(C1-6370)^2);
    R2=sqrt(A2^2+B2^2+(C2-6370)^2);
    R3=sqrt(A3^2+B3^2+(C3-6370)^2);
    R4=sqrt(A4^2+B4^2+(C4-6370)^2);
    
    t1=d+R1/c;
    t2=d+R2/c;
    t3=d+R3/c;
    t4=d+R4/c;
end

if parts==3;
    
    %----------------------------------------------------------------------
    %                         Part 5 a
    %----------------------------------------------------------------------
    
    rho=26570;
    phi1=abs(rand())*pi/2
    phi2=phi1+phi1*rand()*(5/100)
    phi3=phi1+phi1*rand()*(5/100)
    phi4=phi1+phi1*rand()*(5/100)
    
    theta1=abs(rand())*2*pi
    theta2=theta1+theta1*rand()*(5/100)
    theta3=theta1+theta1*rand()*(5/100)
    theta4=theta1+theta1*rand()*(5/100)
    
    A1=rho*cos(phi1)*cos(theta1);
    A2=rho*cos(phi2)*cos(theta2);
    A3=rho*cos(phi3)*cos(theta3);
    A4=rho*cos(phi4)*cos(theta4);
    
    B1=rho*cos(phi1)*sin(theta1);
    B2=rho*cos(phi2)*sin(theta2);
    B3=rho*cos(phi3)*sin(theta3);
    B4=rho*cos(phi4)*sin(theta4);
    
    C1=rho*sin(theta1);
    C2=rho*sin(theta2);
    C3=rho*sin(theta3);
    C4=rho*sin(theta4);
    
    R1=sqrt(A1^2+B1^2+(C1-6370)^2);
    R2=sqrt(A2^2+B2^2+(C2-6370)^2);
    R3=sqrt(A3^2+B3^2+(C3-6370)^2);
    R4=sqrt(A4^2+B4^2+(C4-6370)^2);
    
    t1=d+R1/c;
    t2=d+R2/c;
    t3=d+R3/c;
    t4=d+R4/c;
end

num=10;     % the number of times Newton's method will be run

%--------------------------------------------------------------------------
%                      Multivariate Newton's Method
%--------------------------------------------------------------------------

for times=0:1;
    t1=t1+times*10^(-8);
    
    
for i=1:num
 
F1 = ((x-A1)^2+(y-B1)^2+(z-C1)^2)-(c*(t1-d))^2;
F2 = ((x-A2)^2+(y-B2)^2+(z-C2)^2)-(c*(t2-d))^2;
F3 = ((x-A3)^2+(y-B3)^2+(z-C3)^2)-(c*(t3-d))^2;
F4 = ((x-A4)^2+(y-B4)^2+(z-C4)^2)-(c*(t4-d))^2;

Re = [2*(x-A1), 2*(y-B1), 2*(z-C1), 2*(c*(t1-d))*c, 1, 0, 0, 0;...
      2*(x-A2), 2*(y-B2), 2*(z-C2), 2*(c*(t2-d))*c, 0, 1, 0, 0;...
      2*(x-A3), 2*(y-B3), 2*(z-C3), 2*(c*(t3-d))*c, 0, 0, 1, 0;...
      2*(x-A4), 2*(y-B4), 2*(z-C4), 2*(c*(t4-d))*c, 0, 0, 0, 1];

%-----------------------Gaussian Elimination-------------------------------
  
for col=1:3;
    for row=1+col:4;
        Re(col,:)=Re(col,:)./Re(col,col);
        Re(row,:)=Re(row,:) - Re(row,col)*Re(col,:);
    end
end

Re(4,:)=Re(4,:)/Re(4,4);

for col=4:-1:2
    for row=col-1:-1:1
        Re(row,:)=Re(row,:)-Re(row,col)*Re(col,:);
    end
end

%--------------------------------------------------------------------------

Rinv=Re(:,5:8);

v = [x,y,z,d]' - Rinv*[F1,F2,F3,F4]';

end

if times==0;
xold=v(1);
yold=v(2);
zold=v(3);
dold=v(4);
end

if times==1;
xnew=v(1);
ynew=v(2);
znew=v(3);
dnew=v(4);
end


end

maximum_position_error=max([abs(xnew-xold),abs(ynew-yold),abs(znew-zold)])
error_magnification_factor=maximum_position_error/(c*10^(-8))

end


if parts==4;
    
    %----------------------------------------------------------------------
    %                         Part 5 b
    %----------------------------------------------------------------------
    
    rho=26570;
    phi1=.2*pi/2
    phi2=phi1+phi1*rand()*(5/100)
    phi3=phi1+phi1*rand()*(5/100)
    phi4=phi1+phi1*rand()*(5/100)
    
    theta1=.1*2*pi
    theta2=theta1+theta1*rand()*(5/100)
    theta3=theta1+theta1*rand()*(5/100)
    theta4=theta1+theta1*rand()*(5/100)
    
    
    ttmax=100
    tmax=500;
    
    
    for drawing=1:ttmax
        
        phi1=.2*pi/2;
%         phi2=phi1+1/ttmax*.8*pi/2;
%         phi3=phi1+1/ttmax*.8*pi/2;
%         phi4=phi1+1/ttmax*.8*pi/2;
%         theta1=theta1+1/ttmax*.9*pi*2;
        theta2=theta1+1/ttmax*.9*pi*2;
%         theta3=theta1+1/ttmax*.9*pi*2;
%         theta4=theta1+1/ttmax*.9*pi*2;
        
    for run=1:tmax
    
    phi1=phi1+1/tmax*.8*pi/2;
    phi1_matrix(1,run)=phi1;
    
    A1=rho*cos(phi1)*cos(theta1);
    A2=rho*cos(phi2)*cos(theta2);
    A3=rho*cos(phi3)*cos(theta3);
    A4=rho*cos(phi4)*cos(theta4);
    
    B1=rho*cos(phi1)*sin(theta1);
    B2=rho*cos(phi2)*sin(theta2);
    B3=rho*cos(phi3)*sin(theta3);
    B4=rho*cos(phi4)*sin(theta4);
    
    C1=rho*sin(theta1);
    C2=rho*sin(theta2);
    C3=rho*sin(theta3);
    C4=rho*sin(theta4);
    
    R1=sqrt(A1^2+B1^2+(C1-6370)^2);
    R2=sqrt(A2^2+B2^2+(C2-6370)^2);
    R3=sqrt(A3^2+B3^2+(C3-6370)^2);
    R4=sqrt(A4^2+B4^2+(C4-6370)^2);
    
    t1=d+R1/c;
    t2=d+R2/c;
    t3=d+R3/c;
    t4=d+R4/c;
    
    
        
            
    num=10;     % the number of times Newton's method will be run

%--------------------------------------------------------------------------
%                      Multivariate Newton's Method
%--------------------------------------------------------------------------

        for times=0:1;
    t1=t1+times*10^(-8);
    
    
for i=1:num
 
F1 = ((x-A1)^2+(y-B1)^2+(z-C1)^2)-(c*(t1-d))^2;
F2 = ((x-A2)^2+(y-B2)^2+(z-C2)^2)-(c*(t2-d))^2;
F3 = ((x-A3)^2+(y-B3)^2+(z-C3)^2)-(c*(t3-d))^2;
F4 = ((x-A4)^2+(y-B4)^2+(z-C4)^2)-(c*(t4-d))^2;

Re = [2*(x-A1), 2*(y-B1), 2*(z-C1), 2*(c*(t1-d))*c, 1, 0, 0, 0;...
      2*(x-A2), 2*(y-B2), 2*(z-C2), 2*(c*(t2-d))*c, 0, 1, 0, 0;...
      2*(x-A3), 2*(y-B3), 2*(z-C3), 2*(c*(t3-d))*c, 0, 0, 1, 0;...
      2*(x-A4), 2*(y-B4), 2*(z-C4), 2*(c*(t4-d))*c, 0, 0, 0, 1];

%-----------------------Gaussian Elimination-------------------------------
  
for col=1:3;
    for row=1+col:4;
        Re(col,:)=Re(col,:)./Re(col,col);
        Re(row,:)=Re(row,:) - Re(row,col)*Re(col,:);
    end
end

Re(4,:)=Re(4,:)/Re(4,4);

for col=4:-1:2
    for row=col-1:-1:1
        Re(row,:)=Re(row,:)-Re(row,col)*Re(col,:);
    end
end

%--------------------------------------------------------------------------

Rinv=Re(:,5:8);

v = [x,y,z,d]' - Rinv*[F1,F2,F3,F4]';

end

if times==0;
xold=v(1);
yold=v(2);
zold=v(3);
dold=v(4);
end

if times==1;
xnew=v(1);
ynew=v(2);
znew=v(3);
dnew=v(4);
end


        end
        
maximum_position_error=max([abs(xnew-xold),abs(ynew-yold),abs(znew-zold)]);
error_magnification_factor=maximum_position_error/(c*10^(-8));
maximum_position_error_matrix(1,run)=maximum_position_error;
error_magnification_factor_matrix(1,run)=error_magnification_factor;

        end

    plot(phi1_matrix,error_magnification_factor_matrix)
    str = sprintf('phi2 = %f, phi3 = %f, phi4 = %f, theta1 = %f, theta2 = %f, theta3 = %f, theta4 = %f',phi2,phi3,phi4,theta1,theta2,theta3,theta4);
    title(str);

drawnow;
    end
    
end