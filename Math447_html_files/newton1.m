% Project 2

%--------------------------------------------------------------------------
%                             Define values
%--------------------------------------------------------------------------

coordinate=1; %use A,B,C if 1, rho, phi, theta if 2


c=299792.458;

% x = 0;
% y = 0;
% z = 6370;
% d = 0;

xold = 0;
yold = 0;
zold = 6370;
dold = 0;

x = 0;
y = 0;
z = 6370;
d = 0.0001;



if coordinate==1
A1=15600 ;B1=7540 ;C1=20140 ;t1=0.07074 ;
A2=18760 ;B2=2750 ;C2=18610 ;t2=0.07220 ;
A3=17610 ;B3=14630 ;C3=13480 ;t3=0.07690 ;
A4=19170 ;B4=610 ;C4=18390 ;t4=0.07242 ;
end

if coordinate==2;
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

num=10;     % the number of times Newton's method will be run

%--------------------------------------------------------------------------
%                      Multivariate Newton's Method
%--------------------------------------------------------------------------


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

x=v(1);
y=v(2);
z=v(3);
d=v(4);
   
end

x
y
z
d