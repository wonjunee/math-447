% Project 2 Part 4

clear all;

%--------------------------------------------------------------------------
%                             Define values
%--------------------------------------------------------------------------
bunched_or_unbunched=0; % 0=bunched 1=unbunched


c=299792.458;

x = 0;
y = 0;
z = 6370;
d = 0;

%     x = 0;
%     y = 0;
%     z = 6370;
%     d = 0.0001;

    rho=26570;
    
%     
%     phi1=abs(rand())*pi/2
%     phi2=phi1+phi1*rand()*(5/100)
%     phi3=phi1+phi1*rand()*(5/100)
%     phi4=phi1+phi1*rand()*(5/100)
%     phi5=phi1+phi1*rand()*(5/100)
%     phi6=phi1+phi1*rand()*(5/100)
%     phi7=phi1+phi1*rand()*(5/100)
%     phi8=phi1+phi1*rand()*(5/100)
% 
%     theta1=abs(rand())*2*pi
%     theta2=theta1+theta1*rand()*5/100
%     theta3=theta1+theta1*rand()*5/100
%     theta4=theta1+theta1*rand()*5/100
%     theta5=theta1+theta1*rand()*5/100
%     theta6=theta1+theta1*rand()*5/100
%     theta7=theta1+theta1*rand()*5/100
%     theta8=theta1+theta1*rand()*5/100
%     
    phi1=abs(rand())*pi/2
    phi2=abs(rand())*pi/2
    phi3=abs(rand())*pi/2
    phi4=abs(rand())*pi/2
    phi5=abs(rand())*pi/2
    phi6=abs(rand())*pi/2
    phi7=abs(rand())*pi/2
    phi8=abs(rand())*pi/2

    theta1=abs(rand())*2*pi
    theta2=abs(rand())*2*pi
    theta3=abs(rand())*2*pi
    theta4=abs(rand())*2*pi
    theta5=abs(rand())*2*pi
    theta6=abs(rand())*2*pi
    theta7=abs(rand())*2*pi
    theta8=abs(rand())*2*pi
    
    A1=rho*cos(phi1)*cos(theta1);
    A2=rho*cos(phi2)*cos(theta2);
    A3=rho*cos(phi3)*cos(theta3);
    A4=rho*cos(phi4)*cos(theta4);
    A5=rho*cos(phi5)*cos(theta5);
    A6=rho*cos(phi6)*cos(theta6);
    A7=rho*cos(phi7)*cos(theta7);
    A8=rho*cos(phi8)*cos(theta8);
    
    B1=rho*cos(phi1)*sin(theta1);
    B2=rho*cos(phi2)*sin(theta2);
    B3=rho*cos(phi3)*sin(theta3);
    B4=rho*cos(phi4)*sin(theta4);
    B5=rho*cos(phi5)*sin(theta5);
    B6=rho*cos(phi6)*sin(theta6);
    B7=rho*cos(phi7)*sin(theta7);
    B8=rho*cos(phi8)*sin(theta8);
    
    C1=rho*sin(theta1);
    C2=rho*sin(theta2);
    C3=rho*sin(theta3);
    C4=rho*sin(theta4);
    C5=rho*sin(theta5);
    C6=rho*sin(theta6);
    C7=rho*sin(theta7);
    C8=rho*sin(theta8);
    
    R1=sqrt(A1^2+B1^2+(C1-6370)^2);
    R2=sqrt(A2^2+B2^2+(C2-6370)^2);
    R3=sqrt(A3^2+B3^2+(C3-6370)^2);
    R4=sqrt(A4^2+B4^2+(C4-6370)^2);
    R5=sqrt(A5^2+B5^2+(C5-6370)^2);
    R6=sqrt(A6^2+B6^2+(C6-6370)^2);
    R7=sqrt(A7^2+B7^2+(C7-6370)^2);
    R8=sqrt(A8^2+B8^2+(C8-6370)^2);
    
    
    t1=d+R1/c;
    t2=d+R2/c;
    t3=d+R3/c;
    t4=d+R4/c;
    t5=d+R4/c;
    t6=d+R4/c;
    t7=d+R4/c;
    t8=d+R4/c;
   
    
%--------------------------------------------------------------------------    
    
    
told=t1;

plotpoints=500; %number of points plotted


    for run=1:plotpoints
    
    xxaxis(run,1)=t1-told;
    
%--------------------------------------------------------------------------
%                      Multivariate Newton's Method
%--------------------------------------------------------------------------
x = 0;
y = 0;
z = 6370;

num=20;     % the number of times Newton's method will be run

for i=1:num
 
F1 = ((x-A1)^2+(y-B1)^2+(z-C1)^2)-(c*(t1-d))^2;
F2 = ((x-A2)^2+(y-B2)^2+(z-C2)^2)-(c*(t2-d))^2;
F3 = ((x-A3)^2+(y-B3)^2+(z-C3)^2)-(c*(t3-d))^2;
F4 = ((x-A4)^2+(y-B4)^2+(z-C4)^2)-(c*(t4-d))^2;
F5 = ((x-A5)^2+(y-B5)^2+(z-C5)^2)-(c*(t5-d))^2;
F6 = ((x-A6)^2+(y-B6)^2+(z-C6)^2)-(c*(t6-d))^2;
F7 = ((x-A7)^2+(y-B7)^2+(z-C7)^2)-(c*(t7-d))^2;
F8 = ((x-A8)^2+(y-B8)^2+(z-C8)^2)-(c*(t8-d))^2;

Re = [2*(x-A1), 2*(y-B1), 2*(z-C1), 2*(c*(t1-d))*c;...
      2*(x-A2), 2*(y-B2), 2*(z-C2), 2*(c*(t2-d))*c;...
      2*(x-A3), 2*(y-B3), 2*(z-C3), 2*(c*(t3-d))*c;...
      2*(x-A4), 2*(y-B4), 2*(z-C4), 2*(c*(t4-d))*c;...
      2*(x-A5), 2*(y-B5), 2*(z-C5), 2*(c*(t5-d))*c;...
      2*(x-A6), 2*(y-B6), 2*(z-C6), 2*(c*(t6-d))*c;...
      2*(x-A7), 2*(y-B7), 2*(z-C7), 2*(c*(t7-d))*c;...
      2*(x-A8), 2*(y-B8), 2*(z-C8), 2*(c*(t8-d))*c;];
  
RE=Re'*Re;
REinv=inv(RE);
Rinv=REinv*Re';

%Outputs from Newton's Method

v = [x,y,z,d]' - Rinv*[F1,F2,F3,F4,F5,F6,F7,F8]';
   

x=v(1);
y=v(2);
z=v(3);
d=v(4);

end


if run==1;
    xold=v(1);
    yold=v(2);
    zold=v(3);
    dold=v(4);
end

if run~=1;
    xnew(run-1,:)=v(1);
    ynew(run-1,:)=v(2);
    znew(run-1,:)=v(3);
    dnew(run-1,:)=v(4);
    
    maximum_position_error(run,:)=max([abs(xnew(run-1,:)-xold),abs(ynew(run-1,:)-yold),abs(znew(run-1,:)-zold)]);
    
    error_magnification_factor(run,:)=maximum_position_error(run,:)/(c*xxaxis(run,1));
end
%--------------------------------------------------------------------------

t1=t1+10^(-10);

    end %end of run

%-------------- Making data of error magnification factor -----------------
    
maximum_position_error(1,:)=0;
error_magnification_factor(1,:)=0;

%--------------------------------------------------------------------------


%       Plot delta ti vs. error magnification factor
figure(1)
subplot(1,2,1)
plot(xxaxis(2:end),error_magnification_factor(2:end));

xlabel('delta t_1');
ylabel('Error Magnification Factor');

title('delta t_1 vs. Error Magnification Factor')

%       Plot delta ti vs. error magnification factor
subplot(1,2,2)
plot(xxaxis(2:end),maximum_position_error(2:end));

xlabel('delta t_1');
ylabel('Maximum Position Error');

title('delta t_1 vs. Maximum Position Error')