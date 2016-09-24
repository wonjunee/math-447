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
    
    
    %%%%%%% Bunched %%%%%%%%
    
    if bunched_or_unbunched==0;
    phi1=abs(rand())*pi/2
    phi2=phi1+phi1*rand()*(5/100)
    phi3=phi1+phi1*rand()*(5/100)
    phi4=phi1+phi1*rand()*(5/100)

    theta1=abs(rand())*2*pi
    theta2=theta1+theta1*rand()*5/100
    theta3=theta1+theta1*rand()*5/100
    theta4=theta1+theta1*rand()*5/100
    
    end
    
    
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

%Outputs from Newton's Method

v = [x,y,z,d]' - Rinv*[F1,F2,F3,F4]';
   

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
plot(xxaxis,maximum_position_error);

xlabel('delta t_1');
ylabel('Maximum Position Error');

title('delta t_1 vs. Maximum Position Error')