clear all
%%%%%%%%%%%%%2-D ADVECTION CRANK-NICOLSON METHOD%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create a cyclonic wind field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:25
    for j=1:25
        x=i-13;
        y=j-13;
        if x>=0 && y>0
        theta(j,i)=rad2deg(atan(y/x));
        elseif x<0 && y>=0
        theta(j,i)=180+rad2deg(atan(y/x));
        elseif x<0 && y<0
        theta(j,i)=180+rad2deg(atan(y/x));
        elseif x>=0 && y<0
        theta(j,i)=360+rad2deg(atan(y/x));    
        end
            r=sqrt(x^2 + y^2)*4; %km
            u(j,i)=-sind(theta(j,i)) * 0.08*pi*r; %km/hr
            v(j,i)=cosd(theta(j,i)) * 0.08*pi*r; %km/hr
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('This is the wind field')
% figure
% quiver(u,v)
% title('u-v wind vectors')
% xlabel('x')
% ylabel('y')
% pause(5) %pause for 5 seconds

%set parameters and conditions%%%%%%%%%%%%%%%%%%%
dt=5/60; %5min to hours
dx=4; %km
% %Initial conditions of 2-d pollution plume
x1=zeros(25,25);
x1(3:7,3:7)=100;
matmax=23; %for indexing purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_total=sum(sum(x1));

for t=1:300 %go through time steps
    
    %%%%%COLUMN V-COMPONENT
    for a=2:24
    B= v(1,a)*dt/(4*dx);
    
    knowns_v(1,a-1)=x1(2,a) + (-B*x1(3,a)) + 0;
    for i=3:matmax
    knowns_v(i-1,a-1)=(B*x1(i-1,a)) + x1(i,a) + (-B*x1(i+1,a));
    end
    knowns_v(23,a-1)=(B*x1(23,a)) + x1(24,a) + 0;
    
    for z=1:23
    coeff(z,z)=1;
    end
    clear z
    for z=1:matmax-1
    coeff(z,z+1)=B; %next 2,3
    coeff(z+1,z)=-B; %next 3,2
    end
    clear z
    
    x1(2:24,a)=coeff\knowns_v(:,a-1);
    
    end
    
    %%%%ROW U-COMPONENT
    for b=2:24
    B= u(b,1)*dt/(4*dx);
    
    knowns_u(b-1,1)=x1(b,2) + (-B*x1(b,3)) + 0;   
    for j=3:matmax
    knowns_u(b-1,j-1)=(B*x1(b,j-1)) + x1(b,j) + (-B*x1(b,j+1));
    end    
    knowns_u(b-1,23)=(B*x1(b,23)) + x1(b,24) + 0;
    
    
    for z=1:23
    coeff(z,z)=1;
    end
    clear z
    for z=1:matmax-1
    coeff(z,z+1)=B; %next 2,3
    coeff(z+1,z)=-B; %next 3,2
    end
    clear z
    
    x1(b,2:24)=coeff\knowns_u(b-1,:)';
    end
    
    total=sum(sum(x1));
    mass_conservation_ratio(t)=total/initial_total;
    mass_distribution_ratio(t)=total^2 / initial_total^2;
    
    %create a 3-d surface plot every 50 time steps
    blah=rem(t,50); %remainder of t/50
    if t==1 || blah==0 %or statement
    surf(x1(1:25,1:25));figure(gcf);
    title(['Time step: ' num2str(t)])
    pause(1) %pause for 1 sec for illustrative purposes
    end
    
end

%Mass conservation plots
figure
plot(1:t,mass_conservation_ratio)
hold on
xlabel('Time Steps')
ylabel('Ratio')
title('Mass Conservation Ratio')
grid on

figure
plot(1:t,mass_distribution_ratio)
hold on
xlabel('Time Steps')
ylabel('Ratio')
title('Mass Distribution Ratio')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
