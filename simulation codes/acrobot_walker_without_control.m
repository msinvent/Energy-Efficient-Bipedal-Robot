% Programmer : Manish Sharma
% IIIT- Allahabad

function [] = manish_acrobot()
clear functions % This line is really needed, damn!
close all;
clear all;
clc;

global control;
T_final=500;             %duration of animation  in seconds
animate = 1;
control = 0;

%%%%%%%%% INITIALIZE PARAMETERS %%%%%%
%Mechanical parameters.
global m Mh L a g gamma;
m  =  2;    
Mh =  4;  % masses

L = 1;
a = 0.5;
g   =  9.81;
%gamma = 0.033; % max value of slope for stability without control
%gamma = .2;
%gamma = -0.002
%gamma = 0.030; % at kp = 50 stable
%gamma = 0.020; % at kp = 50 bifurcation observed
%gamma = 0.02; % at kp = 25 stable
%gamma = 0.01; % at kp = 25 stable but probably not converged
%gamma = 0.01; % at kp = 5 stable
%gamma = 0.0; % at kp = 5 stable
%gamma = -0.01; % at kp = -25 stable
%gamma = 0.20
%gamma = 0.02;
gamma = 0.2;


global kp;


global tou_t t_t;
tou_t = [];
t_t = [];
% Initial conditions and other settings.
framespersec=10;  %if view is not speeded or slowed in dbpend_animate
tspan=linspace(0,T_final,T_final*framespersec);

q1    = pi/2+0.2; %angle made by link1 with vertical
if(control == 1)
    q1    = pi/2+0.1; %angle made by link1 with vertical
end
u1    = -1.0;        %abslolute velocity of link1   
q2    = pi-0.4;      %angle made by link2 with vertical
u2    = +.6;        %abslolute velocity of link2



x0=[q1;u1;q2;u2];

options=odeset('abstol',1e-13,'reltol',1e-13,'Events',@events);

%%%%%%% INTEGRATOR or ODE SOLVER %%%%%%%
x_new = 0;
te = [];xe =[];toue=[];
z1 = 0;
z2 = 0;
% steps taking
frequency = 200; % Hz
global counter control_counter_max flag tou;
control_counter_max = 40; %should be an even number
tspan = [0,1/frequency];

counter = 1;
i = 0;
te = 0;
xe = [x0' 0 0];
vel_horiz = [],vel_vert = [];
t22 = [];

% gamma = rand/20
    kp = 2500*(gamma-0.01);
    if(kp >100)
        kp = 100;
    end
    kp
%kp = 5
steps_count = 10;
while (i<steps_count)
    
    
    %flag
    %this flag is used to avoid update of control value while in between
    %the tspan time period
    %control var is edited when flag == 1 and not edited when flag ==0
    %in the first iteration itself flag is set to zero from 1
    
    flag = 1;
    [t,x,TE,YE,IE] = ode45(@db_manish,tspan,x0,options);
    %counter
    %increment counter after each sampling period it keeps the track when
    %the control has to be switched off 
    counter = counter+1;
    toue = [toue; tou];
    
    x = [x(1,:);x(end,:)];
    x = [x,repmat(z1,2,1),repmat(z2,2,1)];
    te = [te ; t(end)];
    xe = [xe ; x(end,:)];
    q1 = x(end,1);
    dq1 = x(end,2);
    horiz_velocity = L*dq1*cos(q1-pi/2);
    vert_velocity = L*dq1*sin(q1-pi/2);
    vel_horiz = [vel_horiz horiz_velocity];
    vel_vert = [vel_vert vert_velocity];
    t22 = [t22 t(end)];
    if(isempty(IE))
        x0 = xe(end,1:4);
        tspan = [te(end),te(end)+1/frequency];
        continue;
    end
    
    if(IE == 2)
        display('Model crashed');
        break;
    end
    
    q1 = xe(end,1);q2 = xe(end,3);
    t_final2 = (ceil(te(end)*frequency))/frequency;
    if(t_final2 == te(end))
        t_final2 = te(end) + 1/frequency;
    end
    tspan=[te(end),t_final2];
    x0 = impact_grizzly(xe(end,:));
    x0 = x0(1:4);
    
    te = te(1:end);
    xe = xe(1:end,:);
    
    %tspan=linspace(te(end),T_final,(T_final)*framespersec);
    z1 = z1 + L*cos(q1) + L*cos(q1+q2); %position of pivot
    z2 = z2 + L*sin(q1) + L*sin(q1+q2);
    
    if(te(end) >= T_final - 1/frequency )
        display('Final time end');
        break
    end
    
    %reinitialize control counter
    counter = 1;
    %save torque values
    i= i+1;
    i
    display('Step complete');
%     gamma = rand/20
%     kp = 5^(gamma*100);
%     if(kp >100)
%         kp = 100;
%     end
%     kp  
%    
end


%%%%%%% POSTPROCESSING %%%%%%%%%%%%%%%%
% A routine to animate the results
% To speed up change the framespersecond
if(animate == 1)
    figure(1)
    
    for i=1:2:length(te)-2
                
        pause(te(i+2)-te(i));
                
        q1 = xe(i,1);
        q2 = xe(i,3);
        z1 = xe(i,5);
        z2 = xe(i,6);
        
        xm1=z1+L*cos(q1);
        ym1=z2+L*sin(q1);
        x_m00 = z1+(L-a)*cos(q1);
        y_m00 = z2+(L-a)*sin(q1);
        xm2=xm1+L*cos(q1+q2);
        ym2=ym1+L*sin(q1+q2);
        x_m11 = xm1 + a*cos(q1+q2);
        y_m11 = ym1 + a*sin(q1+q2);
        
        plot(z1,z2,'ko','MarkerSize',3); %pivot point
        hold on
        plot([z1 xm1],[z2 ym1],'r','LineWidth',2);% first pendulum
        plot(x_m00,y_m00,'ko','MarkerSize',5); %pivot m1 point
        plot(xm1,ym1,'ko','MarkerSize',8); %pivot Mh point
        plot([xm1 xm2],[ym1 ym2],'b','LineWidth',2);% second pendulum
        plot(x_m11,y_m11,'ko','MarkerSize',5); % m2 pivot point
        rr1 = xm1 + L*sin(pi/2-q1);
        rr2 = ym1 - L*cos(pi/2-q1);
        plot([z1-10*cos(gamma),z1+10*cos(gamma)],[z2+10*sin(gamma),z2-10*sin(gamma)]) %plot floor
        %plot([-10*cos(gamma),+10*cos(gamma)],[+10*sin(gamma),-10*sin(gamma)]) %plot floor
        %axis([-1*L 8*L -1*L 1*L]);
        axis([xm1-1.5*L xm1+1.5*L ym1-1.5*L ym1+0.5*L]);
        %axis square
        hold off
        
    end
end
%energy plot
for i=1:length(te)
    [TE_body(i) TE_inertial(i)] = energy(xe(i,:));
end
%TE_diff = diff(TE);
figure(2)
%plot(TE_diff(1:end-1));
plot(te,TE_body,'r',te,TE_inertial,'b','LineWidth',2);
legend('Total Energy(Body Frame)','Total Energy(inertial frame)','Location','SouthWest');
title('Energy Plot');

%limit cycle plot
figure(3);
plot(xe(:,1),xe(:,2),'LineWidth',2);
title('limit cycle 1');
xlabel('q1_{abs}');ylabel('q1_{abs_{dot}}');

figure(4);
plot(xe(:,3),xe(:,4),'LineWidth',2);
title('limit cycle 2');
xlabel('q2_{rel}');ylabel('q2_{rel_{dot}}');

%stride length plot
figure(5);
plot(te(1:end-1),diff(xe(:,5)),'LineWidth',2);
title('stride length');
xlabel('time(s)');ylabel('length(m)');

%joint velocities
figure(6)
plot(te,xe(:,2),te,xe(:,4),'LineWidth',2);
legend('q1_{abs_{dot}}','q2_{rel_{dot}}')
title('Joint Velocities');


figure(7)
plot(t22,vel_horiz,t22,vel_vert,'LineWidth',2)
legend('horizontal velocity','vertical velocity');
title('hip velocities');

%hip trajectory
figure(8)
plot(te,xe(:,5)+L*cos(xe(:,1)),'r','LineWidth',2)
hold on;
plot(te,xe(:,6)+L*sin(xe(:,1)),'b','LineWidth',2)
xlabel('time(t)');ylabel('position(m)')
title('hip position in inertial frame');
legend('horizontal position','vertical position');

hold off
figure(9)
plot(te,L*cos(xe(:,1)),'r','LineWidth',2)
hold on;
plot(te,L*sin(xe(:,1)),'b','LineWidth',2)
xlabel('time(s)');ylabel('position(m)')
title('hip position in body frame');
legend('horizontal position','vertical position','Location','SouthWest');
axis([0 te(end) -0.3 1.1]);
%torque_plot
figure(10)
plot(te(1:length(te)-1,:),toue,'LineWidth',1)
xlabel('time(s)');ylabel('Torque');
title('torque Applied at hip');

%close all;
end


function zdot=db_manish(t,z)

global m Mh L a g;
global counter control_counter_max flag tou control gamma kp;
persistent v tou_copy;

q1 = z(1);dq1 = z(2);q2 = z(3);dq2 = z(4);  

D = ...
     [ (m*(2*L^2 + 4*cos(q2)*L*a + 2*a^2))/2 + L^2*Mh + m*(L - a)^2, a*m*(a + L*cos(q2));...
     a*m*(a + L*cos(q2)),               a^2*m];
M11 = D(1,1);
M12 = D(1,2);
M21 = D(2,1);
M22 = D(2,2);

C = ...
    [ -L*a*dq2*m*sin(q2), -L*a*m*sin(q2)*(dq1 + dq2);...
    L*a*dq1*m*sin(q2),                          0];
C11 = C(1,1);
C12 = C(1,2);
C21 = C(2,1);
C22 = C(2,2);

G = ...
    [g*m*(a*cos(q1 + q2) + L*cos(q1)) + L*Mh*g*cos(q1) + g*m*cos(q1)*(L - a);...
    a*g*m*cos(q1 + q2)];
G1 = G(1);
G2 = G(2);

horiz_velocity = L*dq1*cos(q1-pi/2);

optimal_horiz = -1;
if(counter == 1 && flag == 1)
    err = (horiz_velocity-optimal_horiz);
    %v = (-100+(gamma-0.10)*100)*err;% 
    %v = (-75)*err;% 
    %v = -100*err;
    v = -kp*err;
    tou_copy = -(M21/M11)*(C11*dq1 + C12*dq2 + G1) + C21*(dq1 + C22*dq2 + G2) +(M22 - (M21*M12)/M11)*v;
    flag = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%
if(counter <= control_counter_max/2)
    tou = tou_copy;
end
if(counter > control_counter_max/2 && counter <= control_counter_max)
    tou = -tou_copy;
end
if(counter >control_counter_max)
    tou = 0;
end
   
% to switch of the control
if(control == 0)
    tou = 0;
end

qdot = D \ (-C*[dq1;dq2]-G + [0;1]*tou);
zdot = [dq1;qdot(1);dq2;qdot(2)];

end

function [value,isterminal,direction] = events(t,z,flag)
global L gamma;
q1 = z(1);
q2 = z(3);
xm1=+L*cos(q1);
ym1=L*sin(q1);
xm2=xm1+L*cos(q1+q2);
ym2=ym1+L*sin(q1+q2);
x = (z(4)+z(2)>0); % absolute velocity of link 2 should be negative

value(1) = ym2 + x + xm2*tan(gamma);     % detect height of impact= 0
value(2) = ym1-0.2;               %height of hip

isterminal(1) = 1;   % stop the integration
isterminal(2) = 1;   
direction(1) = -1;
direction(2) = -1;   % negative direction
end


function x_new=impact_grizzly(x)

global m Mh L a g;
q1=x(1); q2=x(3);
De = ...
    [                      L^2*Mh + 2*L^2*m + 2*a^2*m - 2*L*a*m + 2*L*a*m*cos(q2), a*m*(a + L*cos(q2)), - m*(a*sin(q1 + q2) + L*sin(q1)) - m*sin(q1)*(L - a) - L*Mh*sin(q1), m*(a*cos(q1 + q2) + L*cos(q1)) + m*cos(q1)*(L - a) + L*Mh*cos(q1);...
                                                              a*m*(a + L*cos(q2)),               a^2*m,                                                   -a*m*sin(q1 + q2),                                                  a*m*cos(q1 + q2);...
      - (m*(2*a*sin(q1 + q2) + 2*L*sin(q1)))/2 - m*sin(q1)*(L - a) - L*Mh*sin(q1),   -a*m*sin(q1 + q2),                                                            Mh + 2*m,                                                                 0;...
        (m*(2*a*cos(q1 + q2) + 2*L*cos(q1)))/2 + m*cos(q1)*(L - a) + L*Mh*cos(q1),    a*m*cos(q1 + q2),                                                                   0,                                                          Mh + 2*m];


E = ...
    [ - L*sin(q1 + q2) - L*sin(q1), -L*sin(q1 + q2), 1, 0;...
        L*cos(q1 + q2) + L*cos(q1),  L*cos(q1 + q2), 0, 1];

%tmp_vec=inv([De -E';E zeros(2)])*[De*[x(2);x(4);zeros(2,1)];zeros(2,1)];  
tmp_vec=([De -E';E zeros(2)])\[De*[x(2);x(4);zeros(2,1)];zeros(2,1)];  

%position updates
x_new(1)=q1+q2-pi;
x_new(3)=2*pi - q2;
%velocities updates
x_new(2)=tmp_vec(1)+tmp_vec(2);
x_new(4)=-tmp_vec(2);
%forces at point of impact these values are generally
% used to verify the veracity of the impact model
% for more details refer the Eric Westerwellt and grizzle bipedal control book
x_new(5)=tmp_vec(3);
x_new(6)=tmp_vec(4);
                                         
end


function [TE_body TE_inertial] = energy(x)
    global m Mh L a g;
    q1 = x(1);
    q2 = x(3);
    dq1 = x(2);
    dq2 = x(4);
    y = x(6);
    KE = (m*(L^2*dq1^2 + 2*cos(q2)*L*a*dq1^2 + 2*cos(q2)*L*a*dq1*dq2 + a^2*dq1^2 + 2*a^2*dq1*dq2 + a^2*dq2^2))/2 + (L^2*Mh*dq1^2)/2 + (dq1^2*m*(L - a)^2)/2;
    PE_body = g*m*(a*sin(q1 + q2) + L*sin(q1)) + Mh*g*(L*sin(q1)) + g*m*((L - a)*sin(q1));
    PE_inertial = g*m*(a*sin(q1 + q2) + L*sin(q1) + y) + Mh*g*(L*sin(q1) + y) + g*m*((L - a)*sin(q1) + y);
    TE_body = KE + PE_body;
    TE_inertial = KE + PE_inertial;
end