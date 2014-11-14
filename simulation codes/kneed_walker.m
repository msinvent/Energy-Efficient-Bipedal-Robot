% Programmer : Manish Sharma

function [] = model_2_v5()
clear functions % This line is really needed, damn!
close all;
clear all;
clc;

T_final=20;             %duration of animation  in seconds
animate = 1;

%%%%%%%%% INITIALIZE PARAMETERS %%%%%%
%Mechanical parameters.
global m1 m2 Mh L1 L2 a1 a2 g gamma;
m1  =  0.5;    
m2  =  1.5;
Mh =  0.5;  % masses
L1 = 0.4;
L2 = 0.4;
a1 = 0.2;
a2 = 0.2;
g   =  9.81;
gamma = 0.02;

global tou_t t_t;
tou_t = [];
t_t = [];
vel_horiz = [],vel_vert = [];
t22 = [];

% Initial conditions and other settings.
framespersec=20; 
tspan=linspace(0,T_final,T_final*framespersec);
q1    = pi/2+0.2; %angle made by link1 with vertical
u1    = -1.256;        %abslolute velocity of link1   
q2    = pi-0.2;      %angle made by link2 with vertical
u2    = 2;        %abslolute velocity of link2
q3    = -0.1; %angle made by link1 with vertical
u3    = 0;        %abslolute velocity of link1   



x0=[q1;u1;q2;u2;q3;u3];
x0=[1.3775 -1.0383 3.4882 0.7872 0 0]';
x0=[1.4431 -0.7712 3.3925 0.1397 0 0]';

 
%options=odeset('abstol',1e-9,'reltol',1e-9);

%%%%%%% INTEGRATOR or ODE SOLVER %%%%%%%
x_new = 0;
te = 0;xe =[];
z1 = 0;
z2 = 0;
% steps taking
for i = 1:10
    %% 3 link swing till knee impact
    options=odeset('abstol',1e-9,'reltol',1e-9,'Events',@events_knee_impact);
    tspan=linspace(te(end),T_final,(T_final-te(end))*framespersec);
    [t,x,TE,YE,IE] = ode113(@db_manish_3link,tspan,x0,options);
    TE
    x = [x,repmat(z1,size(x,1),1),repmat(z2,size(x,1),1)];
    te = [te ; t(2:end)];
    xe = [xe ; x(2:end,:)];
    q1 = x(end,1);
    dq1 = x(end,2);
    horiz_velocity = (L1+L2)*dq1*cos(q1-pi/2);
    vert_velocity = (L1+L2)*dq1*sin(q1-pi/2);
    vel_horiz = [vel_horiz horiz_velocity];
    vel_vert = [vel_vert vert_velocity];
    t22 = [t22 t(end)];
    
    if(IE == 1)
        display('knee impact detected');
    end
    if(IE == 2)
        display('model crashed');
        break;
    end    
    if(te(end) + (2/framespersec) >= T_final ) %too close to the end
        display('Final time end')
        break
    end
    x0 = impact_manish_knee(xe(end,:));
    
    
    %% 2 link swing till toe impact
    options=odeset('abstol',1e-9,'reltol',1e-9,'Events',@events_toe_impact);
    tspan=linspace(te(end),T_final,(T_final-te(end))*framespersec);
    
    [t,x,TE,YE,IE] = ode113(@db_manish_2link,tspan,x0,options);
    TE
    x = [x,zeros(size(x,1),1),NaN(size(x,1),1),repmat(z1,size(x,1),1),repmat(z2,size(x,1),1)];
    te = [te ; t(2:end)];
    xe = [xe ; x(2:end,:)];
    if(IE == 1)
        display('toe impact detected');
    end
    if(IE == 2)
        display('model crashed');
        break;
    end
    % if t_end is very close dont continue
    if(te(end) + (2/framespersec) >= T_final ) % dont continue loop
        display('Final time end');
        break
    end
    display('iteration complete');
    x0 = impact_grizzly_toe(xe(end,:));
    x0 = [x0(1:4) 0 0];
    
    q1 = xe(end,1);q2 = xe(end,3);
    z1 = z1 + (L1+L2)*cos(q1) + (L1+L2)*cos(q1+q2); %position of pivot
    z2 = z2 + (L1+L2)*sin(q1) + (L1+L2)*sin(q1+q2);
    display(i);
end

te = te(2:end);
for i = 1:length(xe)
    if isnan(xe(i,6))
        xe(i,6) = 0;
    end
end



%%%%%%% POSTPROCESSING %%%%%%%%%%%%%%%%
% A routine to animate the results
% To speed up change the framespersecond
L = (L1 + L2)*0.5;
    figure(1)
if animate == 1
    for i=1:length(te)-1
        %animation at 1x
        pause((te(i+1)-te(i))/1);
                
        q1 = xe(i,1);
        q2 = xe(i,3);
        q3 = xe(i,5);
        z1 = xe(i,7);
        z2 = xe(i,8);
        
        xm1=z1+(L1+L2)*cos(q1);
        ym1=z2+(L1+L2)*sin(q1);
        x_m11 = z1+(L1-a1)*cos(q1);
        y_m11 = z2+(L1-a1)*sin(q1);
        x_m22 = z1+(L1+L2-a2)*cos(q1);
        y_m22 = z2+(L1+L2-a2)*sin(q1);
        x_m33 = xm1 + a2*cos(q1+q2);
        y_m33 = ym1 + a2*sin(q1+q2);
        x_m44 = xm1 + L2*cos(q1+q2) + a1*cos(q1+q2+q3);
        y_m44 = ym1 + L2*sin(q1+q2) + a1*sin(q1+q2+q3);
        
        xm2=xm1+L2*cos(q1+q2);
        ym2=ym1+L2*sin(q1+q2);
        xm3=xm2+L1*cos(q1+q2+q3);
        ym3=ym2+L1*sin(q1+q2+q3);
        
        
        plot(z1,z2,'ko','MarkerSize',3); %pivot point
        hold on
        plot([z1 xm1],[z2 ym1],'r','LineWidth',2);% first pendulum
        plot([xm1 xm2],[ym1 ym2],'b','LineWidth',2);% second pendulum
        plot([xm2 xm3],[ym2 ym3],'g','LineWidth',2);% second pendulum
        
        plot(x_m11,y_m11,'ko','MarkerSize',5); %pivot m1 point
        plot(x_m22,y_m22,'ko','MarkerSize',5); %pivot m2 point
        plot(xm1,ym1,'ko','MarkerSize',8); %pivot Mh point
        plot(x_m33,y_m33,'ko','MarkerSize',5); %pivot m3 point
        plot(x_m44,y_m44,'ko','MarkerSize',5); %pivot m4 point
        
        plot([z1-10*cos(gamma),z1+10*cos(gamma)],[z2+10*sin(gamma),z2-10*sin(gamma)]) %plot floor
        %axis([-3*L 3*L -3*L 3*L]);
        axis([xm1-2.5*L xm1+2.5*L ym1-2.5*L ym1+2.5*L]);
        axis square
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

figure(5);
plot(xe(:,5),xe(:,6),'LineWidth',2);
title('limit cycle 3');
xlabel('q3_{rel}');ylabel('q3_{rel_{dot}}');

%stride length plot
figure(6);
plot(te(1:end-1),diff(xe(:,7)),'LineWidth',2);
title('stride length');
xlabel('time(s)');ylabel('length(m)');

%joint velocities
figure(7)
plot(te,xe(:,2),te,xe(:,4),te,xe(:,6),'LineWidth',2);
legend('q1_{abs_{dot}}','q2_{rel_{dot}}','q3_{rel_{dot}}')
title('Joint Velocities');
figure(8)
plot(te,xe(:,2),'LineWidth',2);
legend('q1_{abs_{dot}}')
title('Joint Velocities');
figure(9)
plot(te,xe(:,4),'LineWidth',2);
legend('q2_{rel_{dot}}')
title('Joint Velocities');
figure(10)
plot(te,xe(:,6),'LineWidth',2);
legend('q3_{rel_{dot}}')
title('Joint Velocities');




%hip velocity
figure(11)
plot(t22,vel_horiz,t22,vel_vert,'LineWidth',2)
legend('horizontal velocity','vertical velocity');
title('hip velocities');

%hip trajectory
figure(12)
plot(te,xe(:,7)+(L1+L2)*cos(xe(:,1)),'r','LineWidth',2)
hold on;
plot(te,xe(:,8)+(L1+L2)*sin(xe(:,1)),'b','LineWidth',2)
xlabel('time(t)');ylabel('position(m)')
title('hip position in inertial frame');
legend('horizontal position','vertical position');

hold off
figure(13)
plot(te,(L1+L2)*cos(xe(:,1)),'r','LineWidth',2)
hold on;
plot(te,(L1+L2)*sin(xe(:,1)),'b','LineWidth',2)
xlabel('time(s)');ylabel('position(m)')
title('hip position in body frame');
legend('horizontal position','vertical position','Location','SouthWest');
axis([0 te(end) -0.3 1.1]);


end


function zdot=db_manish_3link(t,z)

global m1 m2 Mh L1 L2 a1 a2 g;

q1 = z(1);dq1 = z(2);q2 = z(3);dq2 = z(4);q3 = z(5);dq3 = z(6);  

D = ...
     [ L1^2*Mh + L2^2*Mh + 2*L1^2*m1 + 2*L1^2*m2 + 2*L2^2*m1 + 2*L2^2*m2 + 2*a1^2*m1 + 2*a2^2*m2 + 2*L1*L2*Mh + 2*L2^2*m1*cos(q2) + 2*L1*L2*m1 + 4*L1*L2*m2 - 2*L1*a1*m1 - 2*L1*a2*m2 - 2*L2*a2*m2 + 2*L1*a1*m1*cos(q2 + q3) + 2*L2*a1*m1*cos(q2 + q3) + 2*L1*L2*m1*cos(q2) + 2*L1*a2*m2*cos(q2) + 2*L2*a1*m1*cos(q3) + 2*L2*a2*m2*cos(q2), L2^2*m1 + a1^2*m1 + a2^2*m2 + L2^2*m1*cos(q2) + L1*a1*m1*cos(q2 + q3) + L2*a1*m1*cos(q2 + q3) + L1*L2*m1*cos(q2) + L1*a2*m2*cos(q2) + 2*L2*a1*m1*cos(q3) + L2*a2*m2*cos(q2), a1*m1*(a1 + L1*cos(q2 + q3) + L2*cos(q2 + q3) + L2*cos(q3));...
                                                                                                                                                               L2^2*m1 + a1^2*m1 + a2^2*m2 + L2^2*m1*cos(q2) + L1*a1*m1*cos(q2 + q3) + L2*a1*m1*cos(q2 + q3) + L1*L2*m1*cos(q2) + L1*a2*m2*cos(q2) + 2*L2*a1*m1*cos(q3) + L2*a2*m2*cos(q2),                                                                                                                            m1*L2^2 + 2*m1*cos(q3)*L2*a1 + m1*a1^2 + m2*a2^2,                                     a1*m1*(a1 + L2*cos(q3));...
                                                                                                                                                                                                                                                                               a1*m1*(a1 + L1*cos(q2 + q3) + L2*cos(q2 + q3) + L2*cos(q3)),                                                                                                                                                     a1*m1*(a1 + L2*cos(q3)),                                                     a1^2*m1];

                                                                                                                                                               

C = ...
    [                                       - dq2*(L1 + L2)*(a2*m2*sin(q2) + a1*m1*sin(q2 + q3) + L2*m1*sin(q2)) - a1*dq3*m1*(L1*sin(q2 + q3) + L2*sin(q2 + q3) + L2*sin(q3)), - L2^2*dq1*m1*sin(q2) - L2^2*dq2*m1*sin(q2) - L1*a1*dq1*m1*sin(q2 + q3) - L1*a1*dq2*m1*sin(q2 + q3) - L2*a1*dq1*m1*sin(q2 + q3) - L1*a1*dq3*m1*sin(q2 + q3) - L2*a1*dq2*m1*sin(q2 + q3) - L2*a1*dq3*m1*sin(q2 + q3) - L1*L2*dq1*m1*sin(q2) - L1*L2*dq2*m1*sin(q2) - L1*a2*dq1*m2*sin(q2) - L1*a2*dq2*m2*sin(q2) - L2*a2*dq1*m2*sin(q2) - L2*a1*dq3*m1*sin(q3) - L2*a2*dq2*m2*sin(q2), -a1*m1*(dq1 + dq2 + dq3)*(L1*sin(q2 + q3) + L2*sin(q2 + q3) + L2*sin(q3));...
      L2^2*dq1*m1*sin(q2) + L1*a1*dq1*m1*sin(q2 + q3) + L2*a1*dq1*m1*sin(q2 + q3) + L1*L2*dq1*m1*sin(q2) + L1*a2*dq1*m2*sin(q2) + L2*a2*dq1*m2*sin(q2) - L2*a1*dq3*m1*sin(q3),                                                                                                                                                                                                                                                                                                                                                                -L2*a1*dq3*m1*sin(q3),                                       -L2*a1*m1*sin(q3)*(dq1 + dq2 + dq3);...
                                                                                          a1*m1*(L1*dq1*sin(q2 + q3) + L2*dq1*sin(q2 + q3) + L2*dq1*sin(q3) + L2*dq2*sin(q3)),                                                                                                                                                                                                                                                                                                                                                         L2*a1*m1*sin(q3)*(dq1 + dq2),                                                                         0];

G = ...
    [g*m2*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + g*m1*(cos(q1)*(L1 + L2) + L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + g*m1*cos(q1)*(L1 - a1) + g*m2*cos(q1)*(L1 + L2 - a2) + Mh*g*cos(q1)*(L1 + L2);...
                                                                                                                               g*m1*(L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + a2*g*m2*cos(q1 + q2);...
                                                                                                                                                                          a1*g*m1*cos(q1 + q2 + q3)];

B = [-1 0;1 -1;0 1];

tou  = [0;0];
qdot = D \ (-C*[dq1;dq2;dq3]-G + B*tou);
zdot = [dq1;qdot(1);dq2;qdot(2);dq3;qdot(3)];

end

function [value,isterminal,direction] = events_knee_impact(t,xt,flag)
global L1 L2
q1 = xt(1);q3 = xt(5);

ym1=(L1+L2)*sin(q1);

value(1) = q3;
value(2) = ym1-0.5;               %height of hip

isterminal(1) = 1;   % stop the integration
isterminal(2) = 1;   
direction(1) = 1;   % positive direction
direction(2) = -1;   % negative direction
end


function [TE_body,TE_inertial] = energy(x)
    global m1 m2 Mh L1 L2 a1 a2 g;
    q1 = x(1);
    q2 = x(3);
    q3 = x(5);
    dq1 = x(2);
    dq2 = x(4);
    dq3 = x(6);
    y = x(8);
    if(~isnan(dq3))
    KE = (m2*((dq1*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + a2*dq2*cos(q1 + q2))^2 + (dq1*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + a2*dq2*sin(q1 + q2))^2))/2 + (m1*((dq2*(L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)) + dq1*(sin(q1)*(L1 + L2) + L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)) + a1*dq3*sin(q1 + q2 + q3))^2 + (dq2*(L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + dq1*(cos(q1)*(L1 + L2) + L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + a1*dq3*cos(q1 + q2 + q3))^2))/2 + (dq1^2*m1*(L1 - a1)^2)/2 + (dq1^2*m2*(L1 + L2 - a2)^2)/2 + (Mh*dq1^2*(L1 + L2)^2)/2;
    PE_body = g*m1*(sin(q1)*(L1 + L2) + L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)) + g*m2*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + g*m1*sin(q1)*(L1 - a1) + g*m2*sin(q1)*(L1 + L2 - a2) + Mh*g*sin(q1)*(L1 + L2);    
    PE_inertial = g*m1*(sin(q1)*(L1 + L2) + L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)+y) + g*m2*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)+y) + g*m1*(sin(q1)*(L1 - a1)+y) + g*m2*(sin(q1)*(L1 + L2 - a2)+y) + Mh*g*(sin(q1)*(L1 + L2)+y);    
    else
    KE = (m1*((dq1*(cos(q1 + q2)*(L2 + a1) + cos(q1)*(L1 + L2)) + dq2*cos(q1 + q2)*(L2 + a1))^2 + (dq1*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2)) + dq2*sin(q1 + q2)*(L2 + a1))^2))/2 + (m2*((dq1*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + a2*dq2*cos(q1 + q2))^2 + (dq1*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + a2*dq2*sin(q1 + q2))^2))/2 + (dq1^2*m1*(L1 - a1)^2)/2 + (dq1^2*m2*(L1 + L2 - a2)^2)/2 + (Mh*dq1^2*(L1 + L2)^2)/2;
    PE_body = g*m2*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + g*m1*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2)) + g*m1*sin(q1)*(L1 - a1) + g*m2*sin(q1)*(L1 + L2 - a2) + Mh*g*sin(q1)*(L1 + L2);
    PE_inertial = g*m2*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)+y) + g*m1*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2)+y) + g*m1*(sin(q1)*(L1 - a1)+y) + g*m2*(sin(q1)*(L1 + L2 - a2)+y) + Mh*g*(sin(q1)*(L1 + L2)+y);
    end
    TE_body = KE + PE_body;
    TE_inertial = KE + PE_inertial;
end

function zdot=db_manish_2link(t,z)

global m1 m2 Mh L1 L2 a1 a2 g;

q1 = z(1);dq1 = z(2);q2 = z(3);dq2 = z(4);  

D = ...
    [ (m1*(2*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2))^2 + 2*(cos(q1 + q2)*(L2 + a1) + cos(q1)*(L1 + L2))^2))/2 + Mh*(L1 + L2)^2 + m1*(L1 - a1)^2 + (m2*(2*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2))^2 + 2*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2))^2))/2 + m2*(L1 + L2 - a2)^2, L2^2*m1 + a1^2*m1 + a2^2*m2 + L2^2*m1*cos(q2) + 2*L2*a1*m1 + L1*L2*m1*cos(q2) + L1*a1*m1*cos(q2) + L2*a1*m1*cos(q2) + L1*a2*m2*cos(q2) + L2*a2*m2*cos(q2);...
                                                                                                                  L2^2*m1 + a1^2*m1 + a2^2*m2 + L2^2*m1*cos(q2) + 2*L2*a1*m1 + L1*L2*m1*cos(q2) + L1*a1*m1*cos(q2) + L2*a1*m1*cos(q2) + L1*a2*m2*cos(q2) + L2*a2*m2*cos(q2),                                                                                                                  m1*L2^2 + 2*m1*L2*a1 + m1*a1^2 + m2*a2^2];

     
C = ...
    [ -dq2*sin(q2)*(L1 + L2)*(L2*m1 + a1*m1 + a2*m2), -sin(q2)*(dq1 + dq2)*(L1 + L2)*(L2*m1 + a1*m1 + a2*m2);...
       dq1*sin(q2)*(L1 + L2)*(L2*m1 + a1*m1 + a2*m2),                                                      0];

G = ...
    [g*m2*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + g*m1*(cos(q1 + q2)*(L2 + a1) + cos(q1)*(L1 + L2)) + g*m1*cos(q1)*(L1 - a1) + g*m2*cos(q1)*(L1 + L2 - a2) + Mh*g*cos(q1)*(L1 + L2);
                                                                                                                                             g*cos(q1 + q2)*(L2*m1 + a1*m1 + a2*m2)];

qdot = D \ (-C*[dq1;dq2]-G );
zdot = [dq1;qdot(1);dq2;qdot(2)];

end

function [value,isterminal,direction] = events_toe_impact(t,z,flag)
global L1 L2 gamma;
q1 = z(1);
q2 = z(3);
xm1=(L1+L2)*cos(q1);
ym1=(L1+L2)*sin(q1);
xm2=xm1+(L1+L2)*cos(q1+q2);
ym2=ym1+(L1+L2)*sin(q1+q2);
x = (z(4)+z(2)>0); % absolute velocity of link 2 should be negative

value(1) = ym2 + x + xm2*tan(gamma);     % detect height of impact= 0
value(2) = ym1-0.5;               %height of hip

isterminal(1) = 1;   % stop the integration
isterminal(2) = 1;   
direction(1) = -1;
direction(2) = -1;   % negative direction
end

function x_new=impact_grizzly_toe(x)

%impact of 2 link with 5 masses
global m1 m2 Mh L1 L2 a1 a2 g;
q1=x(1); q2=x(3);
De = ...
    [ L1^2*Mh + L2^2*Mh + 2*L1^2*m1 + 2*L1^2*m2 + 2*L2^2*m1 + 2*L2^2*m2 + 2*a1^2*m1 + 2*a2^2*m2 + 2*L1*L2*Mh + 2*L2^2*m1*cos(q2) + 2*L1*L2*m1 + 4*L1*L2*m2 - 2*L1*a1*m1 + 2*L2*a1*m1 - 2*L1*a2*m2 - 2*L2*a2*m2 + 2*L1*L2*m1*cos(q2) + 2*L1*a1*m1*cos(q2) + 2*L2*a1*m1*cos(q2) + 2*L1*a2*m2*cos(q2) + 2*L2*a2*m2*cos(q2), L2^2*m1 + a1^2*m1 + a2^2*m2 + L2^2*m1*cos(q2) + 2*L2*a1*m1 + L1*L2*m1*cos(q2) + L1*a1*m1*cos(q2) + L2*a1*m1*cos(q2) + L1*a2*m2*cos(q2) + L2*a2*m2*cos(q2), - (m1*(2*sin(q1 + q2)*(L2 + a1) + 2*sin(q1)*(L1 + L2)))/2 - (m2*(2*sin(q1)*(L1 + L2) + 2*a2*sin(q1 + q2)))/2 - m1*sin(q1)*(L1 - a1) - m2*sin(q1)*(L1 + L2 - a2) - Mh*sin(q1)*(L1 + L2), (m1*(2*cos(q1 + q2)*(L2 + a1) + 2*cos(q1)*(L1 + L2)))/2 + (m2*(2*cos(q1)*(L1 + L2) + 2*a2*cos(q1 + q2)))/2 + m1*cos(q1)*(L1 - a1) + m2*cos(q1)*(L1 + L2 - a2) + Mh*cos(q1)*(L1 + L2);...
                                                                                                                                                              L2^2*m1 + a1^2*m1 + a2^2*m2 + L2^2*m1*cos(q2) + 2*L2*a1*m1 + L1*L2*m1*cos(q2) + L1*a1*m1*cos(q2) + L2*a1*m1*cos(q2) + L1*a2*m2*cos(q2) + L2*a2*m2*cos(q2),                                                                                                                  m1*L2^2 + 2*m1*L2*a1 + m1*a1^2 + m2*a2^2,                                                                                                                                                  -sin(q1 + q2)*(L2*m1 + a1*m1 + a2*m2),                                                                                                                                                 cos(q1 + q2)*(L2*m1 + a1*m1 + a2*m2);...
                                                                                                                                 - (m1*(2*sin(q1 + q2)*(L2 + a1) + 2*sin(q1)*(L1 + L2)))/2 - (m2*(2*sin(q1)*(L1 + L2) + 2*a2*sin(q1 + q2)))/2 - m1*sin(q1)*(L1 - a1) - m2*sin(q1)*(L1 + L2 - a2) - Mh*sin(q1)*(L1 + L2),                                                                                                                     -sin(q1 + q2)*(L2*m1 + a1*m1 + a2*m2),                                                                                                                                                                       Mh + 2*m1 + 2*m2,                                                                                                                                                                                    0;...
                                                                                                                                   (m1*(2*cos(q1 + q2)*(L2 + a1) + 2*cos(q1)*(L1 + L2)))/2 + (m2*(2*cos(q1)*(L1 + L2) + 2*a2*cos(q1 + q2)))/2 + m1*cos(q1)*(L1 - a1) + m2*cos(q1)*(L1 + L2 - a2) + Mh*cos(q1)*(L1 + L2),                                                                                                                      cos(q1 + q2)*(L2*m1 + a1*m1 + a2*m2),                                                                                                                                                                                      0,                                                                                                                                                                     Mh + 2*m1 + 2*m2];

E = ...
    [ - sin(q1 + q2)*(L1 + L2) - sin(q1)*(L1 + L2), -sin(q1 + q2)*(L1 + L2), 1, 0;...
        cos(q1)*(L1 + L2) + cos(q1 + q2)*(L1 + L2),  cos(q1 + q2)*(L1 + L2), 0, 1];


%tmp_vec=inv([De -E';E zeros(2)])*[De*[x(2);x(4);zeros(2,1)];zeros(2,1)];  
tmp_vec=([De -E';E zeros(2)])\[De*[x(2);x(4);zeros(2,1)];zeros(2,1)];  

%position updates
x_new(1)=q1+q2-pi;
x_new(3)=2*pi - q2;
x_new(5) = 0;
%velocities updates
x_new(2)=tmp_vec(1)+tmp_vec(2);
x_new(4)=-tmp_vec(2);
x_new(6)=0;
%forces at point of impact these values are generally
% used to verify the veracity of the impact model
% for more details refer the grizzle bipedal control book
x_new(7)=tmp_vec(3);
x_new(8)=tmp_vec(4);
                                         
end

function x_new=impact_manish_knee(x)

global m1 m2 Mh L1 L2 a1 a2 g;
syms dq1_p dq2_p

q1 = x(1); q2 = x(3); q3 = x(5);
dq1_n = x(2); dq2_n = x(4); dq3_n = x(6);

sigma_hip_n =...
    L2^2*dq1_n*m1 + L2^2*dq2_n*m1 + a1^2*dq1_n*m1 + a1^2*dq2_n*m1 + a1^2*dq3_n*m1 + a2^2*dq1_n*m2 + a2^2*dq2_n*m2 + L2^2*dq1_n*m1*cos(q2) + L1*a1*dq1_n*m1*cos(q2 + q3) + L2*a1*dq1_n*m1*cos(q2 + q3) + L1*L2*dq1_n*m1*cos(q2) + L1*a2*dq1_n*m2*cos(q2) + 2*L2*a1*dq1_n*m1*cos(q3) + 2*L2*a1*dq2_n*m1*cos(q3) + L2*a2*dq1_n*m2*cos(q2) + L2*a1*dq3_n*m1*cos(q3);
sigma_hip_p =...
    L2^2*dq1_p*m1 + L2^2*dq2_p*m1 + a1^2*dq1_p*m1 + a1^2*dq2_p*m1 + a2^2*dq1_p*m2 + a2^2*dq2_p*m2 + L2^2*dq1_p*m1*cos(q2) + 2*L2*a1*dq1_p*m1 + 2*L2*a1*dq2_p*m1 + L1*L2*dq1_p*m1*cos(q2) + L1*a1*dq1_p*m1*cos(q2) + L2*a1*dq1_p*m1*cos(q2) + L1*a2*dq1_p*m2*cos(q2) + L2*a2*dq1_p*m2*cos(q2);
sigma_pivot_n =...
    m1*((dq2_n*(L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + dq1_n*(cos(q1)*(L1 + L2) + L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + a1*dq3_n*cos(q1 + q2 + q3))*(cos(q1)*(L1 + L2) + L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + (dq2_n*(L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)) + dq1_n*(sin(q1)*(L1 + L2) + L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)) + a1*dq3_n*sin(q1 + q2 + q3))*(sin(q1)*(L1 + L2) + L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3))) + m2*((dq2_n*(L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + dq1_n*(cos(q1)*(L1 + L2) + L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + a1*dq3_n*cos(q1 + q2 + q3))*(cos(q1)*(L1 + L2) + L2*cos(q1 + q2) + a1*cos(q1 + q2 + q3)) + (dq2_n*(L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)) + dq1_n*(sin(q1)*(L1 + L2) + L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3)) + a1*dq3_n*sin(q1 + q2 + q3))*(sin(q1)*(L1 + L2) + L2*sin(q1 + q2) + a1*sin(q1 + q2 + q3))) + Mh*(dq1_n*cos(q1)^2*(L1 + L2)^2 + dq1_n*sin(q1)^2*(L1 + L2)^2) + m1*((dq1_n*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + a2*dq2_n*cos(q1 + q2))*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + (dq1_n*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + a2*dq2_n*sin(q1 + q2))*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2))) + m2*((dq1_n*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + a2*dq2_n*cos(q1 + q2))*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + (dq1_n*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + a2*dq2_n*sin(q1 + q2))*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)));
sigma_pivot_p =...
    Mh*(dq1_p*cos(q1)^2*(L1 + L2)^2 + dq1_p*sin(q1)^2*(L1 + L2)^2) + m1*((dq1_p*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + a2*dq2_p*cos(q1 + q2))*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + (dq1_p*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + a2*dq2_p*sin(q1 + q2))*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2))) + m2*((dq1_p*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + a2*dq2_p*cos(q1 + q2))*(cos(q1)*(L1 + L2) + a2*cos(q1 + q2)) + (dq1_p*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2)) + a2*dq2_p*sin(q1 + q2))*(sin(q1)*(L1 + L2) + a2*sin(q1 + q2))) + m1*((dq1_p*(cos(q1 + q2)*(L2 + a1) + cos(q1)*(L1 + L2)) + dq2_p*cos(q1 + q2)*(L2 + a1))*(cos(q1 + q2)*(L2 + a1) + cos(q1)*(L1 + L2)) + (dq1_p*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2)) + dq2_p*sin(q1 + q2)*(L2 + a1))*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2))) + m2*((dq1_p*(cos(q1 + q2)*(L2 + a1) + cos(q1)*(L1 + L2)) + dq2_p*cos(q1 + q2)*(L2 + a1))*(cos(q1 + q2)*(L2 + a1) + cos(q1)*(L1 + L2)) + (dq1_p*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2)) + dq2_p*sin(q1 + q2)*(L2 + a1))*(sin(q1 + q2)*(L2 + a1) + sin(q1)*(L1 + L2)));

eqn_1 = sigma_hip_n - sigma_hip_p;
eqn_2 = sigma_pivot_n - sigma_pivot_p;
temp = simplify(solve(eqn_1,dq2_p));
eqn_2 = subs(eqn_2,dq2_p,temp);
dq1_p = double(simplify(solve(eqn_2,dq1_p)))
dq2_p = double(simplify(subs(temp,dq1_p)))



x_new = [x(1) dq1_p x(3) dq2_p];

%return x_new of size 1*4
%x_new = x(end,1:4); %have to be replaced by knee impact
end
