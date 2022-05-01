%MECH-4672 Flight Dynamics and Control of UAVs - Course Project
%Atilla Saadat - 104411786
clc;clear;close all;

m = 1.15;
M = 3.57;
l_boom = 0.66;
l_psi = 0.004;
l_theta = 0.014;
l_h = 0.177;
I_xx = 0.036;
I_yy = 0.93;
K_t = 4.25e-3;
theta_rest = deg2rad(-25);
g = 9.81;
L_p = mean([0.02,0.2]);
M_q = mean([0.1,0.9]);

I_zz = 0.93;
K_D = 0;
omega_coll_0 = M*g*l_theta*sin(theta_rest)/(K_t*l_boom);
K_v = 0.0125*omega_coll_0*K_t*l_boom;

A = [0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 0 0 1 0 0;
    -m*g*l_psi/I_xx 0 0 -L_p/I_xx 0 -K_v*l_h/I_xx K_t*l_h/I_xx 0;
    0 -M*g*l_theta*cos(theta_rest)/I_yy 0 0 -M_q/I_yy -K_v*l_boom/I_yy 0 K_t*l_boom/I_yy;
    K_t*l_boom*omega_coll_0/I_zz 0 0 0 0 -K_D*l_boom/I_zz 0 0;
    0 0 0 0 0 0 -6 0;
    0 0 0 0 0 0 0 -6;
    ];

Aroll = [
    0 1 0;
    -m*g*l_psi/I_xx -L_p/I_xx K_t*l_h/I_xx;
    0 0 -6
    ];
        
Apitch = [
    0 1 0;
    -M*g*l_theta*cos(theta_rest)/I_yy -M_q/I_yy K_t*l_boom/I_yy;
    0 0 -6
    ];

Atravel = [0 1; -K_D*l_boom/I_zz 0;];

B = [0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     0   0;
     780 0;
     0 540];
 
Broll = [0;0;780];
Bpitch = [0;0;540];
Btravel = [0;K_t*l_boom*omega_coll_0/I_zz];

C = [1 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0];
Croll = [1 0 0];
Cpitch = [1 0 0];
  
D = 0;

%% LQ Servo - Pitch

%A2pitch = horzcat([Apitch; -Cpitch],zeros(4,1));
A2pitch = [Apitch zeros(3,1); -Cpitch zeros(1,1)];
B2pitch = [Bpitch; 0];
C2pitch = [Cpitch 0];

intP = 300; %integration penalty

Rrho_pitch = 0.2;

dt = 0.001; %define step size for simulation time
maxT = 10.0;
t = 0:dt:maxT; %define simulation time domain
u = deg2rad(20)*ones(size(t,2),1); %define step input wih amplitude 1.6 for all time steps

[y1,t1,x1,Kpitch] = lq_servo(0.1,20,Apitch,Bpitch,Cpitch,D,u,t);
[y2,t2,x2,Kpitch] = lq_servo(0.1,300,Apitch,Bpitch,Cpitch,D,u,t);
[y3,t3,x3,Kpitch] = lq_servo(0.05,20,Apitch,Bpitch,Cpitch,D,u,t);
[y4,t4,x4,Kpitch] = lq_servo(0.05,300,Apitch,Bpitch,Cpitch,D,u,t);

figure
plot(t1,[rad2deg(y1) rad2deg(y2) rad2deg(y3) rad2deg(y4)],'--','LineWidth',2);
xlim([0 maxT]);
title('LQ Servo - Pitch Controller','Interpreter','latex');
legend({'$\rho_{pitch} = 0.1, Int. Penalty = 20$',...
'$\rho_{pitch} = 0.1, Int. Penalty = 300$',...
'$\rho_{pitch} = 0.05, Int. Penalty = 20$',...
'$\rho_{pitch} = 0.05, Int. Penalty = 300$'},...
'Location','southeast','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'Interpreter','latex');
ylabel('Pitch Angle $\theta$ [deg]','FontSize',14,'Interpreter','latex');
grid on;
grid minor

figure
uRoll = -Kpitch.*x4;
Vcoll_pitch = uRoll(:,3);
plot(t1,Vcoll_pitch);
title('$V_{coll}$ Voltage Timeseries','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'FontWeight','bold','Interpreter','latex');
ylabel('$V_{cyc} [V]$','FontSize',14,'FontWeight','bold','Interpreter','latex');
grid on
grid minor

function [y,t,x,Kbar] = lq_servo(Rrho,intP,Apitch,Bpitch,Cpitch,D,u,t)
    
    Aaug = [Apitch zeros(3,1); -Cpitch zeros(1,1)];
    Baug = [Bpitch; 0];
    Caug = [Cpitch 0];
    
    QQ = Cpitch'*Cpitch;
    Qpitch = [QQ zeros(3,1); zeros(1,3) intP];
    RRpitch = Rrho*eye(size(Bpitch,2));
    
    Kbar = lqr(Aaug,Baug,Qpitch,RRpitch);
    lqs_sys = ss(Aaug-Baug*Kbar,[0 0 0 1]',Caug,D);
    [y,t,x] = lsim(lqs_sys,u,t);
end
