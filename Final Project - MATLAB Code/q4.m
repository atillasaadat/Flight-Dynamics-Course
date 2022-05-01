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

sys = ss(A,B,C,D);

%determine stability
stab = eig(Aroll);
disp(stab);
%determine observability
roll_Ob = obsv(Aroll,Croll);
roll_unob = length(Aroll)-rank(roll_Ob);

%determine controllability
roll_Co = ctrb(Aroll,Broll);
roll_unco = length(Aroll) - rank(roll_Co);

pitch_Ob = obsv(Apitch,Cpitch);
pitch_unob = length(Apitch)-rank(pitch_Ob);

%determine controllability
pitch_Co = ctrb(Apitch,Bpitch);
pitch_unco = length(Apitch) - rank(pitch_Co);
%% 

Qroll = 1*(Croll'*Croll);
Rroll = 0.02*eye(size(Broll,2));

[Kroll,Proll,lambda_cl_roll] = lqr(Aroll,Broll,Qroll,Rroll);
ss_orig = ss(Aroll,Broll,Croll,D);
ss_cl = ss((Aroll-Broll*Kroll),Broll,Croll,D);

Nbar = -inv(Croll*(Aroll-Broll*Kroll)^-1*Broll);
ss_clNbar = ss((Aroll-Broll*Kroll),Broll*Nbar,Croll,D);
%knAroll = (Aroll-Broll*Kroll);
%knBroll = Broll*Nbar;
%knAroll = Aroll
%knBroll = Broll

dt = 0.001; %define step size for simulation time
maxT = 5.0;
t = 0:dt:maxT; %define simulation time domain
u = deg2rad(20)*ones(size(t,2),1); %define step input wih amplitude 1.6 for all time steps

[y1,t1,x1] = lsim(ss_orig,u,t); 
[y2,t2,x2] = lsim(ss_cl,u,t); 
[y3,t3,x3] = lsim(ss_clNbar,u,t);

plot(t1,rad2deg(y1),t2,rad2deg(y2),t3,rad2deg(y3));
%plot(t1,rad2deg(y1),t2,rad2deg(y2),tlqe_roll,rad2deg(ylqe_roll));
xlim([0 maxT]);
title('LQR - Roll Controller','Interpreter','latex');
legend({'Raw SS','LQR','LQR w/ NBar'},'Location','southeast','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'FontWeight','bold','Interpreter','latex');
ylabel('Roll Angle $\phi$ [deg]','FontSize',14,'FontWeight','bold','Interpreter','latex');
grid on
grid minor

figure
uRoll = -Kroll.*x3;
Vcyc_roll = uRoll(:,3);
plot(t1,Vcyc_roll);
title('$V_{cyc}$ Voltage Timeseries','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'FontWeight','bold','Interpreter','latex');
ylabel('$V_{cyc} [V]$','FontSize',14,'FontWeight','bold','Interpreter','latex');
grid on
grid minor
