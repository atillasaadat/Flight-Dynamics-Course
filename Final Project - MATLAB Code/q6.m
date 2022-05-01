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
%L_p = mean([0.02,0.2]);
L_p = 0.02;
%M_q = mean([0.1,0.9]);
M_q = 0.1;

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

dt = 0.001; %define step size for simulation time
maxT = 10.0;
t = 0:dt:maxT; %define simulation time domain
u = deg2rad(10)*ones(size(t,2),1); %define step input wih amplitude 1.6 for all time steps

ss_roll = ss(Aroll,Broll,Croll,D);
ss_pitch = ss(Apitch,Bpitch,Cpitch,D);
[yroll,troll,xroll] = lsim(ss_roll,u,t);
[ypitch,tpitch,xpitch] = lsim(ss_pitch,u,t);
%% 
rho_rvv_roll = 0.4;
rho_rww_roll = 0.1;
Rww_roll = rho_rww_roll*eye(size(Broll,2));
Rvv_roll = rho_rvv_roll*eye(size(Broll,2));
Lroll = lqr(Aroll',Croll',Broll*Rww_roll*Broll',Rvv_roll)';
xhat = [0;0;0];
XXhat = xhat;
y_norm = [0];
y_hat = [0];
for i=1:size(troll)-1
    uu = u(i);
    y = Croll*(xroll(i,:))'+D*uu;
    y_norm = [y_norm y];
    yhat = Croll*xhat+D*uu;
    y_hat = [y_hat yhat];    
    xhat = (Aroll*xhat+Broll*uu+Lroll*(y-yhat)')*dt + xhat;
    XXhat = [XXhat,xhat];
end
XXhat = XXhat';
figure
plot(troll,yroll,troll,y_hat,'--','LineWidth',2);
title('LQE - Roll Angle','Interpreter','latex');
legend({'Raw SS Response','LQE Estimation'},'Location','southeast','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'FontWeight','bold','Interpreter','latex');
ylabel('Roll Angle $\phi$ [deg]','FontSize',14,'FontWeight','bold','Interpreter','latex');
grid on
grid minor
figure
plot(troll,xroll-XXhat,'--','LineWidth',2);
title('LQE - Roll States Estimator Error','Interpreter','latex');
roll_legend = {'$\dot{e}_{x_1}$','$\dot{e}_{x_2}$','$\dot{e}_{x_3}$'};
legend(roll_legend,'Location','southeast','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'FontWeight','bold','Interpreter','latex');
ylabel('Value','FontSize',14,'FontWeight','bold','Interpreter','latex');
grid on
grid minor
%% 
rho_rvv_pitch = 0.4;
rho_rww_pitch = 1;
Rww_pitch = rho_rww_pitch*eye(size(Bpitch,2));
Rvv_pitch = rho_rvv_pitch*eye(size(Bpitch,2));
Lpitch = lqr(Apitch',Cpitch',Bpitch*Rww_pitch*Bpitch',Rvv_pitch)';
xhat = [0;0;0];
XXhat = xhat;
y_norm = [0];
y_hat = [0];
for i=1:size(tpitch)-1
    uu = u(i);
    y = Cpitch*(xpitch(i,:))'+D*uu;
    y_norm = [y_norm y];
    yhat = Cpitch*xhat+D*uu;
    y_hat = [y_hat yhat];    
    xhat = (Apitch*xhat+Bpitch*uu+Lpitch*(y-yhat)')*dt + xhat;
    XXhat = [XXhat,xhat];
end
XXhat = XXhat';
figure
plot(tpitch,ypitch,tpitch,y_hat,'--','LineWidth',2);
title('LQE - Pitch Angle','Interpreter','latex');
legend({'Raw SS Response','LQE Estimation'},'Location','southeast','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'FontWeight','bold','Interpreter','latex');
ylabel('Pitch Angle $\theta$ [deg]','FontSize',14,'FontWeight','bold','Interpreter','latex');
grid on
grid minor
figure
plot(troll,xpitch-XXhat,'--','LineWidth',2);
title('LQE - Pitch States Estimator Error','Interpreter','latex');
roll_legend = {'$\dot{e}_{x_1}$','$\dot{e}_{x_2}$','$\dot{e}_{x_3}$'};
legend(roll_legend,'Location','southeast','Interpreter','latex');
xlabel('Time [s]','FontSize',14,'FontWeight','bold','Interpreter','latex');
ylabel('Value','FontSize',14,'FontWeight','bold','Interpreter','latex');
grid on
grid minor
