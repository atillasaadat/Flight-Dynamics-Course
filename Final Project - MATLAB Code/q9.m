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

dt = 0.001; %define step size for simulation time
maxT = 10.0;
t = 0:dt:maxT; %define simulation time domain
u = deg2rad(10)*ones(size(t,2),1); %define step input wih amplitude 1.6 for all time steps

%% 

%Rrho,rho_rww,rho_rvv,intP
[yroll,troll,xroll,Kroll,Lroll,Rroll] = lq_servo(0.1,1,0.01,100,Aroll,Broll,Croll,D,u,t);
KIroll = Kroll(4);
Kroll = Kroll(1:3);
Rroll = [Rroll; 0; 0];
[ypitch,tpitch,xpitch,Kpitch,Lpitch,Rpitch] = lq_servo(0.15,10,0.1,20,Apitch,Bpitch,Cpitch,D,u,t);
KIpitch = Kpitch(4);
Kpitch = Kpitch(1:3);
Rpitch = [Rpitch; 0; 0];

function [y,t,x,Kbar,Lbar,Raug] = lq_servo(Rrho,rho_rww,rho_rvv,intP,Aq,Bq,Cq,D,u,t)
    Aaug = [Aq zeros(3,1); -Cq zeros(1,1)];
    Baug = [Bq; 0];
    
    Qaug = [Cq'*Cq zeros(3,1); zeros(1,3) intP];
    Raug = Rrho*eye(size(Bq,2));
    Kbar = lqr(Aaug,Baug,Qaug,Raug);
    
    Rww = rho_rww*eye(size(Bq,2));
    Rvv = rho_rvv*eye(size(Bq,2));
    Lbar = lqr(Aq',Cq',Bq*Rww*Bq',Rvv)';
    %compensators
    ac1=[Aq-Lbar*Cq-Bq*Kbar(1:size(Aq,1)) -Bq*Kbar(end);zeros(3,1)' 0];
    Br=[Lbar*0;1];
    
    acl = [Aq -Bq*Kbar; [Lbar;-1]*Cq ac1];
    bcl = [Bq*0;Br];
    ccl = [Cq Cq*0 0];
    ss_dofb = ss(acl,bcl,ccl,D);
    [y,t,x] = lsim(ss_dofb,u,t);
end
