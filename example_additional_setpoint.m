pkg load control;

% Plot both x and x_dot
function plotResponse(sys, plotTitle)
figure;
t = 0:0.01:4;
r = 0.2*ones(size(t));
[y,t,x]=lsim(sys,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','x');
set(get(AX(2),'Ylabel'),'String','x_dot');
title(plotTitle)
end

% Plot only x
function plotResponseSingle(sys, plotTitle)
figure;
t = 0:0.01:4;
r = 0.2*ones(size(t));
[y,t,x]=lsim(sys,r,t);
plot(t,y,'-b');
ylabel('x');
title(plotTitle)
end

% State space for simple second-order spring-mass equation
k = 1;
m = 1;

A = [0 1; -k/m 0];
B = [0; 1/m];
C = eye(2);
D = 0;

states = {'x' 'x_dot'};
inputs = {'F'};
outputs = {'x'; 'x_dot'};
sys_ss = ss(A,B,C,D,
    'statename',states,
    'inputname',inputs,
    'outputname',outputs);
plotResponse(sys_ss, 'Open-Loop Step Response');
print -dpng -S"700,300" -F"Helvetia:6" image-ol.png

% Controller using LQR
Q = C'*C;
Q(1,1) = 10;
R = 0.01;
K = lqr(A,B,Q,R);

% Correct position error
Cn = [1 0];
sys_nbar = ss(A,B,Cn,0);
Nbar = rscale(sys_nbar,K);

Ac = [(A-B*K)];
Bc = [B*Nbar];
Cc = [C];
Dc = [D];

sys_cl = ss(Ac,Bc,Cc,Dc,
    'statename',states,
    'inputname',inputs,
    'outputname',outputs);
plotResponse(sys_cl,
    'Closed-Loop Step Response with LQR controller');
print -dpng -S"700,300" -F"Helvetia:6" image-lqr.png

% Use a PID controller
Kp = 100;
Ki = 200;
Kd = 20;

% We need to have SISO, so redefine C to only give us x out
C_siso = [1 0];
outputs_siso = {'x'};
sys_ss_siso = ss(A,B,C_siso,D,
    'statename',states,
    'inputname',inputs,
    'outputname',outputs_siso);
sys_tf = tf(sys_ss_siso);

pid_controller = pid(Kp,Ki,Kd);
sys_cl_pid = feedback(pid_controller*sys_tf,1);
plotResponseSingle(sys_cl_pid,
    'Closed-Loop Step Response with PID controller');
print -dpng -S"700,300" -F"Helvetia:6" image-pid.png

% Now let's use our new PID in SS form controller
%
% Note 1: We're using C_siso since with a PID controller you only have one
% output.
%
% Note 2: bf = [0;0] (and thus Kg = 1) since normally feedback does not couple
% directly and instantaneously into the output, but then it doesn't work
bf = B;
bs = B;
Kg = inv(1 + Kd*C_siso*bf);

% u = [v;s]
Apid = [A zeros(size(A,1),1); C_siso zeros(1,1)];
Bpid = [bf bs; 0 -1];
Cpid = [C_siso 0; C_siso 0; zeros(1,size(C,2)) 1; Kg*C_siso*A 0];
Dpid = [0 0; 0 -1; 0 0; 0 Kg*C_siso*bs];

% u = s, allowing us to run this nicely in lsim
Apid_cl = [A-bf*(Kp*C_siso+Kd*Kg*C_siso*A) -bf*Ki; C_siso 0];
Bpid_cl = [bf*(Kp-Kd*Kg*C_siso*bs)+bs; -1];
Cpid_cl = [C_siso 0];
Dpid_cl = 0;

states_pid = {'x' 'x_dot' 'z'};
inputs_pid = {'s'};
sys_ss_pid = ss(Apid_cl,Bpid_cl,Cpid_cl,Dpid_cl,
    'statename',states_pid,
    'inputname',inputs_pid,
    'outputname',outputs_siso);
plotResponseSingle(sys_ss_pid,
    'Closed-Loop Step Response with PID controller in SS form');
print -dpng -S"700,300" -F"Helvetia:6" image-pid-ss.png
