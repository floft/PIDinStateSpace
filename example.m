pkg load control;

doPlot = false;

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
if doPlot
    plotResponse(sys_ss, 'Open-Loop Step Response');
    print -dpng -S"700,300" -F"Helvetia:6" image-ol.png
end

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
if doPlot
    plotResponse(sys_cl,
        'Closed-Loop Step Response with LQR controller');
    print -dpng -S"700,300" -F"Helvetia:6" image-lqr.png
end

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
if doPlot
    plotResponseSingle(sys_cl_pid,
        'Closed-Loop Step Response with PID controller');
    print -dpng -S"700,300" -F"Helvetia:6" image-pid.png
end

% Now let's use our new PID in SS form controller
%
% Note: We're using C_siso since with a PID controller you only have one
% output.
Kg = inv(1 + Kd*C_siso*B);

% Check to verify that the issue isn't in removing u from the state space
% equations. It's not. This basically is the same as when using lsim.
%
% To check transfer function in sage:
%    factor(matrix([1,0,0])*~(matrix([[s,0,0],[0,s,0],[0,0,s]])-
%      matrix([[0,1,0],[-101,-20,-200],[1,0,0]]))*matrix([[0],[100],[-1]]))
%
% Compare:
%    feedback(pid(Kp,Ki,Kd)*sys_tf,1)
%    tf(sys_ss_pid)
if doPlot && false
    % The open-loop A, B, C, and D
    Apid_ol = [A zeros(size(A,1),1); C_siso zeros(1,1)];
    Bpid_ol = [B zeros(size(B,1),1); 0 -1];
    Cpid_ol = [C_siso 0; C_siso 0; zeros(1,size(C,2)) 1; Kg*C_siso*A 0];
    Dpid_ol = [0 0; 0 -1; 0 0; 0 0];

    % Discretize to have our own lsim-like simulation
    f = 100;
    T = 1/f;
    sys_d = c2d(ss(Apid_ol,Bpid_ol,Cpid_ol,Dpid_ol), T, 'zoh');

    N = 4*f;
    state = zeros(size(C_siso,2)+1, 2);
    output = zeros(N, size(Dpid_ol,1));

    % Constant set point
    s = 0.2;
    input = zeros(N,2);

    for i = 2:N
        input(i,:) = [-[0 Kp Ki Kd]*output(i-1,:)'; s]';
        state(:,1) = sys_d.a*state(:,2) + sys_d.b*input(i,:)';
        output(i,:) = sys_d.c*state(:,2) + sys_d.d*input(i,:)';
        state(:,2) = state(:,1);
    end

    t = 0:T:(size(output,1)-1)/f;
    figure;
    plot(t,output(:,1));
    ylabel('x');
    title('PID in SS - without lsim');
end

Apid = [A-B*(Kp*C_siso+Kd*Kg*C_siso*A) -B*Ki; C_siso 0];
Bpid = [B*Kp; -1];
Cpid = [C_siso 0];
Dpid = 0;

states_pid = {'x' 'x_dot' 'z'};
inputs_pid = {'s'};
sys_ss_pid = ss(Apid,Bpid,Cpid,Dpid,
    'statename',states_pid,
    'inputname',inputs_pid,
    'outputname',outputs_siso);
if doPlot
    plotResponseSingle(sys_ss_pid,
        'Closed-Loop Step Response with PID controller in SS form');
    print -dpng -S"700,300" -F"Helvetia:6" image-pid-ss.png
end

% These should be the same
feedback(pid(Kp,Ki,Kd)*sys_tf,1)
tf(sys_ss_pid)

% Just use a PI controller, which does look the same
Kp = 5;
Ki = 10;
Kd = 0;

pid_controller = pid(Kp,Ki,Kd);
sys_cl_pid = feedback(pid_controller*sys_tf,1);
if doPlot
    plotResponseSingle(sys_cl_pid,
        'Closed-Loop Step Response with PI controller');
    print -dpng -S"700,300" -F"Helvetia:6" image-pi.png
end

Kg = inv(1 + Kd*C_siso*B);
Api = [A-B*(Kp*C_siso+Kd*Kg*C_siso*A) -B*Ki; C_siso 0];
Bpi = [B*Kp; -1];
Cpi = [C_siso 0];
Dpi = 0;

sys_ss_pi = ss(Api,Bpi,Cpi,Dpi,
    'statename',states_pid,
    'inputname',inputs_pid,
    'outputname',outputs_siso);
if doPlot
    plotResponseSingle(sys_ss_pi,
        'Closed-Loop Step Response with PI controller in SS form');
    print -dpng -S"700,300" -F"Helvetia:6" image-pi-ss.png
end
