function Simulation (xa,ya,za,xb,yb,zb)
%% *** Robot (kinematic) model parameters *** 
clear all; 
close all;
l0 = 0.0;  %% in cm 
l1 = 0.0;  
l2 = 5.0;
l3 = 0.0;
l4 = 10.0;
l5 = 10.0;

%% *** sampling period *** 
%% *** for the robot motion, kinematic simulation: 
dt = 0.01; 

%% *** Create (or load from file) reference signals *** 
%% *** DESIRED MOTION PROFILE - TASK SPACE *** 
Tf = 10.0; 	% 10sec duration of motion 
D = 1.0;    
t = 0:dt:2*Tf;

%xd0,td0,yd1: initial/final end-point position --> desired task-space trajectory  
xd0 = 10.0;	% Ax
xd1 = 0.0;  % Bx
yd0 = 0.0;  % Ay
yd1 = 10.0; % By
zd0 = 10; % Az
zd1 = 10; % Bz
g = 10 ; 

% Example of desired trajectory : linear segment (x0,y0)-->(x1,y1); Time duration: Tf; 
disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...'); %% 
disp(' ');   

% Creating trajectory
xAB = Create_Trajectory(xd0, xd1, 0, Tf, D);
xBA = Create_Trajectory(xd1, xd0, Tf, 2*Tf, D);

yAB = Create_Trajectory(yd0, yd1, 0, Tf, D);
yBA = Create_Trajectory(yd1, yd0, Tf, 2*Tf, D);

zAB = Create_Trajectory(zd0, zd1, 0, Tf, D);
zBA = Create_Trajectory(zd1, zd0, Tf, 2*Tf, D);

xd = zeros(1, length(t));
yd = zeros(1, length(t));
zd = zeros(1, length(t));
vx = zeros(1, length(t));
vy = zeros(1, length(t));
vz = zeros(1, length(t));

% position sampling
for i = 1:1:length(0:dt:D)
    tb = i*dt;
    xd(i) = polyval(double(xAB(1:6)),tb);
    yd(i) = polyval(double(yAB(1:6)),tb);
    zd(i) = polyval(double(zAB(1:6)),tb);
    vx(i) = polyval(polyder(double(xAB(1:6))),tb);
    vy(i) = polyval(polyder(double(yAB(1:6))),tb);
    vz(i) = polyval(polyder(double(zAB(1:6))),tb);
end

for i = length(0:dt:D)+1:1:length(0:dt:Tf-D)
    tb = i*dt;
    xd(i) = polyval(double(xAB(7:8)),tb);
    yd(i) = polyval(double(yAB(7:8)),tb);
    zd(i) = polyval(double(zAB(7:8)),tb);
    vx(i) = polyval(polyder(double(xAB(7:8))),tb);
    vy(i) = polyval(polyder(double(yAB(7:8))),tb);
    vz(i) = polyval(polyder(double(zAB(7:8))),tb);
end

for i = length(0:dt:Tf-D)+1:1:length(0:dt:Tf)
    tb = i*dt;
    xd(i) = polyval(double(xAB(9:14)),tb);
    yd(i) = polyval(double(yAB(9:14)),tb);
    zd(i) = polyval(double(zAB(9:14)),tb);
    vx(i) = polyval(polyder(double(xAB(9:14))),tb);
    vy(i) = polyval(polyder(double(yAB(9:14))),tb);
    vz(i) = polyval(polyder(double(zAB(9:14))),tb);
end

for i = length(0:dt:Tf)+1:1:length(0:dt:Tf+D)
    tb = i*dt;
    xd(i) = polyval(double(xBA(1:6)),tb);
    yd(i) = polyval(double(yBA(1:6)),tb);
    zd(i) = polyval(double(zBA(1:6)),tb);
    vx(i) = polyval(polyder(double(xBA(1:6))),tb);
    vy(i) = polyval(polyder(double(yBA(1:6))),tb);
    vz(i) = polyval(polyder(double(zBA(1:6))),tb);
end

for i = length(0:dt:Tf+D)+1:1:length(0:dt:2*Tf-D)
    tb = i*dt;
    xd(i) = polyval(double(xBA(7:8)),tb);
    yd(i) = polyval(double(yBA(7:8)),tb);
    zd(i) = polyval(double(zBA(7:8)),tb);
    vx(i) = polyval(polyder(double(xBA(7:8))),tb);
    vy(i) = polyval(polyder(double(yBA(7:8))),tb);
    vz(i) = polyval(polyder(double(zBA(7:8))),tb);
end

for i = length(0:dt:2*Tf-D)+1:1:length(0:dt:2*Tf)
    tb = i*dt;
    xd(i) = polyval(double(xBA(9:14)),tb);
    yd(i) = polyval(double(yBA(9:14)),tb);
    zd(i) = polyval(double(zBA(9:14)),tb);
    vx(i) = polyval(polyder(double(xBA(9:14))),tb);
    vy(i) = polyval(polyder(double(yBA(9:14))),tb);
    vz(i) = polyval(polyder(double(zBA(9:14))),tb);
end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% ****** KINEMATIC SIMULATION - Main loop ****** 
disp('Kinematic Simulation ...'); %% 
disp(' '); %%  

%% ***** INVESRE KINEMATICS  -->  DESIRED MOTION - JOINT SPACE ***** 
%% compute the reference joint-motion vectors: 
%% {qd(k,i), i=1,...,n (num of degrees of freedom), with k=1,..., kmax,} 
%% and reference joint (angular) velocities {qd_1(k,i)} 
%% and reference joint (angular) velocities {qd_1(k,i)} 
rd3 = (xd(:)-l1).^2 + yd(:).^2 + (zd(:)+l0).^2; 
qd(:,3) = - acos((rd3(:) - l1^2 - l2^2 - l4^2 - l5^2)./(2*l4*l5)); 	

s3 = real(sin(qd(:,3)));
c3 = real(cos(qd(:,3))); 

rd2 = xd(:).^2 + (zd(:)-l0).^2 - l1^2; 
qd(:,2) = -atan2(l5*s3(:), l4 + l5*c3(:)) + asin(yd(:) ./ sqrt(l4^2 + 2*l4*l5*c3(:) + l5^2))

c2 = real(cos(qd(:,2))); 
s2 = real(sin(qd(:,2)));
c23 = real(cos(qd(:,3) + qd(:,2)));
s23 = real(sin(qd(:,3) + qd(:,2)));

qd(:,1) = asin(l2 ./ (sqrt((xd(:) - l1).^2 + (zd(:) + l0).^2))) - asin((xd(:) - l1) ./ (sqrt((xd(:) - l1).^2 + (zd(:) + l0).^2)))

c1 = cos(qd(:,1)); 
s1 = sin(qd(:,1)); 

%% ***** FORWARD KINEMATICS  JOINT MOTION -->  CARTESIAN POSITIONS ***** 
%%(xd1, yd1, zd1) : cartesian position of the 1st link's local reference frame 
xd1 = l1*ones(1, length(t));   
yd1 = zeros(1, length(t));
zd1 = -l0*ones(1, length(t));
%%(xd2, yd2, zd2) : cartesian position of the 2nd link's local reference frame 
xd2 = l1*ones(1, length(t));  
yd2 = zeros(1, length(t));
zd2 = -l0*ones(1, length(t)); 
%%(xd3, yd3, zd3) : cartesian position of the 3rd link's local reference frame
xd3 = l4*s1(:).*c2(:) + l2*c1(:) + l1  
yd3 = l4*s2(:); 
zd3 = -l4*c1(:).*c2(:) + l2*s1(:) - l0 
%%(xdE, ydE, zdE) : cartesian position of endeffector's local reference frame
xdE = l5*s1(:).*c23(:) + l4*s1(:).*c2(:) + l2*c1(:) + l1  
ydE = l4*s2(:) + l5*s23(:); 
zdE = -l5*c1(:).*c23(:) - l4*c1(:).*c2(:) +l2*s1(:) - l0 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% ***** INVERSE DIFFERENTIAL KINEMATICS ***** 
dq1 = (c1(:).*vx(:))./(l4*c2(:)+l5*c23(:)) + ((s1(:).*vz(:)) ./ (l4*c2(:) + l5*c23(:)));
dq2 = ((c23(:).*(l5*s1(:).*c23(:) + l4*s1(:).*c2(:) + l2*c1(:))).*vx(:))./(l4*s3(:).*(l4*c2(:)+l5*c23(:))) + (s23(:).*vy(:))./(l4*s3(:)) - (c23(:).*(l5*c1(:).*c23(:) + l4*c1(:).*c2(:) + l2*s1(:)).*vz(:))./(l4*s3(:).*(l4*c2(:)+l5*c23(:)));
dq3 = -((l5*s1(:).*c23(:) + l4*s1(:).*c2(:) + l2*c1(:)).*vx(:))./(l4*l5*s3(:)) - ((l4*l5*s3(:)+l5^2*s23(:).*c23(:)+l4^2*s2(:).*c2(:)).*vy(:) ./ (l4*l5*s3(:).*(l4*c2(:)+l5*c23(:)))) + (((l5*c1(:).*c23(:) + l4*c1(:).*c2(:) + l2*s1(:)).*vz(:)) ./ (l4*l5*s3(:))) ;

%% *** SAVE and PLOT output data *** %%** use functions plot(...)  
save;  %% --> save data to 'matlab.mat' file   

fig1 = figure;  
subplot(3,1,1);
plot(t,xd); 
ylabel('x_d (cm)'); 
xlabel('time t (sec)');  
axis([0 20 -20 20]);

subplot(3,1,2); 
plot(t,yd); 
ylabel('y_d (cm)'); 
xlabel('time t (sec)');  
axis([0 20 0 20]);

subplot(3,1,3); 
plot(t,zd); 
ylabel('z_d (cm)'); 
xlabel('time t (sec)'); 
axis([0 20 0 20]);

%% sgtitle('p_{Ex} - p_{Ey} - p_{Ez}');

fig2 = figure;  
subplot(3,1,1); 
plot(t,vx); 
ylabel('v_x (cm/sec)'); 
xlabel('time t (sec)');  
axis([0 20 -2 2]);

subplot(3,1,2); 
plot(t,vy); 
ylabel('v_y (cm/sec)'); 
xlabel('time t (sec)');  
axis([0 20 -2 2]);

subplot(3,1,3); 
plot(t,vz); 
ylabel('v_z (cm/sec)'); 
xlabel('time t (sec)');   

%% sgtitle('v_{Ex} - v_{Ey} - v_{Ez}');

fig3 = figure; 
subplot(3,1,1); 
plot(t,qd(:,1)); 
ylabel('q_1 (rad)'); 
xlabel('time t (sec)'); 

subplot(3,1,2); 
plot(t,qd(:,2)); 
ylabel('q_2 (rad)'); 
xlabel('time t (sec)');  

subplot(3,1,3); 
plot(t,qd(:,3)); 
ylabel('q_3 (rad)'); 
xlabel('time t (sec)');  

%% sgtitle('q_1 - q_2 - q_3');

fig4 = figure; 
subplot(3,1,1); 
plot(t,dq1); 
ylabel('dq_1/dt (rad/sec)'); 
xlabel('time t (sec)'); 

subplot(3,1,2); 
plot(t,dq2); 
ylabel('dq_2/dt (rad/sec)'); 
xlabel('time t (sec)');  

subplot(3,1,3); 
plot(t,dq3); 
ylabel('dq_3/dt (rad/sec)'); 
xlabel('time t (sec)'); 

%% sgtitle('dq_1/dt - dq_2/dt - dq_3/dt');

%%*** stick diagram --> animate robot motion ... (**optional**) 
%% within a for (or while) loop, use periodic plot(...) functions to draw the geometry (current pos)  
%% of the robot, and thus animate its motion ...  

fig5 = figure; 
axis([-10 35 -10 35 -10 35]) %%set xyz plot axes (caution: square axes, i.e. dx=dy) 
axis on 
hold on 
xlabel('x (cm)'); 
ylabel('y (cm)'); 
zlabel('z (cm)'); 
plot3(xd,yd,zd,'rs'); 
dtk=100; %% plot robot position every dtk samples, to animate its motion 
plot3([0],[0],[0],'o'); 
kmax=2*Tf/dt + 1; 
for tk=1:dtk:kmax,    %%% 	
   pause(0.1);	%% pause motion to view successive robot configurations    
   plot3([0,xd1(tk)],[0,yd1(tk)],[0,zd1(tk)]);					
   plot3([xd1(tk)],[yd1(tk)],[zd1(tk)],'o');        
   plot3([xd1(tk),xd2(tk)],[yd1(tk),yd2(tk)],[zd1(tk),zd2(tk)]);	
   plot3([xd2(tk)],[yd2(tk)],[zd2(tk)],'o');    
   plot3([xd2(tk),xd3(tk)],[yd2(tk),yd3(tk)],[zd2(tk),zd3(tk)]);
   plot3([xd3(tk)],[yd3(tk)],[zd3(tk)],'o');  
   plot3([xd3(tk),xdE(tk)],[yd3(tk),ydE(tk)],[zd3(tk),zdE(tk)]);
   plot3([xdE(tk)],[ydE(tk)],[zdE(tk)],'y*');   
   plot3([xdE(tk)],[ydE(tk)],[zdE(tk)],'g+');  
end       

title('Robot Motion')

return 
