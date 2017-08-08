%% Week 8, Application of RK to Lorenz 
% Lorenz model
% xdot = beta(y - x)
% ydot = x(rho - z) - y
% zdot = xy - sigma * z

clear all
% Input Lorenz' original parameters that lead to chaos
sigma = 10; beta = 8/3; rho = 28;

%initial cond
y0 = [-8; 8; 27];

%computing trajectory
dt = 0.01;
tspan = [0:dt:15];

Y(:, 1) = y0; %list of trajectory values x,y,z per time step.
y_in = y0;
for i = 1:tspan(end)/dt %this gives the number of steps in tspan
    time = i*dt;
    y_out = rk4_singleStep(@(t,y) lorenz(t, y, sigma, beta, rho), dt, time, y_in);
    Y = [Y y_out];
    y_in = y_out;
end
plot3(Y(1,:), Y(2,:), Y(3,:), 'b');
hold on
% can of course do this with ode 45
[t, y] = ode45(@(t,y) lorenz(t, y, sigma, beta, rho), tspan, y0);
plot3(y(:, 1), y(:,2), y(:,3), 'r')
legend('RK4ss', 'ode45')
%% Second lecture vectorized time step integrators
clear all; clc; close all;
% Input Lorenz' original parameters that lead to chaos
sigma = 10; beta = 8/3; rho = 28;

% Initial cond is a large cube of points, like a 3d mesh grid?
xvec = -20:2:20; yvec = -20:2:20; zvec = -20:2:20;

%called it
[x0, y0, z0] = meshgrid(xvec, yvec, zvec);
%this sort of spits out coordinate reprsentations of planes, so that we
%then have coordinates at every xyz point in a cube in increments of 2.

% Now creating an array of initial conditions where the first array is my x
% coordinates etc. So this is actually an indexed array of all the meshgrid
% coords. I guess geometrically it's a tenser, but computationally it's
% just a list of 3d matrix values.
yIC(1, :, :, :) = x0;
yIC(2, :, :, :) = y0;
yIC(3, :, :, :) = z0;

plot3(yIC(1,:), yIC(2,:), yIC(3,:), 'r.' , 'LineWidth', 2, 'Markersize', 10)
axis([-40 40 -40 40 -40 40]) % This establishes axis values
view(20, 40); % This backs the view up
drawnow

% Computer all trajectories
dt = 0.01;
duration  = 4;
tspan = [0:duration]; 
L = duration / dt;  % number of steps in tspan

yparticles = yIC;

for step = 1:L
    time = step*dt;
    for i = 1:length(xvec)
        for j = 1:length(yvec)
            for k = 1:length(zvec)
                y_in = yparticles(:,i,j,k); %this is all coords for ijk particle
                y_out = rk4_singleStep(@(t,y) lorenz(t,y,sigma,beta,rho), dt, time, y_in);
                yparticles(:,i,j,k) = y_out;
            end
        end
    end
    plot3(yparticles(1,:),yparticles(2,:),yparticles(3,:), 'r.', 'LineWidth', 2, 'Markersize', 10)
    axis([-40 40 -40 40 -40 40]);    % This establishes axis values
    view(20, 40); % This backs the view up
    drawnow
end

%% Unfucked to take advantage of matlab's optimized matrix math
clear all; clc; close all;
% Input Lorenz' original parameters that lead to chaos
sigma = 10; beta = 8/3; rho = 28;

% Initial cond is a large cube of points, like a 3d mesh grid?
%xvec = -20:2:20; yvec = -20:2:20; zvec = -20:2:20;
% Now let's try a smaller cube around another point:
y0 = [-8; 8; 27];
xvec = -1:.1:1; yvec = -1:.1:1; zvec = -1:.1:1;


%called it
[x0, y0, z0] = meshgrid(xvec + y0(1), yvec + y0(2), zvec+ y0(3));
%this sort of spits out coordinate reprsentations of planes, so that we
%then have coordinates at every xyz point in a cube in increments of 2.

% Now creating an array of initial conditions where the first array is my x
% coordinates etc. So this is actually an indexed array of all the meshgrid
% coords. I guess geometrically it's a tenser, but computationally it's
% just a list of 3d matrix values.
yIC(1, :, :, :) = x0;
yIC(2, :, :, :) = y0;
yIC(3, :, :, :) = z0;

plot3(yIC(1,:), yIC(2,:), yIC(3,:), 'r.' , 'LineWidth', 2, 'Markersize', 10)
axis([-40 40 -40 40 -40 40]) % This establishes axis values
view(20, 40); % This backs the view up
drawnow
% Computer all trajectories
dt = 0.01;
duration  = 4;
tspan = [0:duration]; 
L = duration / dt;  % number of steps in tspan

y_in = yIC; % we're just gonna pass the entire meshgrid as initial cond.
for step = 1:L
    time = step*dt;
    y_out = rk4_singleStep(@(t,y) lorenz_vector(t,y,sigma,beta,rho), dt, time, y_in);
    y_in = y_out;
    plot3(y_out(1,:), y_out(2,:), y_out(3,:), 'r.' , 'LineWidth', 2, 'Markersize', 10)
    view(20+360*step/L, 40);
    axis([-40 40 -40 40 -40 40]) % This establishes axis values
    drawnow

end





