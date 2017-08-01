%% Given a function of data f(t)
% df/dt =  lim t--> 0 (f(t+deltat) - f(t)) / deltat
% we are going to taylor expand f(t + deltat). for this equation see
% lecture notes for Monday 7/24

% f(x) = sin(x). We are going to taylor expand around x = 0. If you taylor
% expand around x = 0 it's call the Maclaurin Series. 
% x is the new point, a is the known point:
% f(x) = sin(0) + cos(x) * (x) + 0 - 1* (x - a)^2/2
% more simply: x - x^3/3! + x^5/5! - x^7/7!...
% This generalizes to vector functions when x = [x1; x2; etc] and f = [f1,
% f2, etc] then df/dx^ = [df/dx1 df/dx2 ... df/dxn
%                         df2/dx  df2/dx2 etc.]
%% Actual matlab example
clear all
x = -10:.01:10;
y = sin(x);
plot(x,y,'k', 'Linewidth', 2)
axis([-10 10 -10 10])
grid on, hold on
%% first order taylor xpac
P = [1 0]; %polynomial for x + 0
yt1 = polyval(P, x);
plot(x, yt1, 'b--', 'Linewidth', 1.2)
%% third order taylor (remember our taylor xpac only has odd terms)
P = [-1/factorial(3) 0 1 0]; % this is the polynomial values
yt3 = polyval(P, x);
plot(x, yt3, 'r--', 'Linewidth', 1.2)
%% fifth order
P = [1/factorial(5) 0 -1/factorial(3) 0 1 0];
yt5 = polyval(P, x);
plot(x, yt5, 'g--', 'Linewidth', 1.2)
% This can keep goin and oging. A good rule of thumb is that you need 2
% polynomial expansions for each half period? Iono.
%% Going back to numerical approximation:
% Forward difference: just the same equation with (f(t + deltat) -
% f(t))/deltat = df/dt (which is what we want) + the error. The error is
% the rest of the taylor expansion: deltat/2 * d2f/d2t(t) + etc..
% error is on the order of delta t, or O(deltat)

%backward diff:
% (f(t) - f(t - deltat)) / deltat 
% Talor xpac for this:
% f(t - deltat) = f(t) - deltat * df/dt(t) + deltat^2/2! * d2df/d2dt(t) -
% deltat^3/3! * d3f/d3t(t) ... 

% We then plug this into the expresion for backward diff, since f(t-deltat)
% matches up with part of the backward diff equation.
% I.e. (f(t) - (taylor xpansion here)) / deltat
%% Finally we have central difference

% very logical : (f(t + deltat) - f(t - deltat)) / 2* deltat
% = 2* deltat * df/dt(t) + 2* deltat^2/3! * d3f/d3t(t) + 0(0t^5))
% divide this by 2* delta t = df/dt(t) + deltat^2/3! *d3f/d3t (t) .. 
% As we can see we the error is the order of deltat^2, which in this case
% is a good thing, cause as we cut our deltat our error decreases at a
% square rate.
%% Coding up an example

clear all;close all;

dt = .3;
t = -2:dt:4;
f = sin(t);
dfdt = cos(t); % exact deriv for comparison
plot(t, f, 'k--', 'Linewidth', 1.2)
hold on, grid on
plot(t, dfdt,'k', 'Linewidth', 3)
l1 = legend('Function', 'Exact Derivative');
set(l1, 'FontSize', 14)
axis([-2 4 -1.5 1.5])

%% Forward diff
dfdtF = (sin(t + dt) - sin(t)) / dt;
%backward
dfdtB = (sin(t) - sin(t - dt))/dt;
%central
dfdtC = (sin(t + dt) - sin(t - dt))/(2*dt)  ;
plot(t, dfdtF, 'b', 'Linewidth', 1.2) % forward
plot(t, dfdtB, 'g', 'Linewidth', 1.2) % back
plot(t, dfdtC, 'r', 'Linewidth', 1.2) % central
l1 = legend('Function', 'Exact Derivative', 'Forward', 'Backward', 'Central');
set(l1, 'FontSize', 14)






%% Taylor series aside:
% In general it uses a given point and diff information at that point to
% get an estimate about nearby points