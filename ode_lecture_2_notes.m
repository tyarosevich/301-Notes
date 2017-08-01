%% ODE lecture 2 
% y_k+1 is only appox = y(t_k+1) what we;re trying to find here is y(t_k +
% deltat), or the state one step further forward.
%% Assessing the error of this scheme.

% Let's look at the exact trajectory: y(t_k +deltat) = taylor expans:
% y(t_k) + deltat * dy/dt(t_k) + deltat^2/2! d2y/d2t(c) where c is
% between t_k and t_k+1.
% This then is = y_k + deltat * f(y_k) + O(deltat^2)
% Note that ydot = (y_k), hence the substitution of y_k for dy/dt
% So we see here we have our forward euler plus and error term. So the
% local error, eps_k+1 is = exact trajectory minus our forward euler
% trajectory. or E_k+1 = y(t_k + deltat) - y_k+1 = deltat^2/2 * d2y/d2t(c)
% So every time I take a timestep that timestep if giving an addiional
% deltat^2 error.
% Why this constant and no higher order terms? Because o the mean value
% theorem.
% The global eror is bigE = sum j=0 to k E_k = (b-a)/deltat * deltat^2/2 *
% d2y/d2t(c)
% ==== (b-a)/2 * deltat * d2y/d2t(c)
%% Stability
% suppose ydot = lambda*y, y(0) = y_0
% we know the exact solution is y(t) = e^(lambda * t) * y0
% THis is stable if lambda < 0 because the value gets very small.
% For forward Euler y_k+1 = (1+ lambda*deltat)  (y_k).
% Backward Euler is: y_k+1 = (1 - lambda*deltat)^(-1) (y_k) %Note these
% last y_ks are 'evaluated at' not multiplied.
% How do I tell the stability of a given iteration here?
% I look at the eigenvalues of A, or in this case cause there is just the
% one number and no matrix, therefore the eigenvalues ARE the values
% themselves. THis will be stabe i these eigs are inside the unit circle.
% So in these cases the eigs are 1+lambda*deltat and 1-lambea*deltat.

% So in FE y_0 goes to y_1 etc and y_n = (1+lambda*deltat)^N (y_0)
% This expression is unstable when the modulus or length of abs(1+deltat) >
% 1. 
% Similarly for Backward Euler unstable when the magnitude of
% 1/(1-lambda*deltat) > 1.
% So in the actual system lambda is all that matters, and it just needs to
% be small? negative? But in our numerical scheme the stability involves
% the deltat value. 
% So for z = lambda*deltat there is a complex plane with a unit circle
% centered -1, 0  and only lambda values in this unit circle are stable.
% Consider the same for backward euler. Now it's a unit circle around 1,0
% and this is the only region where anything is UNstable.
%% Suppose no this is a matrix system:

% ydot = Ay     Now it's just that this system is stable i all eigs of A
% are inside the unit circle. Eveything is pretty similar:
% FE is now (I + deltat * A) y_k
% BE is (I-deltat*A)^-1 y_k
% So FE we can say y_N = (I + deltat*A)^N * y_0
   %  This is stable when all eigs of I+deltat*A have mag less than 1.
% FOr BE : y_N = (I - deltat*A)^-N * y_0
%   This is stable if and only if (iff) eig(I - deltat*A)^-1 have magnitude
%   < 1 in the complex plane.
%% Take the example o the single pendulum. M on a L that can pivot with angle th and M is being accelerated by G. 
% KE T = 1/2 ML th'^2
% V = LMG (1-cos(th))====> L = T - V  
% d /d dL/dth - dlL/dth = 0 ===> xdotdot = G/L sin(x)
% make xdot = v and v = G/L (sin(x)):
% d/dt[x;v] = [f_1(x,c); f_2(x,v)] 
%% Coding this up

clear all
% dy = pend(t,y,g,L,d)
t = 0:.1:50;
y0 = [pi/4; 0];
g = -10;
L = 10;
d = .1;
[t,y] = ode45(@(t,y) pend(t,y,g,L,d), t, y0);

figure
plot(t,y(:,1));
figure
plot(y(:,1), y(:,2));
figure
plot3(t,y(:,1), y(:,2));

