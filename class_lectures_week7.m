%% Wednesday Class Lecture 8/2
% Non linear ODEs.
%Example: dx/dt = f(t,x) = -x^2 , x0 = 1 , Exact solution: x=1/t+1
% With backward Euler this becomes problematic because it's implicit, and
% may not even be possible to write an easily solved equation equal to
% x_k+1
% To solve for x_k+1 in this situation we need to set up a root-finding
% problem.
% i.e. the backward euler: xk+1 - dt(-xk+1 ^ 2) = xk
% Root solving problem is:
% g(x) = x - dt * f(tk+1,  x) -xk   
% What we're doing is moving everythng over setting the BE to zero, then
% solving that equation AT zero to get the particular x_k+1 that satisfies
% the equation. 
% we can then newton iterate this:
% DG(x) = I - dt*df/dx * (tk+1, x)
% We'll need to iterate :
% xold - xnew = DG(xold) \ G(xold)
%% 
colors = {'r', 'b', 'g'};

dx_dt = @(t,x) -x^2;
x0 = 1;
t_end = 1;
dt = [.1, .01, .001];
% uses some function.
% t_exact = linspace(t, t_end, 100);
% plot(t_exact, 1./(t_exact + 1), 'k', 'LineWidth', .2)
% for i = 1:length(dt)
%     [t,x] = forward_euler2(dx_dt, t_end, dt(i));
% end
% plot(t,x,colors(i), 'LineWidth', 2)
dx_dt = @(t,x) -x^2;
df_dx = @(t,x) -2*x;
x0=1;
t_end = 1;
dt = .1;
tol = 1e-6;

[t,x] = nonlinear_backward_euler(dx_dt, df_dx, x0, t_end, dt, tol);
t_exact = linspace(0,1,100);
plot(t_exact, 1./(t_exact + 1), 'k' , 'LineWidth', 2);
hold on 
plot(t,x, 'b', 'LineWidth', 2);
