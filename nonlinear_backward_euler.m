function [t, x] = nonlinear_backward_euler(dx_dt, df_dx, x0, t_end, dt, tol)

t = 0:dt:t_end;
n_rows = length(x0);
x = zeros(n_rows, length(t));
x(:,1) = x0;

for i = 1:length(t-1)%this is the backward euler. length -1 cause we already have x0
    G = @(new_x) new_x - dt * dx_dt ( t(i+1), new_x ) - x(:,i);
    DG = @(new_x) eye(n_rows) - dt*df_dx(t(i+1), new_x);
    
    next_x = newton_2(G, DG, x(:,1), tol);
    x(i+1) = next_x(:,end);
end