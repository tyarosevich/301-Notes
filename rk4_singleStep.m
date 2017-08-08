function yout = rk4_singleStep(f, dt, tk, yk) 
% f is my right hand side of ydot etc

f1 = f(tk, yk);
f2 = f(tk+dt/2, yk + (dt/2) * f1);
f3 = f(tk + dt/2, yk + (dt/2) * f2);
f4 = f(tk + dt, yk + dt*f3);
yout = yk + dt/6 * (f1 + 2*f2 + 2*f3 + f4);
