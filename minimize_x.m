function [x_kp1] = minimize_x(z_k,y_k,Qx,Qz,A,B,c,rho)
% INPUTS z_k,y_k,Qx,Qz,A,B,c,rho
% OUTPUTS x_kp1
x=sdpvar(size(Qx,1),1);
z=z_k;
y=y_k;
Lp_xz_y=(x'*Qx*x)+(z'*Qz*z)+y'*(A*x+B*z-c)+(rho/2)*norm((A*x+B*z-c),2)^2;
% DIAGNOSTIC = OPTIMIZE(Constraint,Objective,options)
options = sdpsettings('solver','mosek');
DIAGNOSTIC = optimize([],Lp_xz_y,options);
x_kp1=value(x);
end