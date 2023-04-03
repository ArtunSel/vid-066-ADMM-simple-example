function [z_kp1] = minimize_x(x_kp1,y_k,Qx,Qz,A,B,c,rho)
% INPUTS xkp1,y_k,Qx,Qz,A,B,c,rho
% OUTPUTS z_kp1
z=sdpvar(size(Qz,1),1);
x=x_kp1;
y=y_k;
Lp_xz_y=(x'*Qx*x)+(z'*Qz*z)+y'*(A*x+B*z-c)+(rho/2)*norm((A*x+B*z-c),2)^2;
% DIAGNOSTIC = OPTIMIZE(Constraint,Objective,options)
options = sdpsettings('solver','mosek');
DIAGNOSTIC = optimize([],Lp_xz_y,options);
z_kp1=value(z);
end