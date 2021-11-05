function [dRHO,dOMEGA] = next_step(RHO,THETA,OMEGA,RHO_star,THETA_star,OMEGA_star,Bg,Bh,param_rho,param_omega)

Kp=0.3;
Kw=0.5;

M1_11=-Kp;
M1_12=Gauss_Regression(Bg,RHO,THETA,THETA_star,param_rho,1);
M1_21=-Kp*Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO_star,THETA_star],param_omega,1);
M1_22=Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO_star,THETA_star],param_omega,2)+Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO_star,THETA_star],param_omega,1)*Gauss_Regression(Bg,RHO,THETA,THETA_star,param_rho,1)-Kw;
M1=[M1_11,M1_12;M1_21,M1_22];

M2_11=Kp;
M2_12=0;
M2_21=Kp*Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO_star,THETA_star],param_omega,1);
M2_22=Kw;
M2=[M2_11,M2_12;M2_21,M2_22];

x = [RHO_star;OMEGA_star];
xd = [Gauss_Regression(Bg,RHO,THETA,THETA_star,param_rho);Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO_star,THETA_star],param_omega)];
a=M1*x+M2*xd;
dRHO=a(1);
dOMEGA=a(2);