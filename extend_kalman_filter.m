function [hx,P,v] = extend_kalman_filter(f,h,A,H,Q,R,x0,P0,y,u)
% function [hx, P, v] = extend_kalman_filter(f,h,A,H,Q,R,x0,P0,y,u)
%
% One Prediction and Update step of Classic Kalman Filter for the nominal
%system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  A,H,Q,R -> matrices of the system
%         x0, P0 -> initial conditions
%
% output: hx -> final estimate
%         P -> final Covariance

%
% Author: Marcos Rogério Fernandes
% E-mail: eng.marofe@hotmail.com
% Date: 15/01/2019

%% prediction
hx=f(x0,u);
P=A*P0*A'+Q;
%% update
P=inv(inv(P)+H'*inv(R)*H);
K=P*H'*inv(R);
v=K*(y-h(hx));
hx=hx+v;
end