function [hx,P,v] = kalman_filter(A,H,Q,R,x0,P0,y)
% function [hx, P, v] = kalman_filter(A,H,Q,R,x0,P0,y)
%
% One Prediction and Update step of Classic Kalman Filter for the nominal 
% system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  A,H,Q,R -> matrices of the system
%         x0, P0 -> initial conditions      
%        
% output: hx -> final estimate
%         P -> final Covariance
%         v -> estimation variation
%
%
% Last Update: 01/04/2018
% Author: Marcos Rogério Fernandes
% E-mail: eng.marofe@hotmail.com
% Site: https://marofe.github.com


 %% Prediction Step
    hx=A*x0;
    P=A*P0*A'+Q;  
 %% Update Step
    P=P-P*H'*inv(R+H*P*H')*H*P;
    %P=inv(inv(P)+H'*inv(R)*H);
    K=P*H'*inv(R);
    v=K*(y-H*hx);
    hx=hx+v;
end
