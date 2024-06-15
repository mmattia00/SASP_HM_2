function [b] = antiparallel_diodes(a, Z_diode)
% ANTIPARALLELDIODES This MATLAB function implements the nonlinear
% scattering equation of a pair of diodes in antiparallel configuration
%   a : wave incident to the nonlinear element
%   Z_diode : reference port resistance (set to adapt the junction to which
%             the nonlinear element is connected)

Vt = 25.864e-3;   % Thermal Voltage (assuming 27°C in operation condition)
Is = 2.52e-9;     % Saturation current of 1N914
eta = 1.752;      % Ideality Factor 1N914

W_x = @(a, Z_diode) Lambert_W_Fritsch((Z_diode*Is/(eta*Vt))*exp((Z_diode*Is+abs(a))/(eta*Vt)));

b = sign(a)*(abs(a)+2*Z_diode*Is-2*(eta*Vt)*W_x(a, Z_diode));

end

