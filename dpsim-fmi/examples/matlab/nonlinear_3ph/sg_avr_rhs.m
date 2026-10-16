function [dx, efd_pu, diagnostic] = sg_avr_rhs(x, v_abc, vref_pu, p)
%#codegen
% x = [filtered squared RMS voltage pu; PI integrator; exciter voltage pu].
% p = [VphaseBase Tmeasure Kp Ki Kaw EfdMin EfdMax Texciter].
v2 = (v_abc.'*v_abc)/(3*p(1)^2);
vmeasure = sqrt(max(x(1),1e-12));
error = vref_pu-vmeasure;
raw = x(2)+p(3)*error;
limited = min(max(raw,p(6)),p(7));
dx = [(v2-x(1))/p(2); p(4)*error+p(5)*(limited-raw); ...
      (limited-x(3))/p(8)];
efd_pu=x(3);
diagnostic=[vmeasure; raw; limited];
end
