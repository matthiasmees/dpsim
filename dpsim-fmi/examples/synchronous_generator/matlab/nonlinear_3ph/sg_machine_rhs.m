function [dx, i_gen, diagnostic] = sg_machine_rhs(x, v_abc, efd_pu, pm_pu, p)
%#codegen
% Power-invariant Park coordinates, passive stator current INTO machine.
% x = [psi_d; psi_q; psi_0; psi_f; theta_e; speed_pu].
% p = [Rs Ld Lq L0 Lf M Rf wb Sbase H D Vfbase].
Rs=p(1); Ld=p(2); Lq=p(3); L0=p(4); Lf=p(5); M=p(6);
Rf=p(7); wb=p(8); Sbase=p(9); H=p(10); D=p(11); Vfbase=p(12);
a = x(5) + [0; -2*pi/3; 2*pi/3];
T = sqrt(2/3)*[cos(a).'; -sin(a).'; ones(1,3)/sqrt(2)];
v = T*v_abc;
% Symmetric flux/current matrix: magnetic energy is positive definite.
id = (Lf*x(1)-M*x(4))/(Ld*Lf-M*M);
ifield = (Ld*x(4)-M*x(1))/(Ld*Lf-M*M);
iq=x(2)/Lq; i0=x(3)/L0;
i_passive = [id; iq; i0];
i_gen = -T.'*i_passive; % positive OUT of generator, INTO bus
we = wb*x(6);
% For a power-invariant transform, no 3/2 factor appears in power/torque.
p_gap = -we*(x(1)*iq-x(2)*id);
dx = zeros(6,1);
dx(1) = v(1)-Rs*id+we*x(2);
dx(2) = v(2)-Rs*iq-we*x(1);
dx(3) = v(3)-Rs*i0;
dx(4) = Vfbase*efd_pu-Rf*ifield;
dx(5) = we;
% Kinetic energy = H*Sbase*speed_pu^2; D is damping in power per pu speed.
dx(6) = (pm_pu-p_gap/Sbase-D*(x(6)-1))/(2*H*x(6));
diagnostic = [x(6); x(5); ifield; p_gap; v_abc.'*i_gen];
end
