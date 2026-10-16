function [dx, diagnostic] = sg_network_rhs(x, i_gen, v_grid, g_load, p)
%#codegen
% x = [v_bus_abc; i_line_abc], line current positive from bus to grid.
% p = [Cbus Rline Lline], load conductance in S PER PHASE.
v=x(1:3); iline=x(4:6);
i_load=g_load*v;
dx=[(i_gen-iline-i_load)/p(1); (v-v_grid-p(2)*iline)/p(3)];
diagnostic=[v.'*i_load; v.'*iline; p(2)*(iline.'*iline)];
end
