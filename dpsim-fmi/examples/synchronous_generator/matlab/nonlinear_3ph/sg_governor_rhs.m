function [dx, diagnostic] = sg_governor_rhs(x, speed_pu, pmref_pu, p)
%#codegen
% One first-order prime mover with droop, not a turbine-specific model.
% p = [droop Tgovernor Pmin Pmax].
raw = pmref_pu+(1-speed_pu)/p(1);
limited = min(max(raw,p(3)),p(4));
dx = (limited-x)/p(2);
diagnostic=[raw; limited];
end
