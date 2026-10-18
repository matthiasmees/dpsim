function two_generator_reference
source=fileparts(mfilename('fullpath'));root=fileparts(fileparts(fileparts(fileparts(source))));
addpath(fullfile(root,'dpsim-fmi/examples/synchronous_generator/matlab/nonlinear_3ph'));
b=sg_benchmark_parameters;folder=fullfile(root,'build-macos-check/fmi-validation/two_generators');if ~isfolder(folder),mkdir(folder);end
if isfile(fullfile(folder,'reference_ode15s.csv')),fprintf('TWO_REFERENCE_EXISTS\n');return;end
x0=[b.x0(1:10);b.x0(1:10);b.x0(11:13);2*b.x0(14:16)];
d0=rhs(0,x0,b,[1 1]);
assert(max(abs(d0([1:4 6:10 11:14 16:20])))<1e-8,'Not an operating point');
initial=struct('equilibrium_residual',max(abs(d0([1:4 6:10 11:14 16:20]))),'initial_voltage_V',x0(21:23),'initial_line_current_A',x0(24:26),'generator_power_each_W',b.initialPgen);
fid=fopen(fullfile(folder,'operating_point.json'),'w');fwrite(fid,jsonencode(initial,PrettyPrint=true));fclose(fid);
for h=[50e-6 25e-6]
 t=(0:round(.1/h))'*h;x=zeros(numel(t),26);x(1,:)=x0';maxResidual=0;
 for k=2:numel(t)
  refs=[1+.02*((t(k)+t(k-1))/2>=.05),1-.01*((t(k)+t(k-1))/2>=.075)];
  prev=x(k-1,:)';f0=rhs(t(k-1),prev,b,refs);trial=prev+h*f0;
  for iteration=1:12
   f=rhs(t(k),trial,b,refs);res=trial-prev-h/2*(f0+f);
   scaled=max(abs(res)./(1e-11+1e-11*max(abs(trial),abs(prev))));
   if scaled<1,break;end
   jac=zeros(26);for j=1:26,z=trial;z(j)=z(j)+1i*1e-20;jac(:,j)=imag(rhs(t(k),z,b,refs))/1e-20;end
   trial=trial-(eye(26)-h/2*jac)\res;
  end
  assert(scaled<1,'Reference Newton did not converge');maxResidual=max(maxResidual,scaled);x(k,:)=trial';
 end
 saveResult(fullfile(folder,sprintf('reference_trap_%gus.csv',h*1e6)),t,x,b);
 fprintf('TWO_REFERENCE_DONE h=%g residual=%g\n',h,maxResidual);
end
% Independently integrated, segmented high-accuracy ode15s cross-check.
t=(0:2000)'*50e-6;x=zeros(numel(t),26);state=x0;bounds=[0 .05 .075 .1];
opts=odeset('RelTol',2e-10,'AbsTol',1e-11,'MaxStep',10e-6);
for k=1:3
 a=bounds(k);z=bounds(k+1);refs=[1+.02*(a>=.05),1-.01*(a>=.075)];
 sol=ode15s(@(time,y)rhs(time,y,b,refs),[a z],state,opts);mask=t>=a & t<=z;x(mask,:)=deval(sol,t(mask))';state=sol.y(:,end);
end
saveResult(fullfile(folder,'reference_ode15s.csv'),t,x,b);
end
function dx=rhs(t,x,b,refs)
d1=sg_full_rhs_reference(t,[x(1:10);x(21:26)],b,2*b.G0,refs(1));
d2=sg_full_rhs_reference(t,[x(11:20);x(21:26)],b,2*b.G0,refs(2));
[~,i1]=sg_machine_rhs(x(1:6),x(21:23),x(9),x(10),b.machine);
[~,i2]=sg_machine_rhs(x(11:16),x(21:23),x(19),x(20),b.machine);
v=x(21:23);il=x(24:26);grid=sqrt(2)*b.Vphase*cos(b.wb*t+b.phase);
dx=[d1(1:10);d2(1:10);(i1+i2-il-2*b.G0*v)/(2*b.network(1));(v-grid-b.network(2)/2*il)/(b.network(3)/2)];
end
function saveResult(path,t,x,b)
columns={'time'};data=t;
names={'i_a','i_b','i_c','speed_pu','theta_e_rad','field_current_A','p_gap_W','p_gen_W','v_measured_pu','avr_raw_pu','avr_limited_pu','pm_pu','governor_raw_pu','governor_limited_pu','v2_filter_pu','avr_integral_pu','efd_pu','psi_d_Wb','psi_q_Wb','psi_0_Wb','psi_f_Wb','theta_state_rad','speed_state_pu'};
for n=1:2
 offset=(n-1)*10;y=zeros(numel(t),numel(names));
 for k=1:numel(t)
  state=x(k,offset+(1:10))';voltage=x(k,21:23)';
  [~,current,machine]=sg_machine_rhs(state(1:6),voltage,state(9),state(10),b.machine);
  reference=1;if n==1 && t(k)>=.05-1e-12,reference=1.02;elseif n==2 && t(k)>=.075-1e-12,reference=.99;end
  [~,~,avr]=sg_avr_rhs(state(7:9),voltage,reference,b.avr);
  [~,governor]=sg_governor_rhs(state(10),state(6),b.pm0,b.governor);
  y(k,:)=[current;machine;avr;state(10);governor;state(7:9);state(1:6)]';
 end
 prefixes={'first','second'};prefix=prefixes{n};
 for j=1:numel(names),columns{end+1}=[prefix '_' names{j}];end
 data=[data,y];
end
for phase=1:3,columns{end+1}=['v_pcc_' char('a'+phase-1)];end
for phase=1:3,columns{end+1}=['i_line_' char('a'+phase-1)];end
data=[data,x(:,21:26)];writetable(array2table(data,'VariableNames',columns),path);
end
