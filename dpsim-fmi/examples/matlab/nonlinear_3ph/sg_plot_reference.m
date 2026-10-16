function sg_plot_reference(b,t,y,columns,outputDir)
fig=figure('Name','Nonlinear synchronous-generator grid','Position',[50 50 1250 1100]);
tiledlayout(4,2,'TileSpacing','compact');
nexttile; plot(t,y(:,24:26)); xlim([.97,1.09]); grid on;
title('Load step: bus EMT voltages'); ylabel('V'); legend('a','b','c');
nexttile; plot(t,y(:,1:3)); xlim([.97,1.09]); grid on;
title('Load step: generator currents'); ylabel('A'); legend('a','b','c');
nexttile; plot(t,y(:,9),t,y(:,31),'--'); grid on;
ylabel('Voltage [pu]'); legend('Measured','Reference');
nexttile; plot(t,b.frequency*(y(:,4)-1)); grid on;
ylabel('Rotor speed deviation [Hz]');
nexttile; plot(t,y(:,8)/1000,t,y(:,32)/1000,t,y(:,33)/1000, ...
    t,y(:,12)*b.Sbase/1000); grid on;
ylabel('Instantaneous power [kW]'); xlabel('Time [s]');
legend('Generator','Load','Line at bus','Mechanical','Location','best');
nexttile; plot(t,y(:,17),t,y(:,11),'--',t,y(:,16),':'); grid on;
ylabel('Excitation [pu]'); xlabel('Time [s]'); legend('Field supply','Limited command','PI state');
nexttile; plot(t,(y(:,5)-b.wb*t+pi/2)*180/pi); grid on;
ylabel('Rotor angle relative to grid [deg]'); xlabel('Time [s]');
nexttile; plot(t,y(:,12),t,y(:,13),'--',t,y(:,14),':'); grid on;
ylabel('Prime mover [pu]'); xlabel('Time [s]');
legend('Actual power','Droop command','Limited command');
sgtitle('Nonlinear EMT grid: load +6 kW at 1 s, AVR reference +2% at 2 s');
savefig(fig,fullfile(outputDir,'sg_3ph_reference.fig'));
exportgraphics(fig,fullfile(outputDir,'sg_3ph_reference.png'),'Resolution',160);
% Preserve readable signal metadata separately from the numeric CSV.
f=fopen(fullfile(outputDir,'signal_names.json'),'w');
assert(f~=-1,'Cannot write signal metadata.');
c=onCleanup(@() fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(columns));
end
