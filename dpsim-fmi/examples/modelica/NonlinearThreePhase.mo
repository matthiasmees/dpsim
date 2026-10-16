within ;
package NonlinearThreePhase
  // Port adaptation of MSL 4.0.0 Electrical.Analog.Examples.ShowSaturatingInductor.
  // Same parallel inductors, parameters, initial currents, and source waveform.
  model PhasePort
    input Real v;
    output Real i;
    Modelica.Electrical.Analog.Sources.SignalVoltage source;
    Modelica.Electrical.Analog.Basic.Ground ground;
    Modelica.Electrical.Analog.Basic.SaturatingInductor saturated(
      Lzer=2, Lnom=1, Inom=1, Linf=0.5, i(start=0, fixed=true));
    Modelica.Electrical.Analog.Basic.Inductor linear(L=1, i(start=0, fixed=true));
  equation
    source.v=v;
    connect(source.n, ground.p);
    connect(source.p, saturated.p);
    connect(source.n, saturated.n);
    connect(source.p, linear.p);
    connect(source.n, linear.n);
    i=saturated.i+linear.i;
  end PhasePort;
  model Ported
    input Real v_a; input Real v_b; input Real v_c;
    output Real i_a; output Real i_b; output Real i_c;
    PhasePort a, b, c;
  equation
    a.v=v_a; b.v=v_b; c.v=v_c;
    i_a=a.i; i_b=b.i; i_c=c.i;
  end Ported;
  model Reference
    Modelica.Electrical.Analog.Examples.ShowSaturatingInductor a;
    Modelica.Electrical.Analog.Examples.ShowSaturatingInductor b(phase=Modelica.Constants.pi/2-2*Modelica.Constants.pi/3);
    Modelica.Electrical.Analog.Examples.ShowSaturatingInductor c(phase=Modelica.Constants.pi/2+2*Modelica.Constants.pi/3);
    output Real i_a=a.SaturatingInductance1.i+a.Inductance1.i;
    output Real i_b=b.SaturatingInductance1.i+b.Inductance1.i;
    output Real i_c=c.SaturatingInductance1.i+c.Inductance1.i;
    annotation(experiment(StopTime=6.2832, Interval=0.0001, Tolerance=1e-10));
  end Reference;
end NonlinearThreePhase;
