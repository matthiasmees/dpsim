# Two generators and operating-point initialization

`FMI_EMT_Ph3_TwoGenerators.cpp` connects two separately exported nonlinear FMUs to one shared EMT ABC bus. `sources/` contains the second model/FMU/XML and references the first retained generator export and configuration. `matlab/` contains the independent full-state trapezoidal and ode15s references and saved results. `Comparison.ipynb` is the executed comparison, including independent controls, singleton rejection, model-order isolation, tightened tolerances and timestep refinement. Generated outputs remain under build.
