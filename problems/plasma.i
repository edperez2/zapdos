[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 10
  xmax = 5e-3 # Length of test chamber
  ymax = 5e-4 # Test chamber radius
  uniform_refine = 3
[]

[Variables]
  active = 'electron_density ion_density potential'
  [./electron_density]
    #Adds a Linear Lagrange variable by default
    #scaling = 1e0
  [../]
  [./ion_density]
  [../]
  [./potential]
    #scaling = 1e0
    initial_condition = 0
  [../]
[]

[AuxVariables]
  [./e_field_mag]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./e_field_mag_kernel]
    type = EFieldMag
    variable = e_field_mag
    potential = potential
    execute_on = timestep_end
    #component = x
  [../]
[]

[Functions]
  [./parsed_function]
    type = ParsedFunction
    #value = '1e13 + 1e19*exp(-pow(x-2.5e-3,2)/pow(100.e-6,2))*exp(-pow(y-2.5e-4,2)/pow(10.e-6,2))'
    value = '1.0*exp(-pow(x-2.5e-3,2)/pow(100.e-6,2))'
  [../]
[]

[ICs]
  active = 'electron_density_ic ion_density_ic'
  [./electron_density_ic]
    type = FunctionIC
    variable = 'electron_density'
    function = parsed_function
  [../]
  [./ion_density_ic]
    type = FunctionIC
    variable = 'ion_density'
    function = parsed_function
  [../]
[]

[Kernels]
  active = 'electrons_time_deriv electron_diffusion poisson_diffusion electron_convection electron_src ions_src ions_time_deriv poisson_src' 
  [./poisson_diffusion]
    type = Diffusion
    variable = potential
  [../]
  [./poisson_src]
    type = PoissonSource
    variable = potential
    #permittivity = 8.85e-12
    ion_density = ion_density
    electron_density = electron_density
  [../]
  [./electrons_time_deriv]
    type = TimeDerivative
    variable = electron_density
  [../]
  [./electron_convection]
    type = BovineConvection
    variable = electron_density
    some_variable = potential
  [../]
  [./electron_diffusion]
    type = ConstTimesDiffusion
    variable = electron_density
    diffusion_coeff = 0.18
  [../]
  [./electron_src]
    type = IonizationSource
    variable = electron_density
    potential = potential
  [../]
  [./ions_time_deriv]
    type = TimeDerivative
    variable = ion_density
  [../]
  [./ions_src]
    type = CoupledIonizationSource
    variable = ion_density
    potential = potential
    electron_density = electron_density
    #ionization_coeff = 0.35e-6
  [../]
[]

[BCs] 
  [./anode]
    type = DirichletBC # Simple u=value BC
    variable = potential
    boundary = left # Name of a sideset in the mesh
    value = 20.0e-1 #
  [../]
  [./cathode]
    type = DirichletBC
    variable = potential
    boundary = right
    value = 0
  [../]
[]

[Problem]   
  type = FEProblem # This is the "normal" type of Finite Element Problem in MOOSE   
  # coord_type = XYZ # Cartesian   
[]   

[Preconditioning]
  [./SMP_jfnk_full]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'lu'
  [../]
[]

[Executioner]
  type = Transient
  dt = 1.0e-7
  #solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type'
  #petsc_options_value = 'hypre boomeramg'
  end_time = 1.0e-4
#  trans_ss_check = true
#  ss_check_tol = 1e-8
  nl_rel_tol = 1e-2
  l_tol = 1e-1
#  nl_abs_tol = 1e-3
  l_max_its = 500
#  [./TimeStepper]
#    type = IterationAdaptiveDT
#    linear_iteration_ratio = 10
#    cutback_factor = 0.67
#    dt = 1.0e-7
#    growth_factor = 1.5
#    optimal_iterations = 4
#  [../]
[]

[Adaptivity]
  marker = error_frac
  max_h_level = 3
  [./Indicators]
    [./temp_jump]
      type = GradientJumpIndicator
      variable = potential
      scale_by_flux_faces = true
    [../]
  [../]
  [./Markers]
    [./error_frac]
      type = ErrorFractionMarker
      coarsen = 0.1
      indicator = temp_jump
      refine = 0.6
    [../]
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  output_initial = true # Output initial condiiton
  exodus = true # Output Exodus format
  print_perf_log = true # Show performance log information on screen
  print_linear_residuals = true # Display linear residuals on screen
[]