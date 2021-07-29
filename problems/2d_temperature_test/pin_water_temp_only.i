# Sets up a pin-to-water geometry, with gas flowing directly down from an 
# inlet tube at the axis of symmetry.
# A pin electrode is at the center of the tube with a tip temperature of 
# 1440 K. 
# The electrode tip is 1 mm away from the water surface.
# 
# Gas flows through the tube and is incident on the water surface.
# Gas flow will form a small convective cell in the liquid, which carries 
# heat away from the water surface. 

dom0Scale=1.0
dom1Scale=1.0
dom2Scale=1.0

[GlobalParams]
  offset = 20
  #supg = true
  #pspg = true
  potential_units = kV
  use_moles = true
  convective_term = true
  transient_term = true
  integrate_p_by_parts = true
  #supg = true
  #integrate_p_by_parts = false
  gravity = '0 0 0'
[]

[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = 'pin_water_full_out_04.e' 
    use_for_exodus_restart = true
  []
  #[geo]
  #  type = FileMeshGenerator
  #  file = 'pin_water_full_mesh.msh'
  #[]
  #
  #[interface1]
  #  type = SideSetsBetweenSubdomainsGenerator
  #  primary_block = '0'
  #  paired_block = '1'
  #  new_boundary = 'gas_bottom'
  #  input = geo
  #[]
  #[interface2]
  #  type = SideSetsBetweenSubdomainsGenerator
  #  primary_block = '1'
  #  paired_block = '0'
  #  new_boundary = 'water_top'
  #  input = interface1
  #[]
  #
  ## Creating a corner node to pin pressure (p = 0)
  #[corner_node]
  #  type = ExtraNodesetGenerator
  #  new_boundary = top_right_corner
  #  coord = '0.03334 0.01872729'
  #  input = interface2
  #[]
[]

[Problem]
  type = FEProblem
  coord_type = RZ
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 1e6
  automatic_scaling = true
  compute_scaling_once = false
  line_search = 'basic'
  petsc_options = '-snes_converged_reason'
  solve_type = newton
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -snes_stol'
  petsc_options_value = 'lu NONZERO 1.e-10 0'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  dtmin = 1e-17
  l_max_its = 100
  nl_max_its = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-6
  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-4
    growth_factor = 1.4
    optimal_iterations = 10
  []
[]

[Outputs]
  #[full_out_01]
  #[out_20mA_01]
  [out_25mA_01]
  #[out_30mA_01]
  #[out_35mA_01]
  #[out_40mA_01]
  #[out_45mA_01]
  #[out_50mA_01]
    type = Exodus
    interval = 1
    output_material_properties = true
    show_material_properties = 'thermal_conductivity emissivity'
  []
[]

[Variables]
  #[vel_x]
  #  # Velocity in radial (r) direction
  #  family = LAGRANGE
  #  order = SECOND
  #  block = '0 1'
  #  initial_condition = 0
  #[]
  #[vel_y]
  #  # Velocity in axial (z) direction
  #  family = LAGRANGE
  #  order = SECOND
  #  block = '0 1'
  #  initial_condition = 0
  #[]
  #[p]
  #  family = LAGRANGE
  #  order = FIRST
  #  block = '0 1'
  #  initial_condition = 0
  #[]
  [T]
    family = LAGRANGE
    order = SECOND
    #order = FIRST
    initial_condition = 300
  []
[]

[Kernels]
  ## Navier-Stokes problem test
  #[x_momentum_time]
  #  type = INSMomentumTimeDerivative
  #  variable = vel_x
  #  block = '0 1'
  #[]
  #[y_momentum_time]
  #  type = INSMomentumTimeDerivative
  #  variable = vel_y
  #  block = '0 1'
  #[]
  #[mass]
  #  type = INSMassRZ
  #  variable = p
  #  u = vel_x
  #  v = vel_y
  #  p = p
  #  block = '0 1'
  #[]
  #[x_momentum_space]
  #  type = INSMomentumLaplaceFormRZ
  #  variable = vel_x
  #  u = vel_x
  #  v = vel_y
  #  p = p
  #  component = 0
  #  block = '0 1'
  #[]
  #[y_momentum_space]
  #  type = INSMomentumLaplaceFormRZ
  #  variable = vel_y
  #  u = vel_x
  #  v = vel_y
  #  p = p
  #  component = 1
  #  block = '0 1'
  #[]

  [temperature_time]
    type = INSTemperatureTimeDerivative
    variable = T
    block = '0 1'
  []
  [temperature_space]
    type = INSTemperature
    variable = T
    u = vel_x
    v = vel_y
    block = '0 1'
  []

  
  # Heat conduction through solid regions
  [HeatTdot_dielectric]
    type = ADHeatConductionTimeDerivative
    variable = T
    block = '2 3'
  []
  [HeatDiff_dielectric]
    type = ADHeatConduction
    variable = T
    block = '2 3'
  []
[]

[AuxVariables]
  [vel_x]
    initial_from_file_var = aux_vel_x
    initial_from_file_timestep = 'LATEST'
    block = '0 1'
  []
  [vel_y]
    initial_from_file_var = aux_vel_y
    initial_from_file_timestep = 'LATEST'
    block = '0 1'
  []
  [p]
    initial_from_file_var = aux_p 
    initial_from_file_timestep = 'LATEST'
    block = '0 1'
  []

  # Heat flux
  #[heat_x]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[]
  #[heat_y]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  block = 0
  #[]
[]

#[AuxKernels]
#  [heat_x_calc]
#    type = HeatFlux
#    variable = heat_x
#    u = vel_x 
#    v = vel_y
#    temperature = T
#    component = 0
#    position_units = ${dom0Scale}
#    block = 0
#  []
#  [heat_y_calc]
#    type = HeatFlux
#    variable = heat_y
#    u = vel_x 
#    v = vel_y
#    temperature = T
#    component = 1
#    position_units = ${dom0Scale}
#    block = 0
#  []
#[]

  # Total heat flux to the water surface is just a surface integral
#[Postprocessors]
#  [flux_down]
#    type = SideDiffusiveFluxIntegral
#    variable = T
#    diffusivity = k
#    boundary = 'gas_bottom'
#  []
#
#  # h calculation
#  [h_tip]
#    type = HeatTransferCoefficientAverage
#    variable = T
#    diffusivity = k
#    bulk_temperature = 300
#    boundary = 'electrode_tip electrode_needle'
#  []
#  [h_side_short]
#    type = HeatTransferCoefficientAverage
#    variable = T
#    diffusivity = k
#    bulk_temperature = 300
#    boundary = 'electrode_side_short'
#  []
#  [h_side_long]
#    type = HeatTransferCoefficientAverage
#    variable = T
#    diffusivity = k
#    bulk_temperature = 300
#    boundary = 'electrode_side_long'
#  []
#[]

[BCs]
  [tip_temperature]
    type = NeumannBC
    variable = T
    value = 2060350000
    #value = 875.8333
    #value = 1063.3 
    #value = 1092.70017
    #value = 1212.62379
    #value = 1257.68526
    #value = 1326.56407
    #value = 1498.98188
    boundary = 'electrode_tip'
  []
  #[needle_top]
  #  type = DirichletBC
  #  variable = T
  #  value = 300
  #  boundary = 'needle_top'
  #[]

  #############################
  # Convective cooling
  #############################
  [T_cooling]
    type = ConvectiveHeatFluxBC
    variable = T
    T_infinity = 300
    heat_transfer_coefficient = h_mat
    boundary = 'ground right top_right top_left' 
  []

  ############################
  # Radiative cooling of electrode
  ############################
  #[bot_right]
  #  type = FunctionRadiativeBC
  #  variable = T
  #  Tinfinity = 300
  #  boundary = 'electrode_needle electrode_side_short electrode_side_long electrode_side'
  #  # htc/(stefan-boltzmann*4*T_inf^3)
  #  #emissivity_function = '3/(5.670367e-8*4*300*300*300)'
  #  emissivity_function = '1'
  #[]
  #[bot_right]
  #  type = MatRadiativeBC
  #  variable = T
  #  Tinfinity = 300
  #  boundary = 'electrode_needle electrode_side_short electrode_side_long electrode_side'
  #  #emissivity_function = '1'
  #[]

  
  ############################
  # Heat flux to water...?
  ############################
  #[heat_flux_water]
  #  type = NeumannBC
  #  variable = T
  #  #value = 0.03006   # 0.9 W * 0.0334 m
  #  #value = 300.6
  #  value = 900
  #  boundary = 'water_top' 
  #[]
  #[heat_flux_water_gaussian]
  #  type = FunctionNeumannBC
  #  variable = T
  #  function = q_gaussian
  #  boundary = 'water_top'
  #[]

  #############################
  # WATER REGION
  #############################
  [temperature_water_bulk]
    type = DirichletBC
    variable = T 
    value = 300
    boundary = 'water_bulk'
  []
[]

[Functions]
  [inlet_function]
    type = ParsedFunction
    vars = 'vmax k r2'
    vals = '-0.145 0.125 0.002'
    value = '(vmax / (1.0 - (sqrt(0.5*((1.0 - k^2.)/(log(1.0/k))))^2.0)*(1.0 - log(sqrt(0.5*((1.0 - k^2.)/(log(1.0/k))))^2.0))))*(1.0 - (x/r2)^2.0 - ((1.0 - k^2.)/(log(1.0/k))) * log(r2/x))'
  []

  [water_ic_func]
    type = ParsedFunction
    value = 'log(8.6949e23/6.022e23)'
  []

  [potential_bc_func]
    type = ParsedFunction
    value = 3
  []

  # a = 5.32473409e3
  # mu = 0
  # sigma = 1.60734825e-3 
  # Gaussian heat flux
  [q_gaussian]
    type = ParsedFunction
    #vars = 's m A q'
    #vals = '
    #vals = '0.008 0 159.57691216445178 256.802769519971354'
    #vals = '0.008 0 159.57691216445178 256.802769519971354'
    #value = 'q*A*exp(-(x-m)^2/(2*(s^2)))'

    vars = 'a mu sigma'
    vals = '5.32473409e3 0 1.60734825e-3' 
    value = 'a*exp(-(x-mu)^2/(2*(sigma^2)))'
    #val = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
  []
[]

[Materials]
  [heat_transfer_coefficient]
    type = GenericConstantMaterial
    prop_names = 'h_mat'
    prop_values = '10'
    block = '0 1'
  []
  [test_0]
    type = ADGenericConstantMaterial
    block = 0
    prop_names = 'diffp' 
    prop_values = '1'
  []
  [test_1]
    type = ADGenericConstantMaterial
    block = 1
    prop_names = 'diffw_p' 
    prop_values = '1'
  []
  # INS stuff
  [const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu cp k'
    prop_values = '1.784 3.77e-5 312.2 0.016'
  []

  [const_w]
    type = GenericConstantMaterial
    block = 1
    prop_names = 'rho mu cp k'
    prop_values = '977 8.9e-4 4182 0.6'
  []

  # Gas phase
  [density_argon]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.784' # kg m^-3
    block = 0
  []
  [thermal_conductivity_argon]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '0.016' # W m^-1 K^-1
    block = 0
  []
  [specific_heat_argon]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '312.2' # J kg^-1 K^-1
    block = 0
  []


  # Water phase
  [density_water]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '997' # kg m^-3
    block = 1
  []
  [thermal_conductivity_water]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '0.6' # W m^-1 K^-1
    block = 1
  []
  [specific_heat_water]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '4182' # J kg^-1 K^-1
    block = 1
  []

  # Dielectric
  [density_dielectric]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    #prop_values = '2500' # kg m^-3
    prop_values = '2230'
    block = 2
  []
  [thermal_conductivity_dielectric]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '1' # W m^-1 K^-1
    block = 2
  []
  [specific_heat_dielectric]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '800' # J kg^-1 K^-1
    block = 2
  []

  # Cathode
  [density_stainless]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '7850' # kg m^-3
    block = 3
  []
  #[thermal_conductivity_stainless]
  #  type = ADGenericConstantMaterial
  #  prop_names = 'thermal_conductivity'
  #  prop_values = '16.26' # W m^-1 K^-1
  #  block = 3
  #[]
  [thermal_conductivity_ss420]
    type = ADDerivativeParsedMaterial
    f_name = 'thermal_conductivity'
    args = 'T'
    function = '19.867+6.34e-3*T'
    block = 3
  []
  [specific_heat_stainless]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '502.1' # J kg^-1 K^-1
    block = 3
  []
  [emissivity_ss420]
    type = ADDerivativeParsedMaterial
    f_name = 'emissivity'
    args = 'T'
    function = '0.102+3.071e-5*T+(1.885e-7*T^2)'
    block = 3
  []
[]
