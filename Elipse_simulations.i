###
### Rotman A. Criollo Manjarrez
### CO2SINk project
### Evolution of superhot reservoir from a magmatic intrusion
###
#
# PARAMETERS USED
# automatic_scaling = false
# report_peak_value = true
# perf_graph = false
# first_ICs = true

print_linear_residuals = true
checkpoint = true
exodus = true

BCtop = 'BCtop'
BCintrusion = 'intrusion_nodes'
Tintrusion = 823.15 #Kelvin
R = 8e9
w = 5.5e9

block_intrusion = 'BCintrusion'
block_reservoir = 'BCreservoir1 BCreservoir2'

# [Debug]
#   show_var_residual_norms = false
#   show_actions = false
# []

[Mesh]
  # uniform_refine = 1 
  [file]
    type = FileMeshGenerator
    file = '28_50x7x3_ellipse.e' #for shallow intrusion at 3km depth
	# file = '28_50x10x6_elipse.e' #for depth intrusion at 6km depth	
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 -9.81'
[]

[Variables]
  [pliquid]
    # scaling = 1E-5
  []
  [sgas]
    initial_condition = 0
  []
  [temp]
    scaling = 1E-6
  []
[]

[ICs]
  [pliquid]
    type = FunctionIC
    function = 'pres_hydrostatic'
    variable = 'pliquid'
    # ignore_uo_dependency = ${first_ICs}
  []  
  [temp]
    type = FunctionIC
    function = 'tempic'
    variable = 'temp'
    # ignore_uo_dependency = ${first_ICs}    
  []
[]

[Functions]
  [pres_hydrostatic]
    type = ParsedFunction
    expression = '101325 - 1005.88*9.81*(z)'
  []
  [tempic]
    type = ParsedFunction
    expression = '288.15 - (0.045*z)' 
  []
  [intrusion_func]
    type = ParsedFunction
    symbol_names = 'start_time time_mid Tini Tintrusion R w PI'
    symbol_values = '0.0 1e+10 662.924 ${Tintrusion} ${R} ${w} 3.14159265358979323846' #973.15 (700C) according with Violay et al., 2012
    # expression = (Tintrusion-Tini)/2*(tanh(PI*(t-R)/w)+1)+Tini
    expression = ${Tintrusion}
  []
[]

[BCs]
  [top_pressure]
    type = DirichletBC
    variable = 'pliquid'
    value = 101325 # Atmospheric pressure, Pa
    boundary = ${BCtop}
  []
  [top_temperature]
    type = DirichletBC
    variable = 'temp'
    value = 288.15 #Kelvin
    boundary = ${BCtop}
  []
  [pres_liq_lateral]
    type = PorousFlowPiecewiseLinearSink
    variable = 'pliquid'
    boundary = 'BCright BCleft BCfront BCback'
    pt_vals = '-1e9 1e9' # x coordinates defining g
    multipliers = '-1e9 1e9' # y coordinates defining g
    flux_function = 1e-11  # Variable C
    fluid_phase = 0
    save_in = 'pres_liq_lateral'
    use_mobility = true
    use_relperm = true
    use_displaced_mesh = false
    # mass_fraction_component = 0
  []
  [temp_liq_lateral]
    type = PorousFlowPiecewiseLinearSink
    variable = 'temp'
    boundary = 'BCright BCleft BCfront BCback'
    pt_vals = '-1e9 1e9' # x coordinates defining g
    multipliers = '-1e9 1e9' # y coordinates defining g
    flux_function = 1e-11 #'tempic' # Variable C
    fluid_phase = 0
    save_in = 'temp_liq_lateral'
    use_mobility = true
    use_relperm = true
    use_displaced_mesh = false
    # mass_fraction_component = 0
  []
  [pres_gas_lateral]
    type = PorousFlowPiecewiseLinearSink
    variable = 'pliquid'
    boundary = 'BCright BCleft BCfront BCback'
    pt_vals = '-1e9 1e9' # x coordinates defining g
    multipliers = '-1e9 1e9' # y coordinates defining g
    PT_shift = 'pgas' # 20E6   # BC pressure
    flux_function =  1e-11 # Variable C
    fluid_phase = 1
    save_in = 'pres_gas_lateral'
    use_mobility = true
    use_relperm = true
    use_displaced_mesh = false
    # mass_fraction_component = 0
  []
  [temp_gas_lateral]
    type = PorousFlowPiecewiseLinearSink
    variable = 'temp'
    boundary = 'BCright BCleft BCfront BCback'
    pt_vals = '-1e9 1e9' # x coordinates defining g
    multipliers = '-1e9 1e9' # y coordinates defining g
    PT_shift = 'temp' #288.15   # BC pressure
    flux_function = 1e-11 # Variable C
    fluid_phase = 1
    save_in = 'temp_gas_lateral'
    use_mobility = true
    use_relperm = true
    use_displaced_mesh = false
    # mass_fraction_component = 0
  []  
  [bottom_temperature_func]
    type = FunctionDirichletBC
    variable = 'temp'
    function = intrusion_func
    boundary = ${BCintrusion}
  []
[]

[AuxVariables]
  [pres_liq_lateral]
  []
  [temp_liq_lateral]
  []
  [temp_gas_lateral]
  []
  [pres_gas_lateral]
  []
  [massfrac_ph0_sp0]
    initial_condition = 1
  []
  [massfrac_ph1_sp0]
    initial_condition = 0
  []
  [massfrac_ph0_sp1]
  []
  [massfrac_ph1_sp1]
  []
  [pgas]
    family = MONOMIAL
    order = FIRST
  []
  [sliq]
    family = MONOMIAL
    order = FIRST
  []
  [density_liq]
    order = CONSTANT
    family = MONOMIAL
  []
  [density_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_liq]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [enthalpy_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [enthalpy_liq]
    order = CONSTANT
    family = MONOMIAL
  []
  [permeability]
    order = CONSTANT
    family = MONOMIAL
  []
  [effective_fluid_pressure]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [vel_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [vel_z]
    order = CONSTANT
    family = MONOMIAL
  []  
[]

[Kernels]
  [mass_liq_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = 'pliquid'
    multiply_by_density = true
  []
  [flux_liq]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = 'pliquid'
  []
  [mass_co2_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = 'sgas'
    multiply_by_density = true
  []
  [flux_co2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = 'sgas'
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = 'temp'
  []
  [advection]
    type = PorousFlowHeatAdvection
    variable = 'temp'
  []
  [conduction]
    type = PorousFlowHeatConduction
    variable = 'temp'
  []
  [heatcapacity]
    type = SpecificHeatConductionTimeDerivative
    variable = 'temp'
    density = 2700.0
  []
[]

[AuxKernels]
  [pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = 'pgas'
  []
  [sliq]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = 'sliq'
  []
  [enthalpy_gas]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 1
    variable = 'enthalpy_gas'
  []
  [enthalpy_liq]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 0
    variable = 'enthalpy_liq'
  []
  [density_liq]
    type = PorousFlowPropertyAux
    property = density
    phase = 0
    variable = 'density_liq'
  []
  [density_gas]
    type = PorousFlowPropertyAux
    property = density
    phase = 1
    variable = 'density_gas'
  []
  [viscosity_liq]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 0
    variable = 'viscosity_liq'
  []
  [viscosity_gas]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 1
    variable = 'viscosity_gas'
  []
  [permeability]
    type = ParsedAux
    coupled_variables = 'temp'
    variable = 'permeability'
    expression = 'if(temp <= 723.15,pow(10,log10(1e-15)),
              if(723.15 <= temp & temp < 773.15, pow(10,log10(1e-15) * (773.15 - temp) / (773.15 - 723.15) + log10(1e-16) * (temp - 723.15) / (773.15 - 723.15)),
              if(773.15 < temp & temp <= ${Tintrusion}, pow(10,log10(1e-16) * (${Tintrusion} - temp) / (${Tintrusion} - 773.15) + log10(1e-19) * (temp - 773.15) / (${Tintrusion} - 773.15)), 1e-19)))'
    execute_on ='INITIAL TIMESTEP_BEGIN'    
  []
  [effective_fluid_pressure]
    type = ParsedAux
    coupled_variables = 'pliquid pgas sliq sgas'
    expression = 'pliquid * sliq + pgas * sgas'
    variable = effective_fluid_pressure 
  []
  [vel_x]
    type = PorousFlowDarcyVelocityComponent
    variable = 'vel_x'
    component = x
    fluid_phase = 0
    # scaling = 1E15
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = PorousFlowDarcyVelocityComponent
    variable = 'vel_y'
    component = y
    fluid_phase = 0
    # scaling = 1E15
    execute_on = 'TIMESTEP_END'
  []
  [vel_z]
    type = PorousFlowDarcyVelocityComponent
    variable = 'vel_z'
    component = z
    fluid_phase = 0
    # scaling = 1E15
    execute_on = 'TIMESTEP_END'
  []  
[]

[UserObjects]
  [produced_mass_h2o]
    type = PorousFlowSumQuantity
  []
  [produced_mass_co2]
    type = PorousFlowSumQuantity
  []
  [produced_heat_h2o]
    type = PorousFlowSumQuantity
  []
  [produced_heat_co2]
    type = PorousFlowSumQuantity
  []
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'temp pliquid sgas'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
   type = PorousFlowCapillaryPressureVG
   alpha = 1e-5
   m = 0.5
  []
[]

[FluidProperties]
  [./h2o]
    type = Water97FluidProperties
  [../]
  # [./tabulated_water]
  #   type = TabulatedBicubicFluidProperties
  #   fp = h2o
  #   fluid_property_file = 'water_fluid_properties_Coolprop_old.csv'
  #   interpolated_properties = 'density enthalpy viscosity c'
  #   temperature_min = 273
  #   temperature_max = 1400
  #   pressure_min = 15000000
  #   pressure_max = 401000000
  # [../]
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 'temp'
  []
  [thermal_intrusion]
    type = HeatConductionMaterial
    block = ${block_intrusion}
    specific_heat = 850.0
    thermal_conductivity = 1.937  #this makes alpha 9.74e-5 m^2/s
          # above conductivity arbitrarily increased by 2 decades to make the
          #   object soak faster for the present purposes
    temp = 'temp'
  []
  [thermal_reservoir]
    type = HeatConductionMaterial
    block = ${block_reservoir}
    specific_heat = 850.0
    thermal_conductivity = 1.937  #this makes alpha 9.74e-5 m^2/s
    temp = 'temp'
  []   
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = 'pliquid'
    phase1_saturation = 'sgas'
    capillary_pressure = 'pc'
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  []
  [h2o_liquid]
    type = PorousFlowSingleComponentFluid
    fp = h2o
    phase = 0
  []
  [h2o_gas]
    type = PorousFlowSingleComponentFluid
    fp = h2o
    phase = 1
  []
  [relperm_liquid]
    type = PorousFlowRelativePermeabilityVG
    phase = 0
    m = 0.5 #VG exponent phase
    s_res = 0.2
    sum_s_res = 0.5
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityVG
    phase = 1
    m = 0.5  #VG exponent phase
    s_res = 0.2  # The residual saturation of the phase j. Must be between 0 and 1
    sum_s_res = 0.5  #Sum of residual saturations over all phases. Must be between 0 and 1
    # If true, use the van Genuchten form appropriate for a wetting (liquid) phase.
    # If false, use the non-wetting (gas) expression 
    wetting = false
  []
  [porosity_intrusion]
    type = PorousFlowPorosityConst # only the initial value of this is ever used
    porosity = 0.01
    block = ${block_intrusion}
  []
  [porosity_reservoir]
    type = PorousFlowPorosityConst # only the initial value of this is ever used
    porosity = 0.1
    block = ${block_reservoir}
  []  
  [biot_modulus]
    type = PorousFlowConstantBiotModulus #according with Parisio et al., 2019 Biot=0.5, K/Ks=0.5
    solid_bulk_compliance = 1E-10 #according with Parisio et al., 2019 #K
    fluid_bulk_modulus = 2E9  #according with  et al., 2019 #Ks
  []
  [permeability]
    type = PorousFlowPermeabilityTensorFromVar
    perm = permeability  
  []
  [thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    fluid_coefficient = 5E-6
    drained_coefficient = 2E-4
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1.5 0 0  0 1.5 0  0 0 1.5'  #'3 0 0  0 3. 0  0 0 3' #
    wet_thermal_conductivity = '3 0 0  0 3. 0  0 0 3'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    density = 2700.0 #according with Parisio et al., 2019
    specific_heat_capacity = 850.0 #according with Parisio et al., 2019
  []
  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
[]

[Preconditioning]
  active = 'mumps'
  [smp_faster]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -pc_factor_mat_solver_package -ksp_gmres_restart -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = ' gmres     lu       asm          superlu_dist                  1000               nonzero                   2               1e-1       1e-5       500'
  []
  [smp]
    type = SMP
    full = true
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    # petsc_options_value = 'lu        superlu_dist                  NONZERO                   2               1e1       1e-5        500'
    # petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -pc_factor_mat_solver_package -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      lu      asm           superlu_dist                  NONZERO                   2               1e1       1e-5        500'
  []
  [mumps]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    # petsc_options_value = 'lu        superlu_dist                  NONZERO                   2               1e1       1e-5        500'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
    petsc_options_value = 'gmres      lu       mumps                         NONZERO               1e1       1e1       500'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  # automatic_scaling = ${automatic_scaling}
  start_time = '-3.15576E+9' #-100y for stabilization before the intrusion is activated
  end_time = '3.15576e+12' #100kyr
  dtmax = '3.15576E+9' #100y
  dtmin = '1750.0'

  nl_rel_tol = '1e-6'
  nl_abs_tol = '1e-6'
  nl_max_its = '40' # max Non linear iterations before cutback is applied
  l_tol = '1e-6'
  l_abs_tol = '1e-6'
  l_max_its = '500'

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 5
    dt = 31557600 # 1y
    growth_factor = 2 # if iterations is less than nl max --> if iteration runs smoothly
    cutback_factor = 0.5 # if iteration exceeds nl max --> runs poorly
  []
[]

# [Dampers]
#   [./limit]
#     type = BoundingValueNodalDamper
#     variable = 'pliquid'
#     max_value = 4e8
#     min_value = 1e5
#   [../]
#   [./limit2]
#       type = BoundingValueNodalDamper
#       variable = 'temp'
#       max_value = 1400
#       min_value = 283
#   [../]
# []

# [Problem]
#   previous_nl_solution_required = true
# []

[Outputs]
  print_linear_residuals = ${print_linear_residuals}
  # perf_graph = ${perf_graph}
  exodus = ${exodus}
  checkpoint = ${checkpoint}
[]

[Controls]
  [./enable_volcan] #intrusion activation
    type = TimePeriod
    enable_objects = 'BCs::bottom_temperature_func'
    start_time = '0.0' #starts after 100y of stabilization
    end_time = '1.57788e+11' #5000y
    execute_on ='INITIAL TIMESTEP_BEGIN'
  [../]
[]
