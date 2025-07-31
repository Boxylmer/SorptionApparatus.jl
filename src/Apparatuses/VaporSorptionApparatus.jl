struct VaporSorptionApparatus end

struct VaporSorptionSetup{T, CCIP, CCFP, SCFP, CCV, SCV, PMW, AA, AB, AC, PD, PM, VB, TSS, PLQMVT}
    temperature::T
    system_initial_pressures::Vector{CCIP}
    system_final_pressures::Vector{CCFP}
    charge_chamber_final_pressures::Vector{SCFP}
    
    charge_chamber_volume::CCV
    sampling_chamber_volume::SCV
    pen_mw::PMW
    antoine_a::AA
    antoine_b::AB
    antoine_c::AC

    polymer_density::PD
    polymer_mass::PM
    polymer_name::String
    num_steps::Int64

    num_beads::Int64
    vol_bead::VB

    transient_sorption_setup::TSS
    penetrant_liquid_phase_molar_volume::PLQMVT
end

VaporSorptionSetup(path::AbstractString; kwargs...) = readtemplate(VaporSorptionApparatus(), path; kwargs...)

struct VaporSorptionSystem{VSS, MSAS, ISO, VP, A, TFT, TSS}
    setup::VSS  # need to fix with actual VaporSorptionSystem type
    moles_sorbed_at_step::Vector{MSAS}
    isotherm::ISO
    vapor_pressure::VP
    activities::A
    thermo_factor_analysis::TFT
    transient_system::TSS
end # end VaporSorptionSystem

function VaporSorptionSystem(template_filepath::AbstractString; verbose=false, skip_transients=skip_transients)
    vss = VaporSorptionSetup(template_filepath; skip_transients=skip_transients)
    return VaporSorptionSystem(vss; verbose=verbose)
end

"""
Core math for a single vapor sorption step.
    units:
        # todo

"""

function moles_sorbed_during_step(vss::VaporSorptionSetup, step::Integer; transient_pressure=nothing)
    vc = vss.charge_chamber_volume * 0.01^3 # convert cm3 to m3
    empty_vs = vss.sampling_chamber_volume * 0.01^3 # convert cm3 to m3
    pol_vol = (vss.polymer_mass / vss.polymer_density) * 0.01^3 # convert cm3 to m3
    beads_vol = (vss.num_beads * vss.vol_bead) * 0.01^3 # convert cm3 to m3
    vs = empty_vs - pol_vol - beads_vol
    t = vss.temperature
    p_i = vss.system_initial_pressures[step]
    
    if isnothing(transient_pressure)
        p_f = vss.system_final_pressures[step]
    else
        p_f = transient_pressure    
    end

    if (step==1)
        pch_fm1 = 0.0
        p_fm1 = 0.0
    else
        pch_fm1 = vss.charge_chamber_final_pressures[step-1]
        p_fm1 = vss.system_final_pressures[step-1]
    end
    term_1 = p_i * vc          # initial moles in charge chamber
    term_2 = p_f * (vc + vs)   # final moles in the entire system
    term_3 = p_fm1 * (vc + vs) # final moles in the previous step
    term_4 = pch_fm1 * vc      # final moles in the charge chamber after closing the valve
    total_moles_sorbed_during_step = (term_1 - term_2 + term_3 - term_4) / t / MembraneBase.R_PA_M3_K_MOL
    return total_moles_sorbed_during_step
end
moles_sorbed_during_step(system::VaporSorptionSystem, step; transient_pressure=nothing) = moles_sorbed_during_step(system.setup, step; transient_pressure)

function dimensionless_mass_sorbed_during_step(vss::VaporSorptionSetup, step::Integer, equilibrium_moles_sorbed::Number, transient_pressure::Number)
    return moles_sorbed_during_step(vss, step; transient_pressure) / equilibrium_moles_sorbed  # if you multiply both by molecular weight, they cancel out
end

function VaporSorptionSystem(vss::VaporSorptionSetup; verbose=false)
    nps_at_each_step = [moles_sorbed_during_step(vss, i) for i in 1:vss.num_steps]
    nps = [sum(nps_at_each_step[1:i]) for i in 1:vss.num_steps]
    concs = [nps[i] / vss.polymer_mass * MembraneBase.CC_PER_MOL_STP * vss.polymer_density for i in 1:vss.num_steps]

    # define what equilibrium pressures are
    final_pressures = vss.system_final_pressures * MembraneBase.MPA_PER_PA
    
    # calculate vapor pressure
    vap_pres = MembraneBase.antoines_vapor_pressure(vss.antoine_a, vss.antoine_b, vss.antoine_c, vss.temperature) * MembraneBase.MPA_PER_BAR

    # approximate the activities
    activities = final_pressures / vap_pres

    # construct the isotherm
    isotherm = IsothermData(;
        partial_pressures_mpa=final_pressures, concentrations_cc=concs, activities=activities, temperature_k=vss.temperature,
        pen_mws_g_mol=[vss.pen_mw], rho_pol_g_cm3=vss.polymer_density
    )

    # calculate the thermodynamic factors
    tfa = ThermodynamicFactorAnalysis(isotherm)

    # deal with transients
    if isnothing(vss.transient_sorption_setup)
        transient_sorption_system = nothing
    else
        transient_sorption_system = TransientSorptionSystem(vss.transient_sorption_setup; verbose=verbose)
    end
        
    # construct the system struct
    system = VaporSorptionSystem(vss, nps, isotherm, vap_pres, activities, tfa, transient_sorption_system)
    return system
end

module VSAHelper
    ###### define templating system ######
    # inputs
    const default_sheet_title = "Vapor Sorption Input"
    const default_file_name = "Vapor Sorption Template (rename this!).xlsx"

    const overall_header, overall_val, overall_err = "D1", "E1", "F1"
    const t_header, t, t_err = "D2", "E2", "F2"
    const pres_apparatus_header, pres_apparatus_dummy, pres_apparatus_err = "D3", "E3", "F3"
    const vc_header, vc, vc_err = "D4", "E4", "F4"
    const vs_header, vs, vs_err = "D5", "E5", "F5"
    const num_beads_header, num_beads, num_beads_err = "D6", "E6", "F6"
    const vol_beads_header, vol_beads, vol_beads_err = "D7", "E7", "F7"

    const penetrant_property_header, penetrant_property_var, penetrant_property_err = "A1", "B1", "C1"
    const pen_name_header, pen_name, pen_name_dummy = "A2", "B2", "C2"
    const antoine_a_header, antoine_a, antoine_a_err  = "A3", "B3", "C3"
    const antoine_b_header, antoine_b, antoine_b_err  = "A4", "B4", "C4"
    const antoine_c_header, antoine_c, antoine_c_err = "A5", "B5", "C5"
    const mw_header, mw, mw_err =  "A6", "B6", "C6"

    const pol_property_header, polymer_property_var, polymer_property_err = "A7", "B7", "C7"
    const pol_name_header, pol_name, pol_name_dummy = "A8", "B8", "C8"
    const pol_dens_header, pol_dens, pol_dens_err = "A9", "B9", "C9"
    const pol_mass_header, pol_mass, pol_mass_err = "A10", "B10", "C10"

        # pressure data information
    const step_header, p_sys_in_header, p_sys_fin_header, p_ch_fin_header = "A11", "B11", "C11", "D11"
    const step_start, p_sys_in_start, p_sys_fin_start, p_ch_fin_start = "A12", "B12", "C12", "D12"
    const step_search_end = "D112"  # should be the index that would make a rectangle over the above constants
                                    # and encompass a large number of steps (more than anyone would ever use)
    const p_sys_in_col, p_sys_fin_col, p_ch_fin_col = 2, 3, 4  # in the above "rectangle" what columns correspond to what values?

    const vapor_pressure_header, vapor_pressure, vapor_pressure_err = "D9", "E9", "F9"

        # optional analyses
    const optional_analyses_header, optional_analyses_var, optional_analyses_err = "L1", "M1", "N1"
        # zimm-lundberg
    const zimm_lundberg_header, zim_lundberg_dummy_val, zim_lundberg_dummy_err = "L2", "M2", "N2"
    const liq_phase_mol_vol_header, liq_phase_mol_vol_val, liq_phase_mol_vol_err = "L3", "M3", "N3"


    # results output
    const pres_result_header, pres_result_err_header, pres_result_start, pres_result_err_start = "F11", "G11", "F12", "G12"
    const conc_result_header, conc_result_err_header, conc_result_start, conc_result_err_start = "H11", "I11", "H12", "I12"
    const conc_result_g_g_header, conc_result_g_g_err_header, conc_result_g_g_start, conc_result_g_g_err_start = "J11", "K11", "J12", "K12"
    const act_result_header, act_result_err_header, act_result_start, act_result_err_start = "L11", "M11", "L12", "M12"
    const mass_frac_result_header, mass_frac_result_err_header, mass_frac_result_start, mass_frac_result_err_start = "N11", "O11", "N12", "O12" # todo include mass fracs (ctrl+f the vars)
    const dual_mode_predictions_header, dual_mode_predictions_err_header, dual_mode_predictions_start, dual_mode_predictions_err_start = "P11", "Q11", "P12", "Q12"
    const gab_predictions_header, gab_predictions_err_header, gab_predictions_start, gab_predictions_err_start = "R11", "S11", "R12", "S12"
    const thermo_factor_header, thermo_factor_err_header, thermo_factor_start, thermo_factor_err_start = "T11", "U11", "T12", "U12"


    # parameters output
    const param_header, param_result_header, param_err_header, param_unit_header = "G1", "H1", "I1", "J1"
    const dual_mode_ch_header, dual_mode_ch_val, dual_mode_ch_err, dual_mode_ch_units= "G2", "H2", "I2", "J2"
    const dual_mode_b_header, dual_mode_b_val, dual_mode_b_err, dual_mode_b_units = "G3", "H3", "I3", "J3"
    const dual_mode_kd_header, dual_mode_kd_val, dual_mode_kd_err, dual_mode_kd_units = "G4", "H4", "I4", "J4"
    const gab_cp_header, gab_cp_val, gab_cp_err, gab_cp_units = "G5", "H5", "I5", "J5"
    const gab_k_header, gab_k_val, gab_k_err, gab_k_units = "G6", "H6", "I6", "J6"
    const gab_a_header, gab_a_val, gab_a_err, gab_a_units = "G7", "H7", "I7", "J7"
    
    function _add_vapor_sorption_results_headers(sheet)
        # add future results headers
        sheet[VSAHelper.pres_result_header], sheet[VSAHelper.pres_result_err_header] = "Pressure (MPa)", "Pres. err (MPa)"
        sheet[VSAHelper.conc_result_header], sheet[VSAHelper.conc_result_err_header] = "Concentration (CC(STP)/CC(Pol))", "Conc. err (CC(STP)/CC(Pol))"
        sheet[VSAHelper.conc_result_g_g_header], sheet[VSAHelper.conc_result_g_g_err_header] = "Concentration (g/g(Pol))", "Conc. err (g/g(Pol))"
        sheet[VSAHelper.act_result_header], sheet[VSAHelper.act_result_err_header] = "Activity", "act. err"
        sheet[VSAHelper.mass_frac_result_header], sheet[VSAHelper.mass_frac_result_err_header] = "Mass frac.", "Mass frac. err"
        sheet[VSAHelper.vapor_pressure_header], sheet[VSAHelper.vapor_pressure], sheet[VSAHelper.vapor_pressure_err] = "Vapor Pressure (MPa)", "calculated", "calculated"
            # params
        sheet[VSAHelper.param_header], sheet[VSAHelper.param_result_header], sheet[VSAHelper.param_err_header], sheet[VSAHelper.param_unit_header] = "Parameters", "Value", "Uncertainty (Jackknife)", "Units"
                # dual mode
        sheet[VSAHelper.dual_mode_ch_header], sheet[VSAHelper.dual_mode_ch_val], sheet[VSAHelper.dual_mode_ch_err], sheet[VSAHelper.dual_mode_ch_units] = "Dual Mode (CH')", "calculated", "calculated", "(CC/CC)"
        sheet[VSAHelper.dual_mode_b_header], sheet[VSAHelper.dual_mode_b_val], sheet[VSAHelper.dual_mode_b_err], sheet[VSAHelper.dual_mode_b_units] = "Dual Mode (b)", "calculated", "calculated", "1 / MPa"
        sheet[VSAHelper.dual_mode_kd_header], sheet[VSAHelper.dual_mode_kd_val], sheet[VSAHelper.dual_mode_kd_err], sheet[VSAHelper.dual_mode_kd_units] = "Dual Mode (kd)", "calculated", "calculated", "(CC/CC) / MPa"
        sheet[VSAHelper.dual_mode_predictions_header], sheet[VSAHelper.dual_mode_predictions_err_header] = "Dual Mode Pred (CC/CC)", "Dual Mode Pred Err (CC/CC)"
                # gab
        sheet[VSAHelper.gab_predictions_header], sheet[VSAHelper.gab_predictions_err_header] = "GAB Pred (CC/CC)", "GAB Pred Err (CC/CC)"
        sheet[VSAHelper.thermo_factor_header], sheet[VSAHelper.thermo_factor_err_header] = "Thermodynamic Factor α", "α err"
        sheet[VSAHelper.gab_cp_header], sheet[VSAHelper.gab_cp_val], sheet[VSAHelper.gab_cp_err], sheet[VSAHelper.gab_cp_units] = "GAB (Cp)", "calculated", "calculated", "(CC/CC)"
        sheet[VSAHelper.gab_k_header], sheet[VSAHelper.gab_k_val], sheet[VSAHelper.gab_k_err], sheet[VSAHelper.gab_k_units] = "GAB (k)", "calculated", "calculated", "Unitless"
        sheet[VSAHelper.gab_a_header], sheet[VSAHelper.gab_a_val], sheet[VSAHelper.gab_a_err], sheet[VSAHelper.gab_a_units] = "GAB (a)", "calculated", "calculated", "Unitless"

    end

    # Error analysis locations and data
    const default_error_analysis_name = "Uncertainty Analysis"
    const contribution_table_start_col = 6
    const contribution_table_start_row = 2
    const relative_antoine_a_contribution_col = 1
    const relative_antoine_b_contribution_col = 2
    const relative_antoine_c_contribution_col = 3
    const relative_pen_mw_contribution_col = 4
    const relative_pol_dens_contribution_col = 5
    const relative_pol_mass_contribution_col = 6
    const relative_temperature_contribution_col = 7
    const relative_vc_contribution_col = 8
    const relative_vs_contribution_col = 9
    const relative_sum_p_i_contribution_col = 10
    const relative_sum_p_f_contribution_col = 11
    const relative_sum_pcf_contribution_col = 12

end

# vapor sorption error analysis
struct VaporSorptionIsothermStepErrorAnalysis 
    # contributions (in percent)
    antoine_a; antoine_b; antoine_c
    pen_mw; pol_dens; pol_mass
    temperature; vc; vs
    p_i; p_f; p_cf
    sum_p_i; sum_p_f; sum_p_cf
end

function create_error_analysis(::VaporSorptionApparatus, system::VaporSorptionSystem)
    # println("Acquired analysis function")
    isotherm_concs = concentration(system.isotherm)
    # println("Isotherm ccs:", isotherm_concs)
    setup = system.setup
    analyses = []

    for (step_idx, conc) in enumerate(isotherm_concs)
        antoine_a = fractional_uncertainty_contribution(conc, setup.antoine_a)
        antoine_b = fractional_uncertainty_contribution(conc, setup.antoine_b)
        antoine_c = fractional_uncertainty_contribution(conc, setup.antoine_c)
        pen_mw = fractional_uncertainty_contribution(conc, setup.pen_mw)
        pol_dens = fractional_uncertainty_contribution(conc, setup.polymer_density)
        pol_mass = fractional_uncertainty_contribution(conc, setup.polymer_mass)
        temperature = fractional_uncertainty_contribution(conc, setup.temperature)
        vc = fractional_uncertainty_contribution(conc, setup.charge_chamber_volume)
        vs = fractional_uncertainty_contribution(conc, setup.sampling_chamber_volume)
        
        pis =  setup.system_initial_pressures[1:step_idx]
        pfs =  setup.system_final_pressures[1:step_idx]
        pcfs = setup.charge_chamber_final_pressures[1:step_idx]
        sum_p_i  = sum(fractional_uncertainty_contribution.(conc, pis))
        sum_p_f  = sum(fractional_uncertainty_contribution.(conc, pfs))
        sum_p_cf = sum(fractional_uncertainty_contribution.(conc, pcfs))
        p_i = fractional_uncertainty_contribution.(conc, pis[end])
        p_f = fractional_uncertainty_contribution.(conc, pfs[end])
        p_cf = fractional_uncertainty_contribution.(conc, pcfs[end])

        step_analysis = VaporSorptionIsothermStepErrorAnalysis(
            antoine_a, antoine_b, antoine_c, 
            pen_mw, pol_dens, pol_mass,
            temperature, vc, vs, 
            p_i, p_f, p_cf,
            sum_p_i, sum_p_f, sum_p_cf)

        push!(analyses, step_analysis)

    end
    return analyses

end
function write_error_analysis(::VaporSorptionApparatus, system::VaporSorptionSystem, xlsxfile)
    analyses = create_error_analysis(VaporSorptionApparatus(), system)
    sheet = XLSX.addsheet!(xlsxfile, VSAHelper.default_error_analysis_name)
    setup = system.setup

    # write out the initial things
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col, dim=1] = ["Variable", "Value", "Uncertainty", "Units"]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_antoine_a_contribution_col, dim=1] = ["Antoine A", setup.antoine_a.val, setup.antoine_a.err] 
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_antoine_b_contribution_col, dim=1] = ["Antoine B", setup.antoine_b.val, setup.antoine_b.err]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_antoine_c_contribution_col, dim=1] = ["Antoine C", setup.antoine_c.val, setup.antoine_c.err]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_pen_mw_contribution_col, dim=1] = ["Pen MW", setup.pen_mw.val, setup.pen_mw.err, "g/mol"]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_pol_dens_contribution_col, dim=1] = ["Pol Dens", setup.polymer_density.val, setup.polymer_density.err, "g/cm3"]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_pol_mass_contribution_col, dim=1] = ["Pol Mass", setup.polymer_mass.val, setup.polymer_mass.err, "g"]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_temperature_contribution_col, dim=1] = ["Temp", setup.temperature.val, setup.temperature.err, "K"]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_vc_contribution_col, dim=1] = ["Vc", setup.charge_chamber_volume.val, setup.charge_chamber_volume.err, "cm3"]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_vs_contribution_col, dim=1] = ["Vs", setup.sampling_chamber_volume.val, setup.sampling_chamber_volume.err, "cm3"]
    quick_and_dirty_pressure_err_fraction = string(setup.system_initial_pressures[1].err/setup.system_initial_pressures[1].val * 100) * "%"
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_sum_p_i_contribution_col, dim=1] = ["∑pi", "N/a", quick_and_dirty_pressure_err_fraction]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_sum_p_f_contribution_col, dim=1] = ["∑pf", "N/a", quick_and_dirty_pressure_err_fraction]
    sheet[VSAHelper.contribution_table_start_row, VSAHelper.contribution_table_start_col + VSAHelper.relative_sum_pcf_contribution_col, dim=1] = ["∑pcf", "N/a", quick_and_dirty_pressure_err_fraction] 

    # write out the uncertainty table
    for (analysis_idx, analysis) in enumerate(analyses)
        start_row = VSAHelper.contribution_table_start_row + 3
        start_col = VSAHelper.contribution_table_start_col
        sheet[start_row + analysis_idx, start_col] = analysis_idx
        
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_antoine_a_contribution_col] = analysis.antoine_a
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_antoine_b_contribution_col] = analysis.antoine_b
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_antoine_c_contribution_col] = analysis.antoine_c
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_pen_mw_contribution_col] = analysis.pen_mw
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_pol_dens_contribution_col] = analysis.pol_dens
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_pol_mass_contribution_col] = analysis.pol_mass
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_temperature_contribution_col] = analysis.temperature
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_vc_contribution_col] = analysis.vc
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_vs_contribution_col] = analysis.vs
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_sum_p_i_contribution_col] = analysis.sum_p_i
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_sum_p_f_contribution_col] = analysis.sum_p_f
        sheet[start_row + analysis_idx, start_col + VSAHelper.relative_sum_pcf_contribution_col] = analysis.sum_p_cf
    end
end

# generate, read, and process (read, evaluate, and write) templates
function generatetemplate(::VaporSorptionApparatus, filepath = VSAHelper.default_file_name)
    XLSX.openxlsx(filepath, mode="w"; enable_cache=true) do xf
        sheet = xf[1]
        XLSX.rename!(sheet, VSAHelper.default_sheet_title)
        sheet[VSAHelper.overall_header] = "Overall Properties"; sheet[VSAHelper.overall_val] = "Value"; sheet[VSAHelper.overall_err] = "Uncertainty"
        sheet[VSAHelper.t_header] = "Temperature (K)"
        sheet[VSAHelper.pres_apparatus_header] = "Pressure Transducer Err (e.g., 1% -> 0.01)"; sheet[VSAHelper.pres_apparatus_dummy] = "-----"
        sheet[VSAHelper.vc_header] = "Charge Chamber Vol. (cm3)"
        sheet[VSAHelper.vs_header] = "Sampling Chamber Vol. (empty) (cm3)"
        sheet[VSAHelper.num_beads_header] = "Num steel beads"
        sheet[VSAHelper.vol_beads_header] = "Volume of steel bead (cm3)"

        sheet[VSAHelper.penetrant_property_header] = "Penetrant Properties"; sheet[VSAHelper.penetrant_property_var] = "Value"; sheet[VSAHelper.penetrant_property_err] = "Uncertainty"
        sheet[VSAHelper.pen_name_header] = "Penetrant Name"; sheet[VSAHelper.pen_name_dummy] = ["-----"]
        
        sheet[VSAHelper.antoine_a_header] = "Antoine A (K -> Bar)"
        sheet[VSAHelper.antoine_b_header] = "Antoine B (K -> Bar)"
        sheet[VSAHelper.antoine_c_header] = "Antoine C (K -> Bar)"
        sheet[VSAHelper.mw_header] = "Pen. MW (g/mol)"
        
        sheet[VSAHelper.pol_property_header] = "Polymer Properties"; sheet[VSAHelper.polymer_property_var] = "Value"; sheet[VSAHelper.polymer_property_err] = "Uncertainty"
        sheet[VSAHelper.pol_name_header] = "Polymer Name"; sheet[VSAHelper.pol_name_dummy] = "-----"
        sheet[VSAHelper.pol_dens_header] = "Pol. dens (g/cm3)"
        sheet[VSAHelper.pol_mass_header] = "Pol. mass (g)"

        sheet[VSAHelper.step_header] = "Step"
        sheet[VSAHelper.p_sys_in_header] = "Initial System Pressure (Pa)" 
        sheet[VSAHelper.p_sys_fin_header] ="Final System Pressure (Pa)"
        sheet[VSAHelper.p_ch_fin_header] = "Final Charge Pressure (after valve is closed, Pa)"

        # add optional section
        sheet[VSAHelper.optional_analyses_header] = "Optional Analyses (add values to activate)"; sheet[VSAHelper.optional_analyses_var] = "Value"; sheet[VSAHelper.optional_analyses_err] = "σ"
        sheet[VSAHelper.zimm_lundberg_header] = "Zimm Lundberg Analysis"
        sheet[VSAHelper.zim_lundberg_dummy_val] = "---"
        sheet[VSAHelper.zim_lundberg_dummy_err] = "---"
        sheet[VSAHelper.liq_phase_mol_vol_header] = "Penetrant Liquid Phase Molar Volume (cm3/mol)"

        # add step values
        sheet[VSAHelper.step_start, dim=1] = collect(1:7)

        # add default values for some 
        sheet[VSAHelper.antoine_a_err] = 0; sheet[VSAHelper.antoine_b_err] = 0; sheet[VSAHelper.antoine_c_err] = 0; sheet[VSAHelper.mw_err] = 0;
        sheet[VSAHelper.t_err] = 0.5; sheet[VSAHelper.pres_apparatus_err] = 0.0025
        sheet[VSAHelper.vc], sheet[VSAHelper.vc_err] = 29.47783133, 0.098002
        sheet[VSAHelper.vs], sheet[VSAHelper.vs_err] = 7.630675074, 0.023176
        sheet[VSAHelper.num_beads], sheet[VSAHelper.num_beads_err] = 0, "---"
        sheet[VSAHelper.vol_beads], sheet[VSAHelper.vol_beads_err] = 0.13399839, 0
        
        VSAHelper._add_vapor_sorption_results_headers(sheet)
    end

    # add transient input
    generatetemplate(TransientSorptionApparatus(), filepath; standalone=false)
    return nothing
end

function readtemplate(::VaporSorptionApparatus, path::String; skip_transients=false)  # read a template into a sorption setup struct
    xf = XLSX.readxlsx(path)
    sheet = xf[VSAHelper.default_sheet_title]
    
    # _penetrant_name = string(sheet[VSAHelper.pen_name])
    if any([ismissing(sheet[VSAHelper.antoine_a]), ismissing(sheet[VSAHelper.antoine_b]), ismissing(sheet[VSAHelper.antoine_c]), ismissing(sheet[VSAHelper.mw])])
        throw(MissingException("Some required penetrant parameters are missing from the Excel sheet.")); 
    end
    
    _antoine_a = sheet[VSAHelper.antoine_a] ± (ismissing(sheet[VSAHelper.antoine_a_err]) ? 0 : sheet[VSAHelper.antoine_a_err])
    _antoine_b = sheet[VSAHelper.antoine_b] ± (ismissing(sheet[VSAHelper.antoine_b_err]) ? 0 : sheet[VSAHelper.antoine_b_err])
    _antoine_c = sheet[VSAHelper.antoine_c] ± (ismissing(sheet[VSAHelper.antoine_c_err]) ? 0 : sheet[VSAHelper.antoine_c_err])
    _mw = sheet[VSAHelper.mw] ± (ismissing(sheet[VSAHelper.mw_err]) ? 0 : sheet[VSAHelper.mw_err])
    
    if any([ismissing(sheet[VSAHelper.pol_dens]), ismissing(sheet[VSAHelper.pol_mass])])
        throw(MissingException("Some required polymer parameters are missing from the Excel sheet.")); 
    end
    _polymer_name = string(sheet[VSAHelper.pol_name])
    _pol_dens = sheet[VSAHelper.pol_dens] ± (ismissing(sheet[VSAHelper.pol_dens_err]) ? 0 : sheet[VSAHelper.pol_dens_err])
    _pol_mass = sheet[VSAHelper.pol_mass] ± (ismissing(sheet[VSAHelper.pol_mass_err]) ? 0 : sheet[VSAHelper.pol_mass_err])

    if any([ismissing(sheet[VSAHelper.t]), ismissing(sheet[VSAHelper.pres_apparatus_err]), ismissing(sheet[VSAHelper.vc]), ismissing(sheet[VSAHelper.vs])])
        throw(MissingException("Some required overall properties are missing from the Excel sheet.")); 
    end
    _t = sheet[VSAHelper.t] ± (ismissing(sheet[VSAHelper.t_err]) ? 0 : sheet[VSAHelper.t_err])
    _p_err = sheet[VSAHelper.pres_apparatus_err]
    _vc = sheet[VSAHelper.vc] ± (ismissing(sheet[VSAHelper.vc_err]) ? 0 : sheet[VSAHelper.vc_err][])
    _vs = sheet[VSAHelper.vs] ± (ismissing(sheet[VSAHelper.vs_err]) ? 0 : sheet[VSAHelper.vs_err])
    _nbeads = sheet[VSAHelper.num_beads]
    _vbead =  sheet[VSAHelper.vol_beads] ± (ismissing(sheet[VSAHelper.vol_beads_err]) ? 0 : sheet[VSAHelper.vol_beads_err])

    # read steps and format them correctly
    pressure_info_matrix = sheet[VSAHelper.step_start * ":" * VSAHelper.step_search_end]
    _initial_system_pressures = chop_vector_at_first_missing_value(pressure_info_matrix[:, VSAHelper.p_sys_in_col])
    _final_system_pressures = chop_vector_at_first_missing_value(pressure_info_matrix[:, VSAHelper.p_sys_fin_col])
    _final_charge_pressures = chop_vector_at_first_missing_value(pressure_info_matrix[:, VSAHelper.p_ch_fin_col])
    
    # ensure chopped steps are of the correct length
    n_isp = length(_initial_system_pressures)
    n_fsp = length(_final_system_pressures)
    n_fcp = length(_final_charge_pressures)
    if !((n_isp == n_fsp) && (n_isp == n_fcp || n_isp == n_fcp + 1))
        throw(MissingException("The number of steps specified didn't match up. 
        Was a value left blank, or was the sheet shifted away from the expected step start point?: " * VSAHelper.step_start *"?"))
    end  
 
    _initial_system_pressures = add_uncertainty_to_values(_initial_system_pressures, _p_err)
    _final_system_pressures = add_uncertainty_to_values(_final_system_pressures, _p_err)
    _final_charge_pressures = add_uncertainty_to_values(_final_charge_pressures, _p_err)
     
    _nsteps = length(_initial_system_pressures)
    transient_sorption_setup = nothing
    
    # copy of the setup to add to the transient apparatus (it'll grab some of the data to convert pressures to dimensionless sorption)
    temporary_setup = VaporSorptionSetup(
        _t, _initial_system_pressures, _final_system_pressures, _final_charge_pressures, _vc, _vs, _mw,
        _antoine_a, _antoine_b, _antoine_c, _pol_dens, _pol_mass, _polymer_name, _nsteps, _nbeads, _vbead, transient_sorption_setup,
        nothing) 
    
    if has_template(TransientSorptionApparatus(), xf) && !skip_transients
        transient_sorption_setup = readtemplate(TransientSorptionApparatus(), xf; apparatus_setup=temporary_setup) 
    end

    if ismissing(sheet[VSAHelper.liq_phase_mol_vol_val])
        pen_liq_phase_mol_vol = missing
    else
        pen_liq_phase_mol_vol = sheet[VSAHelper.liq_phase_mol_vol_val] ± (ismissing(sheet[VSAHelper.liq_phase_mol_vol_err]) ? 0 : sheet[VSAHelper.liq_phase_mol_vol_err])
    end
    
    setup = VaporSorptionSetup(
        _t, _initial_system_pressures, _final_system_pressures, _final_charge_pressures, _vc, _vs, _mw,
        _antoine_a, _antoine_b, _antoine_c, _pol_dens, _pol_mass, _polymer_name, _nsteps, _nbeads, _vbead, 
        transient_sorption_setup, pen_liq_phase_mol_vol
        )
    return setup
end

function processtemplate(::VaporSorptionApparatus, template_path::String, results_path::Union{String, Nothing}; sheet_name=VSAHelper.default_sheet_title, overwrite=false, verbose=false, skip_transients=false)  # read a template and write it out with options
    # load and calculate the system
    system = VaporSorptionSystem(template_path; verbose=verbose, skip_transients=skip_transients)
    if isnothing(results_path) return system end 
    isotherm = system.isotherm

    # copy the template (this is the only time we will not force by default)
    if template_path == results_path
        @error("The vapor sorption template and results file paths were the same.")
        return # not really dealing with this behavior yet, it would directly modify the file
    end
    cp(template_path, results_path; force=overwrite)

    # open the results file and start adding in the calculations done
    XLSX.openxlsx(results_path, mode="rw") do xf
    # xf = XLSX.open_xlsx_template(results_path)
        sheet = xf[sheet_name]
        # put the results headers back in, in case we added more columns
        VSAHelper._add_vapor_sorption_results_headers(sheet)

        sheet[VSAHelper.pres_result_start, dim=1] = [p.val for p in partial_pressures(isotherm, component=1)]
        sheet[VSAHelper.pres_result_err_start, dim=1] = [p.err for p in partial_pressures(isotherm, component=1)]
        sheet[VSAHelper.conc_result_start, dim=1] = [c.val for c in concentration(isotherm, component=1)]
        sheet[VSAHelper.conc_result_err_start, dim=1] = [c.err for c in concentration(isotherm, component=1)]
        sheet[VSAHelper.conc_result_g_g_start, dim=1] = [c.val for c in concentration(isotherm; component=1, gas_units=:g, pol_units=:g)]
        sheet[VSAHelper.conc_result_g_g_err_start, dim=1] = [c.err for c in concentration(isotherm; component=1, gas_units=:g, pol_units=:g)]

        sheet[VSAHelper.act_result_start, dim=1] = [a.val for a in system.activities]
        sheet[VSAHelper.act_result_err_start, dim=1] = [a.err for a in system.activities]

        sheet[VSAHelper.vapor_pressure], sheet[VSAHelper.vapor_pressure_err] = system.vapor_pressure.val, system.vapor_pressure.err
        
        # dual mode
        dmmodel = fit_model(DualMode(), isotherm; uncertainty_method=:JackKnife)
        sheet[VSAHelper.dual_mode_ch_val], sheet[VSAHelper.dual_mode_ch_err] = dmmodel.ch.val, dmmodel.ch.err
        sheet[VSAHelper.dual_mode_b_val], sheet[VSAHelper.dual_mode_b_err] = dmmodel.b.val, dmmodel.b.err
        sheet[VSAHelper.dual_mode_kd_val], sheet[VSAHelper.dual_mode_kd_err] = dmmodel.kd.val, dmmodel.kd.err
        
        dmpredictions = predict_concentration(dmmodel, partial_pressures(isotherm, component=1))
        sheet[VSAHelper.dual_mode_predictions_start, dim=1] = [dm.val for dm in dmpredictions]
        sheet[VSAHelper.dual_mode_predictions_err_start, dim=1] = [dm.err for dm in dmpredictions]

        # GAB
        gmmodel = fit_gab_model(system.activities, concentration(isotherm, component=1); uncertainty_method=:JackKnife)
        sheet[VSAHelper.gab_cp_val], sheet[VSAHelper.gab_cp_err] = gmmodel.cp.val, gmmodel.cp.err
        sheet[VSAHelper.gab_k_val], sheet[VSAHelper.gab_k_err] = gmmodel.k.val, gmmodel.k.err
        sheet[VSAHelper.gab_a_val], sheet[VSAHelper.gab_a_err] = gmmodel.a.val, gmmodel.a.err

        gmpredictions = a_predict_concentration(gmmodel, system.activities)
        sheet[VSAHelper.gab_predictions_start, dim=1] = [gm.val for gm in gmpredictions]
        sheet[VSAHelper.gab_predictions_err_start, dim=1] = [gm.err for gm in gmpredictions]

        # α factor

        αs = measurement.(system.thermo_factor_analysis.thermodynamic_factors)
        sheet[VSAHelper.thermo_factor_start, dim=1] = [α.val for α in αs]
        sheet[VSAHelper.thermo_factor_err_start, dim=1] = [α.err for α in αs]

        # write error analysis 
        write_error_analysis(VaporSorptionApparatus(), system, xf)

        # handle transient analysis if data is present
        if !isnothing(system.transient_system)
            write_transient_sorption_system_to_sheet(
                system.transient_system, 
                xf[TSAHelper.default_sheet_name] 
            )
            if system.setup.num_steps == system.transient_system.setup.num_steps
                write_kinetic_analysis(xf, isotherm, system.transient_system)
            end
        end

        # deal with zimm-lundberg
        if !ismissing(system.setup.penetrant_liquid_phase_molar_volume)
            zimm_lundberg_analysis = ZimmLundbergAnalysis(gmmodel, system.activities, system.setup.penetrant_liquid_phase_molar_volume)
            write_analysis(zimm_lundberg_analysis, xf)
        end

    end
    # XLSX.closefile(xf)
    return system

end

function processtemplate(::VaporSorptionApparatus, template_path::String; kwargs...)
    return processtemplate(VaporSorptionApparatus(), template_path, nothing; kwargs...)
end