struct GasSorptionApparatus end

struct GasSorptionSetup{T, CCIP, CCFP, SCFP, CCV, SCV, PTC, PPC, PO, PMW, PD, PM, NB, VB, MM, AD} 
    temperature::T                                      # kelvin
    charge_chamber_initial_pressures::Vector{CCIP}      # 
    charge_chamber_final_pressures::Vector{CCFP}        #
    sampling_chamber_final_pressures::Vector{SCFP}      #

    charge_chamber_volume::CCV                          #
    sampling_chamber_volume::SCV                        # 

    penetrant_name::String                              
    pen_tc::PTC                                         #                          
    pen_pc::PPC                                         #
    pen_omega::PO                                       #
    pen_mw::PMW                                         #

    polymer_density::PD                                 #
    polymer_mass::PM                                    #
    polymer_name::String                                #
    num_steps::Int64                                    # 

    n_beads::NB                                         # 
    vol_bead::VB                                        # 
    mesh_mass::MM                                       # 
    alum_dens::AD                                       # 
    # transient_sorption_setup::TSS  
end

GasSorptionSetup(path::AbstractString) = readtemplate(GasSorptionApparatus(), path)

struct GasSorptionSystem{GSS, MSAS, ISO, MVT}
    setup::GSS  # need to fix with actual GasSorptionSystem type todo
    moles_sorbed_at_step::AbstractVector{MSAS}
    # model::MembraneEOS.CubicModel
    isotherm::ISO
    molar_volumes::MVT

    # transient_system::TSS
end 

function moles_sorbed_during_step(gss::GasSorptionSetup, step::Integer, model::MembraneEOS.CubicModel; transient_pressure=nothing)
    polymer_volume = gss.polymer_mass / gss.polymer_density
    bead_volume = gss.vol_bead * gss.n_beads
    mesh_volume = gss.mesh_mass / gss.alum_dens
    real_sampling_chamber_volume = (gss.sampling_chamber_volume - polymer_volume - bead_volume - mesh_volume) * MembraneBase.L_PER_CC

    real_charge_chamber_volume = gss.charge_chamber_volume * MembraneBase.L_PER_CC    
    
    v_i_ch = volume(model, MembraneBase.ATM_PER_PA * gss.charge_chamber_initial_pressures[step], gss.temperature)
    n_i_ch = real_charge_chamber_volume / v_i_ch
    
    
    v_f_ch = volume(model, MembraneBase.ATM_PER_PA * gss.charge_chamber_final_pressures[step], gss.temperature)
    n_f_ch = real_charge_chamber_volume / v_f_ch
    
    if isnothing(transient_pressure)
        v_f_s = volume(model, MembraneBase.ATM_PER_PA * gss.sampling_chamber_final_pressures[step], gss.temperature)
    else
        v_f_s = volume(model, MembraneBase.ATM_PER_PA * transient_pressure, gss.temperature)
    end
    n_f_s = real_sampling_chamber_volume / v_f_s
    
    if step == 1
        n_fm1_s = 0
    else
        v_fm1_s = volume(model, MembraneBase.ATM_PER_PA * gss.sampling_chamber_final_pressures[step - 1], gss.temperature)
        n_fm1_s = real_sampling_chamber_volume / v_fm1_s
    end

    moles_sorbed = n_i_ch + n_fm1_s - n_f_ch - n_f_s # mol
    return moles_sorbed
end
moles_sorbed_during_step(system::GasSorptionSystem, step, model; kwargs...) = moles_sorbed_during_step(system.setup, step, model; kwargs...)  

# function dimensionless_mass_sorbed_during_step(gss::GasSorptionSetup, step::Integer, model::MembraneEOS.CubicModel, equilibrium_moles_sorbed::Number, transient_pressure::Number)
#     return moles_sorbed_during_step(gss, step, model; transient_pressure) / equilibrium_moles_sorbed  # if you multiply both by molecular weight, they cancel out
# end

function GasSorptionSystem(template_filepath::AbstractString;  kwargs...)
    gss = GasSorptionSetup(template_filepath)
    return GasSorptionSystem(gss;  kwargs...)
end

function GasSorptionSystem(gss::GasSorptionSetup)
    model = PR(CubicParameters(gss.pen_tc, gss.pen_pc, gss.pen_omega, gss.pen_mw))  # todo add options for other models?
    nps_at_each_step = [moles_sorbed_during_step(gss, i, model) for i in 1:gss.num_steps]
    nps = [sum(nps_at_each_step[1:i]) for i in 1:gss.num_steps]
    concs = [nps[i] / gss.polymer_mass * MembraneBase.CC_PER_MOL_STP * gss.polymer_density for i in 1:gss.num_steps]
    
    # define what equilibrium pressures are
    final_pressures = gss.sampling_chamber_final_pressures * MembraneBase.MPA_PER_PA
    molar_volumes = [volume(model, pres_mpa * MembraneBase.ATM_PER_MPA, gss.temperature) for pres_mpa in final_pressures]
    equilibrium_fugacities = [fugacity(model, pres ./ MembraneBase.MPA_PER_ATM, gss.temperature)[1] for pres in final_pressures] * MembraneBase.MPA_PER_ATM # MPa


    # create an isotherm
    isotherm = IsothermData(
        partial_pressures_mpa=final_pressures, concentrations_cc=concs, fugacities_mpa=equilibrium_fugacities, 
        temperature_k=gss.temperature, rho_pol_g_cm3=gss.polymer_density, pen_mws_g_mol=gss.pen_mw)
    system = GasSorptionSystem(gss, nps, isotherm, molar_volumes)
    return system
end

module GSAHelper
    ###### define templating system ######

    # input fields
    const default_sheet_title = "Gas Sorption Input"
    const default_file_name = "Gas Sorption Template (rename this!).xlsx"

    const overall_header, overall_val, overall_err = "D1", "E1", "F1"
    const t_header, t, t_err = "D2", "E2", "F2"
    const pres_apparatus_header, pres_apparatus_dummy, pres_apparatus_err = "D3", "E3", "F3"

    const vc_header, vc, vc_err = "D4", "E4", "F4"
    const vs_header, vs, vs_err = "D5", "E5", "F5"
    const nb_header, nb, nb_dummy = "D6", "E6", "F6"
    const nb_vol_header, nb_vol, nb_vol_err = "D7", "E7", "F7"
    const mesh_header, mesh_mass, mesh_mass_err = "D8", "E8", "F8"
    const alum_header, alum_dens, alum_dens_err = "D9", "E9", "F9"

    const penetrant_property_header, penetrant_property_var, penetrant_property_err = "A1", "B1", "C1"
    const pen_name_header, pen_name, pen_name_dummy = "A2", "B2", "C2"
    const omega_header, omega, omega_err = "A3", "B3", "C3"
    const pc_header, pc, pc_err  = "A4", "B4", "C4"
    const tc_header, tc, tc_err = "A5", "B5", "C5"
    const mw_header, mw, mw_err =  "A6", "B6", "C6"

    const pol_property_header, polymer_property_var, polymer_property_err = "A7", "B7", "C7"
    const pol_name_header, pol_name, pol_name_dummy = "A8", "B8", "C8"
    const pol_dens_header, pol_dens, pol_dens_err = "A9", "B9", "C9"
    const pol_mass_header, pol_mass, pol_mass_err = "A10", "B10", "C10"

        # pressure data information
    const step_header, p_ch_in_header, p_ch_fin_header, p_samp_header = "A11", "B11", "C11", "D11"
    const step_start, p_ch_in_start, p_ch_fin_start, p_samp_start = "A12", "B12", "C12", "D12"
    const step_search_end = "D112"  # should be the index that would make a rectangle over the above constants
                                    # and encompass a large number of steps (more than anyone would ever use)
    const p_ch_in_col, p_ch_fin_col, p_samp_col = 2, 3, 4  # in the above "rectangle" what columns correspond to what values?
   

    # fitted parameters headers
    const param_header, param_result_header, param_err_header, param_unit_header = "G1", "H1", "I1", "J1"
    const dual_mode_ch_header, dual_mode_ch_val, dual_mode_ch_err, dual_mode_ch_units= "G2", "H2", "I2", "J2"
    const dual_mode_b_header, dual_mode_b_val, dual_mode_b_err, dual_mode_b_units = "G3", "H3", "I3", "J3"
    const dual_mode_kd_header, dual_mode_kd_val, dual_mode_kd_err, dual_mode_kd_units = "G4", "H4", "I4", "J4"

    # result headers
    const pres_result_header, pres_result_err_header, pres_result_start, pres_result_err_start = "F11", "G11", "F12", "G12"
    const conc_result_header, conc_result_err_header, conc_result_start, conc_result_err_start = "H11", "I11", "H12", "I12"
    const conc_result_g_g_header, conc_result_g_g_err_header, conc_result_g_g_start, conc_result_g_g_err_start = "J11", "K11", "J12", "K12"
    const fug_result_header, fug_result_err_header, fug_result_start, fug_result_err_start = "L11", "M11", "L12", "M12"
    const molar_vol_result_header, molar_vol_result_err_header, molar_vol_result_start, molar_vol_result_err_start = "N11", "O11", "N12", "O12"
    const mass_frac_result_header, mass_frac_result_err_header, mass_frac_result_start, mass_frac_result_err_start = "P11", "Q11", "P12", "Q12"
    const dual_mode_predictions_header, dual_mode_predictions_err_header, dual_mode_predictions_start, dual_mode_predictions_err_start = "R11", "S11", "R12", "S12"

    function add_headers_to_sheet(sheet)
        # input headers
        sheet[overall_header] = "Overall Properties"; sheet[overall_val] = "Value"; sheet[overall_err] = "Uncertainty"
        sheet[t_header] = "Temperature (K)"
        sheet[pres_apparatus_header] = "Pressure Transducer Err (e.g., 1% -> 0.01)"; sheet[pres_apparatus_dummy] = "-----"
        sheet[vc_header] = "Charge Chamber Vol. (cm3)"
        sheet[vs_header] = "Sampling Chamber Vol. (empty) (cm3)"
        sheet[nb_header] = "Num. beads"; sheet[nb_dummy] = "-----"; sheet[nb_vol_header] = "Volume of steel bead (cm3)"
        sheet[mesh_header] = "Mesh mass (g)"; sheet[alum_header] = "ρ mesh (g/cm3)"
        
        sheet[penetrant_property_header] = "Penetrant Properties"; sheet[penetrant_property_var] = "Value"; sheet[penetrant_property_err] = "Uncertainty"
        sheet[pen_name_header] = "Penetrant Name"; sheet[pen_name_dummy] = ["-----"]
        sheet[omega_header] = "ω"
        sheet[pc_header] = "Pc (atm)"
        sheet[tc_header] = "Tc (K)"
        sheet[mw_header] = "Pen. MW (g/mol)"
        
        sheet[pol_property_header] = "Polymer Properties"; sheet[polymer_property_var] = "Value"; sheet[polymer_property_err] = "Uncertainty"
        sheet[pol_name_header] = "Polymer Name"; sheet[pol_name_dummy] = "-----"
        sheet[pol_dens_header] = "Pol. dens (g/cm3)"
        sheet[pol_mass_header] = "Pol. mass (g)"

        sheet[step_header] = "Step"; sheet[p_ch_in_header] = "Init. Charge P (Pa)"; sheet[p_ch_fin_header] ="Fin. Charge P (Pa)"; sheet[p_samp_header] = "Fin. Sampling P (Pa)"

        # fitted parameter headers
        sheet[GSAHelper.param_header], sheet[GSAHelper.param_result_header], sheet[GSAHelper.param_err_header], sheet[GSAHelper.param_unit_header] = "Parameters", "Value", "Uncertainty (Jackknife)", "Units"
            # dual mode
        sheet[GSAHelper.dual_mode_ch_header], sheet[GSAHelper.dual_mode_ch_val], sheet[GSAHelper.dual_mode_ch_err], sheet[GSAHelper.dual_mode_ch_units] = "Dual Mode (CH')", "calculated", "calculated", "(CC/CC)"
        sheet[GSAHelper.dual_mode_b_header], sheet[GSAHelper.dual_mode_b_val], sheet[GSAHelper.dual_mode_b_err], sheet[GSAHelper.dual_mode_b_units] = "Dual Mode (b)", "calculated", "calculated", "1 / MPa"
        sheet[GSAHelper.dual_mode_kd_header], sheet[GSAHelper.dual_mode_kd_val], sheet[GSAHelper.dual_mode_kd_err], sheet[GSAHelper.dual_mode_kd_units] = "Dual Mode (kd)", "calculated", "calculated", "(CC/CC) / MPa"
            # todo NELF
        
        # result headers
        sheet[GSAHelper.pres_result_header], sheet[GSAHelper.pres_result_err_header] = "Pressure (MPa)", "Pres. err (MPa)"
        sheet[GSAHelper.conc_result_header], sheet[GSAHelper.conc_result_err_header] = "Concentration (CC(STP)/CC(Pol))", "Conc. err (CC(STP)/CC(Pol))"
        sheet[GSAHelper.conc_result_g_g_header], sheet[GSAHelper.conc_result_g_g_err_header] = "Concentration (g/g(Pol))", "Conc. err (g/g(Pol))"
        sheet[GSAHelper.fug_result_header], sheet[GSAHelper.fug_result_err_header] = "Fugacity (MPa)", "Fug. err (MPa)"
        sheet[GSAHelper.molar_vol_result_header], sheet[GSAHelper.molar_vol_result_err_header] = "Molar Vol (L/mol)", "M. Vol. err (L/mol)"
        sheet[GSAHelper.mass_frac_result_header], sheet[GSAHelper.mass_frac_result_err_header] = "Mass frac.", "Mass frac. err"
        sheet[GSAHelper.dual_mode_predictions_header], sheet[GSAHelper.dual_mode_predictions_err_header] = "Dual Mode Pred (CC/CC)", "Dual Mode Pred Err (CC/CC)"

    end

 
    # Error analysis locations and data
    const default_error_analysis_name = "Uncertainty Analysis"
    const contribution_table_start_col = 6
    const contribution_table_start_row = 2
    const relative_omega_contribution_col = 1
    const relative_pc_contribution_col = 2
    const relative_tc_contribution_col = 3
    const relative_pen_mw_contribution_col = 4
    const relative_pol_dens_contribution_col = 5
    const relative_pol_mass_contribution_col = 6
    const relative_temperature_contribution_col = 7
    const relative_vc_contribution_col = 8
    const relative_vs_contribution_col = 9
    const relative_vol_bead_contribution_col = 10
    const relative_mesh_mass_contribution_col = 11
    const relative_alum_dens_contribution_col = 12
    const relative_sum_p_ch_i_contribution_col = 13
    const relative_sum_p_ch_f_contribution_col = 14
    const relative_sum_p_s_f_contribution_col = 15

end

# gas sorption error analysis
struct GasSorptionIsothermStepErrorAnalysis
    omega; pc; tc
    pen_mw; pol_dens; pol_mass
    temperature; vc; vs
    vol_bead; mesh_mass; alum_dens
    sum_p_ch_i; sum_p_ch_f; sum_p_s_f
    p_ch_i; p_ch_f; p_s_f 
end

function create_error_analysis(::GasSorptionApparatus, system::GasSorptionSystem)
    setup = system.setup
    # component mode will use finite differencing via the @uncertain macro to fully determine the uncertainty_components of the isotherm
    isotherm_concs = concentration(system.isotherm)
    analyses = []

    for (step_idx, conc) in enumerate(isotherm_concs)

        omega = fractional_uncertainty_contribution(conc, system.setup.pen_omega)
        pc = fractional_uncertainty_contribution(conc, system.setup.pen_pc)
        tc = fractional_uncertainty_contribution(conc, system.setup.pen_tc)
        pen_mw = fractional_uncertainty_contribution(conc, system.setup.pen_mw)
        pol_dens = fractional_uncertainty_contribution(conc, setup.polymer_density)
        pol_mass = fractional_uncertainty_contribution(conc, setup.polymer_mass)
        temperature = fractional_uncertainty_contribution(conc, setup.temperature)
        vs = fractional_uncertainty_contribution(conc, setup.sampling_chamber_volume)
        vc = fractional_uncertainty_contribution(conc, setup.charge_chamber_volume)
        
        vol_bead = fractional_uncertainty_contribution(conc, setup.vol_bead)
        mesh_mass = fractional_uncertainty_contribution(conc, setup.mesh_mass)
        alum_dens = fractional_uncertainty_contribution(conc, setup.alum_dens)

        p_ch_is =  setup.charge_chamber_initial_pressures[1:step_idx]
        p_ch_fs =  setup.charge_chamber_final_pressures[1:step_idx]
        p_s_fs  =  setup.sampling_chamber_final_pressures[1:step_idx]
        sum_p_ch_i  = sum(fractional_uncertainty_contribution.(conc, p_ch_is))
        sum_p_ch_f  = sum(fractional_uncertainty_contribution.(conc, p_ch_fs))
        sum_p_s_f   = sum(fractional_uncertainty_contribution.(conc, p_s_fs))
        p_ch_i = fractional_uncertainty_contribution.(conc, p_ch_is[end])
        p_ch_f = fractional_uncertainty_contribution.(conc, p_ch_fs[end])
        p_s_f =  fractional_uncertainty_contribution.(conc, p_s_fs[end])

        step_analysis = GasSorptionIsothermStepErrorAnalysis(
            omega, pc, tc,
            pen_mw, pol_dens, pol_mass,
            temperature, vc, vs,
            vol_bead, mesh_mass, alum_dens,
            sum_p_ch_i, sum_p_ch_f, sum_p_s_f,
            p_ch_i, p_ch_f, p_s_f)

        push!(analyses, step_analysis)
    end
    return analyses
end

function write_error_analysis(::GasSorptionApparatus, system::GasSorptionSystem, xlsxfile)
    analyses = create_error_analysis(GasSorptionApparatus(), system)
    sheet = XLSX.addsheet!(xlsxfile, GSAHelper.default_error_analysis_name)
    setup = system.setup

    # write out the initial things
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col, dim=1] = ["Variable", "Value", "Uncertainty", "Units"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_omega_contribution_col, dim=1] = ["ω", setup.pen_omega.val, setup.pen_omega.err] 
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_pc_contribution_col, dim=1] = ["Pc", setup.pen_pc.val, setup.pen_pc.err, "atm"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_tc_contribution_col, dim=1] = ["Tc", setup.pen_tc.val, setup.pen_tc.err, "K"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_pen_mw_contribution_col, dim=1] = ["Pen MW", setup.pen_mw.val, setup.pen_mw.err, "g/mol"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_pol_dens_contribution_col, dim=1] = ["Pol Dens", setup.polymer_density.val, setup.polymer_density.err, "g/cm3"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_pol_mass_contribution_col, dim=1] = ["Pol Mass", setup.polymer_mass.val, setup.polymer_mass.err, "g"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_temperature_contribution_col, dim=1] = ["Temp", setup.temperature.val, setup.temperature.err, "K"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_vc_contribution_col, dim=1] = ["Vc", setup.charge_chamber_volume.val, setup.charge_chamber_volume.err, "cm3"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_vs_contribution_col, dim=1] = ["Vs", setup.sampling_chamber_volume.val, setup.sampling_chamber_volume.err, "cm3"]
    
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_vol_bead_contribution_col, dim=1] = ["Bead volume", setup.vol_bead.val, setup.vol_bead.err, "cm3"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_mesh_mass_contribution_col, dim=1] = ["Mesh mass", setup.mesh_mass.val, setup.mesh_mass.err, "g"]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_alum_dens_contribution_col, dim=1] = ["Mesh ρ", setup.alum_dens.val, setup.alum_dens.err, "g/cm3"]
    

    quick_and_dirty_pressure_err_fraction = string(setup.charge_chamber_initial_pressures[1].err/setup.charge_chamber_initial_pressures[1].val * 100) * "%"
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_sum_p_ch_i_contribution_col, dim=1] = ["∑pi charge", "N/a", quick_and_dirty_pressure_err_fraction]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_sum_p_ch_f_contribution_col, dim=1] = ["∑pf charge", "N/a", quick_and_dirty_pressure_err_fraction]
    sheet[GSAHelper.contribution_table_start_row, GSAHelper.contribution_table_start_col + GSAHelper.relative_sum_p_s_f_contribution_col, dim=1] =  ["∑pf sample", "N/a", quick_and_dirty_pressure_err_fraction] 

    # write out the uncertainty table
    for (analysis_idx, analysis) in enumerate(analyses)
        start_row = GSAHelper.contribution_table_start_row + 3
        start_col = GSAHelper.contribution_table_start_col
        sheet[start_row + analysis_idx, start_col] = analysis_idx
        
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_omega_contribution_col] = analysis.omega
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_pc_contribution_col] = analysis.pc
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_tc_contribution_col] = analysis.tc
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_pen_mw_contribution_col] = analysis.pen_mw
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_pol_dens_contribution_col] = analysis.pol_dens
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_pol_mass_contribution_col] = analysis.pol_mass
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_temperature_contribution_col] = analysis.temperature
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_vc_contribution_col] = analysis.vc
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_vs_contribution_col] = analysis.vs

        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_vol_bead_contribution_col] = analysis.vol_bead
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_mesh_mass_contribution_col] = analysis.mesh_mass
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_alum_dens_contribution_col] = analysis.alum_dens

        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_sum_p_ch_i_contribution_col] = analysis.sum_p_ch_i
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_sum_p_ch_f_contribution_col] = analysis.sum_p_ch_f
        sheet[start_row + analysis_idx, start_col + GSAHelper.relative_sum_p_s_f_contribution_col] = analysis.sum_p_s_f
    end
end

function generatetemplate(::GasSorptionApparatus, filepath = GSAHelper.default_file_name)
    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, GSAHelper.default_sheet_title)
        
        GSAHelper.add_headers_to_sheet(sheet)

        # add step values
        sheet[GSAHelper.step_start, dim=1] = collect(1:7)

        # add default values upon sheet creation
        sheet[GSAHelper.omega_err], sheet[GSAHelper.pc_err], sheet[GSAHelper.tc_err], sheet[GSAHelper.mw_err] = 0, 0, 0, 0

        sheet[GSAHelper.alum_dens] = 2.71; sheet[GSAHelper.alum_dens_err] = 0.00
        sheet[GSAHelper.t_err] = 0.5
        sheet[GSAHelper.pres_apparatus_err] = 0.0005
        
        sheet[GSAHelper.nb], sheet[GSAHelper.nb_vol], sheet[GSAHelper.nb_vol_err] = 0, 0.13399839, 0
        sheet[GSAHelper.mesh_mass], sheet[GSAHelper.mesh_mass_err] = 0, 0

        # sheet[GSAHelper.vs], sheet[GSAHelper.vs_err] = 15.707, 0.06 # todo add values
        # sheet[GSAHelper.vc], sheet[GSAHelper.vc_err] = 18.442, 0.04

    end 
    # # add transient template
    # generatetemplate(TransientSorptionApparatus(), filepath; standalone=false)
end

function readtemplate(::GasSorptionApparatus, path::AbstractString; verify::Bool=true)
    xf = XLSX.readxlsx(path)
    sheet = xf[1]
    if verify
        if any([ismissing(sheet[GSAHelper.omega]), ismissing(sheet[GSAHelper.pc]), ismissing(sheet[GSAHelper.tc]), ismissing(sheet[GSAHelper.mw])])
            throw(MissingException("Some required penetrant parameters are missing from the Excel sheet.")); 
        end
        if any([ismissing(sheet[GSAHelper.pol_dens])])
            throw(MissingException("Some required polymer parameters are missing from the Excel sheet.")); 
        end
        if any([ismissing(sheet[GSAHelper.t]), ismissing(sheet[GSAHelper.pres_apparatus_err]), ismissing(sheet[GSAHelper.vc]), ismissing(sheet[GSAHelper.vs]),
            ismissing(sheet[GSAHelper.nb_vol]), ismissing(sheet[GSAHelper.mesh_mass]), ismissing(sheet[GSAHelper.alum_dens])])
            throw(MissingException("Some required overall properties are missing from the Excel sheet.")); 
        end
    end
    _penetrant_name = string(sheet[GSAHelper.pen_name])

    _tc = sheet[GSAHelper.tc] ± (ismissing(sheet[GSAHelper.tc_err]) ? 0 : sheet[GSAHelper.tc_err])
    _pc = sheet[GSAHelper.pc] ± (ismissing(sheet[GSAHelper.pc_err]) ? 0 : sheet[GSAHelper.pc_err])
    _omega = sheet[GSAHelper.omega] ± (ismissing(sheet[GSAHelper.omega_err]) ? 0 : sheet[GSAHelper.omega_err])
    _mw = sheet[GSAHelper.mw] ± (ismissing(sheet[GSAHelper.mw_err]) ? 0 : sheet[GSAHelper.mw_err])
 
    _polymer_name = string(sheet[GSAHelper.pol_name])
    _pol_dens = sheet[GSAHelper.pol_dens] ± (ismissing(sheet[GSAHelper.pol_dens_err]) ? 0 : sheet[GSAHelper.pol_dens_err])
    _pol_mass = sheet[GSAHelper.pol_mass] ± (ismissing(sheet[GSAHelper.pol_mass_err]) ? 0 : sheet[GSAHelper.pol_mass_err])

  
    _t = sheet[GSAHelper.t] ± (ismissing(sheet[GSAHelper.t_err]) ? 0 : sheet[GSAHelper.t_err])
    _p_err = sheet[GSAHelper.pres_apparatus_err]
    _vc = sheet[GSAHelper.vc] ± (ismissing(sheet[GSAHelper.vc_err]) ? 0 : sheet[GSAHelper.vc_err][])
    _vs = sheet[GSAHelper.vs] ± (ismissing(sheet[GSAHelper.vs_err]) ? 0 : sheet[GSAHelper.vs_err])
    _nb = ismissing(sheet[GSAHelper.nb]) ? 0 : sheet[GSAHelper.nb]
    _nb_vol = sheet[GSAHelper.nb_vol] ± (ismissing(sheet[GSAHelper.nb_vol_err]) ? 0 : sheet[GSAHelper.nb_vol_err])
    _mesh_mass = sheet[GSAHelper.mesh_mass] ± (ismissing(sheet[GSAHelper.mesh_mass_err]) ? 0 : sheet[GSAHelper.mesh_mass_err]) 
    _alum_dens = sheet[GSAHelper.alum_dens] ± (ismissing(sheet[GSAHelper.alum_dens_err]) ? 0 : sheet[GSAHelper.alum_dens_err])

    # read steps and format them correctly
    pressure_info_matrix = sheet[GSAHelper.step_start * ":" * GSAHelper.step_search_end]
    _initial_charge_pressures = chop_vector_at_first_missing_value(pressure_info_matrix[:, GSAHelper.p_ch_in_col])
    _final_charge_pressures = chop_vector_at_first_missing_value(pressure_info_matrix[:, GSAHelper.p_ch_fin_col])
    _final_sampling_pressure = chop_vector_at_first_missing_value(pressure_info_matrix[:, GSAHelper.p_samp_col])
    
    # ensure chopped steps are of the correct length
    # TODO: verify that we will always need information at every step, etc
    if verify
        if !(length(_initial_charge_pressures) == length(_final_charge_pressures) == length(_final_sampling_pressure))
            throw(MissingException("The number of steps specified didn't match up. 
            Was a value left blank or was the sheet shifted away from the expected step start point: " * GSAHelper.step_start *"?"))
        end
    end

    _initial_charge_pressures = add_uncertainty_to_values(_initial_charge_pressures, _p_err)
    _final_charge_pressures = add_uncertainty_to_values(_final_charge_pressures, _p_err)
    _final_sampling_pressure = add_uncertainty_to_values(_final_sampling_pressure, _p_err)
     
    _nsteps = length(_initial_charge_pressures)

    # transient_sorption_setup = nothing
    
    # # copy of the setup to add to the transient apparatus (it'll grab some of the data to convert pressures to dimensionless sorption)
    # temporary_setup = GasSorptionSetup(_t, _initial_charge_pressures, _final_charge_pressures, _final_sampling_pressure, 
    #     _vc, _vs, _penetrant, _model, _pol_dens, _pol_mass, _polymer_name, _nsteps, _nb, _nb_vol, _mesh_mass, _alum_dens, transient_sorption_setup) 
        
    # if is_template_valid(TransientSorptionApparatus(), path; apparatus_setup=temporary_setup)
    #     transient_sorption_setup = readtemplate(TransientSorptionApparatus(), path; apparatus_setup=temporary_setup)
    # end

    setup = GasSorptionSetup(_t, _initial_charge_pressures, _final_charge_pressures, _final_sampling_pressure, 
            _vc, _vs, _penetrant_name, _tc, _pc, _omega, _mw, _pol_dens, _pol_mass, _polymer_name, _nsteps, _nb, _nb_vol, _mesh_mass, _alum_dens) 
    return setup
end

function processtemplate(::GasSorptionApparatus, template_path::String, results_path::Union{String, Nothing}; sheet_name=GSAHelper.default_sheet_title, overwrite=false)  # read a template and write it out with options

    # load and calculate the system
    system = GasSorptionSystem(template_path)
    if isnothing(results_path) return system end 
    isotherm = system.isotherm

    # copy the template (this is the only time we will not force by default)
    if template_path == results_path
        @warn("The gas sorption template and results file paths were the same")
        return # not really dealing with this behavior yet, it would directly modify the file
    end
    cp(template_path, results_path; force=overwrite)

    # open the results file and start adding in the calculations done
    XLSX.openxlsx(results_path, mode="rw") do xf
        sheet = xf[sheet_name]
        # put the results headers back in, in case we added more columns or otherwise changed the way we process the template
        GSAHelper.add_headers_to_sheet(sheet)

        # output the results of the isotherm
        sheet[GSAHelper.pres_result_start, dim=1] = [p.val for p in partial_pressures(isotherm, component=1)]
        sheet[GSAHelper.pres_result_err_start, dim=1] = [p.err for p in partial_pressures(isotherm, component=1)]
        sheet[GSAHelper.conc_result_start, dim=1] = [c.val for c in concentration(isotherm, component=1)]
        sheet[GSAHelper.conc_result_err_start, dim=1] = [c.err for c in concentration(isotherm, component=1)]
        sheet[GSAHelper.conc_result_g_g_start, dim=1] = [c.val for c in concentration(isotherm; component=1, gas_units=:g, pol_units=:g)]
        sheet[GSAHelper.conc_result_g_g_err_start, dim=1] = [c.err for c in concentration(isotherm; component=1, gas_units=:g, pol_units=:g)]
        sheet[GSAHelper.fug_result_start, dim=1] = [f.val for f in fugacities(isotherm; component=1)]
        sheet[GSAHelper.fug_result_err_start, dim=1] = [f.err for f in fugacities(isotherm; component=1)]
        sheet[GSAHelper.molar_vol_result_start, dim=1] = [molar_vol.val for molar_vol in system.molar_volumes]
        sheet[GSAHelper.molar_vol_result_err_start, dim=1] = [molar_vol.err for molar_vol in system.molar_volumes]
        sheet[GSAHelper.mass_frac_result_start, dim=1] = [mf.val for mf in penetrant_mass_fractions(isotherm; component=1)]
        sheet[GSAHelper.mass_frac_result_err_start, dim=1] = [mf.err for mf in penetrant_mass_fractions(isotherm; component=1)]

        # dual mode
        dmmodel = fit_dualmode_model(isotherm; uncertainty_method=:JackKnife)
        sheet[GSAHelper.dual_mode_ch_val], sheet[GSAHelper.dual_mode_ch_err] = dmmodel.ch.val, dmmodel.ch.err
        sheet[GSAHelper.dual_mode_b_val], sheet[GSAHelper.dual_mode_b_err] = dmmodel.b.val, dmmodel.b.err
        sheet[GSAHelper.dual_mode_kd_val], sheet[GSAHelper.dual_mode_kd_err] = dmmodel.kd.val, dmmodel.kd.err
        
        dmpredictions = predict_concentration(dmmodel, partial_pressures(isotherm, component=1))
        sheet[GSAHelper.dual_mode_predictions_start, dim=1] = [dm.val for dm in dmpredictions]
        sheet[GSAHelper.dual_mode_predictions_err_start, dim=1] = [dm.err for dm in dmpredictions]

        # write error analysis 
        write_error_analysis(GasSorptionApparatus(), system, xf)

        # if !isnothing(system.transient_system)
        #     write_transient_sorption_system_to_sheet(
        #         system.transient_system, 
        #         xf[TSAHelper.default_sheet_name] 
        #     )
            
        #     # write_kinetic_analysis(xf, isotherm, system.transient_system) # need to #todo implement fugacity in the kinetic analysis
        # else
        #     # println("Not running combined apparatus as no transient data was present")
        # end

    end

    return system

end

function processtemplate(::GasSorptionApparatus, template_path::String; kwargs...)
    return processtemplate(GasSorptionApparatus(), template_path, nothing; kwargs...)
end

function savetemplate(setup::GasSorptionSetup, filepath::String, overwrite=false)
    if !overwrite
        # todo check if file exists, error if it does
    end

    XLSX.openxlsx(filepath, mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, GSAHelper.default_sheet_title)
        
        GSAHelper.add_headers_to_sheet(sheet)

        maybe_missing_val(obj) = ismissing(obj) ? missing : measurement(obj).val
        maybe_missing_err(obj) = ismissing(obj) ? missing : measurement(obj).err

        # add temperature
        sheet[GSAHelper.t] = maybe_missing_val(setup.temperature)
        sheet[GSAHelper.t_err] = maybe_missing_err(setup.temperature)

        # add step pressures
        sheet[GSAHelper.step_start, dim=1] = collect(1:setup.num_steps)
        sheet[GSAHelper.p_ch_in_start, dim=1] = strip_measurement_to_value(setup.charge_chamber_initial_pressures) 
        sheet[GSAHelper.p_ch_fin_start, dim=1] = strip_measurement_to_value(setup.charge_chamber_final_pressures) 
        sheet[GSAHelper.p_samp_start, dim=1] = strip_measurement_to_value(setup.sampling_chamber_final_pressures)
        if !ismissing(setup.charge_chamber_initial_pressures)
            if length(setup.charge_chamber_initial_pressures) > 0
                sheet[GSAHelper.pres_apparatus_err] = maybe_missing_err(setup.charge_chamber_initial_pressures[1])/maybe_missing_val(setup.charge_chamber_initial_pressures[1])
            else
                sheet[GSAHelper.pres_apparatus_err] = 0
            end
        else
            sheet[GSAHelper.pres_apparatus_err] = 0
        end

        # chamber volumes
        sheet[GSAHelper.vc] = maybe_missing_val(setup.charge_chamber_volume)
        sheet[GSAHelper.vc_err] = maybe_missing_err(setup.charge_chamber_volume)
        sheet[GSAHelper.vs] = maybe_missing_val(setup.sampling_chamber_volume)
        sheet[GSAHelper.vs_err] = maybe_missing_err(setup.sampling_chamber_volume)               
        
        # penetrant name
        sheet[GSAHelper.pen_name] = setup.penetrant_name

        # peng robinson parameters
        sheet[GSAHelper.tc] = maybe_missing_val(setup.pen_tc)
        sheet[GSAHelper.tc_err] = maybe_missing_err(setup.pen_tc)
        sheet[GSAHelper.pc] = maybe_missing_val(setup.pen_pc)
        sheet[GSAHelper.pc_err] = maybe_missing_err(setup.pen_pc)
        sheet[GSAHelper.omega] = maybe_missing_val(setup.pen_omega)
        sheet[GSAHelper.omega_err] = maybe_missing_err(setup.pen_omega)
        sheet[GSAHelper.mw] = maybe_missing_val(setup.pen_mw)
        sheet[GSAHelper.mw_err] = maybe_missing_err(setup.pen_mw)
            

        # polymer density
        sheet[GSAHelper.pol_dens] = maybe_missing_val(setup.polymer_density)
        sheet[GSAHelper.pol_dens_err] = maybe_missing_err(setup.polymer_density)
        
        # polymer mass
        sheet[GSAHelper.pol_mass] = maybe_missing_val(setup.polymer_mass)
        sheet[GSAHelper.pol_mass_err] = maybe_missing_err(setup.polymer_mass)

        # polymer name
        sheet[GSAHelper.pol_name] = setup.polymer_name
            
        # num_steps::Int64 # this should be the number of *complete* steps
        # todo not really sure if this is needed yet         

        # beads
        sheet[GSAHelper.nb] = setup.n_beads
        sheet[GSAHelper.nb_vol] = maybe_missing_val(setup.vol_bead)
        sheet[GSAHelper.nb_vol_err] = maybe_missing_err(setup.vol_bead)

        # aluminum mesh
        sheet[GSAHelper.mesh_mass] = maybe_missing_val(setup.mesh_mass)
        sheet[GSAHelper.mesh_mass_err] = maybe_missing_err(setup.mesh_mass)
        sheet[GSAHelper.alum_dens] = maybe_missing_val(setup.alum_dens)
        sheet[GSAHelper.alum_dens_err] = maybe_missing_err(setup.alum_dens)                                 


    end
end