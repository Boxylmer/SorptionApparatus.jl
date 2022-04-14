"""
Holds methods to write out data that can only be found by combining two or more apparatus systems. 

For example, transient and equilibrium data are required to get kinetic and thermodynamic factors of diffusion
"""
module KineticAnalysisHelper
    const default_kinetic_sheet_name = "Sorption Kinetics Analysis"

    const title_col, step_start_col = 1, 2

    const activity_col, activity_err_col = 3, 4
    const conc_col, conc_err_col = 5, 6 
    const conc_g_col, conc_g_err_col = 7, 8 
    const mass_frac_col, mass_frac_err_col = 9, 10  
    const ln_activity_col, ln_activity_err_col = 11, 12  
    const ln_mass_frac_col, ln_mass_frac_err_col = 13, 14  
    const slope_col, slope_err_col = 15, 16  
    const diffusivity_col, diffusivity_err_col = 17, 18   
    const kinetic_factor_col, kinetic_factor_err_col = 19, 20 
    const alpha_col, alpha_err_col = 21, 22
    const step_start_row = 2  
    function _write_kinetic_analysis_header_row(sheet, row_num, title)
        sheet[row_num, title_col] = title
        sheet[row_num, step_start_col] = "Step"
        sheet[row_num, activity_col], sheet[row_num, activity_err_col] = "Activity", "Activity err"
        sheet[row_num, conc_col], sheet[row_num, conc_err_col] = "Conc. (cc/cc)", "Conc. err"
        sheet[row_num, conc_g_col], sheet[row_num, conc_g_err_col] = "Conc. (g)", "Conc. err"
        sheet[row_num, mass_frac_col], sheet[row_num, mass_frac_err_col] = "ω", "ω err"
        sheet[row_num, ln_activity_col], sheet[row_num, ln_activity_err_col] = "ln(Activity)", "ln(a) err"
        sheet[row_num, ln_mass_frac_col], sheet[row_num, ln_mass_frac_err_col] = "ln(ω)", "ln(ω) err"
        sheet[row_num, slope_col], sheet[row_num, slope_err_col] = "d(ln(a))/d(ln(ω))", "d(ln(a))/d(ln(ω)) err"
        sheet[row_num, diffusivity_col], sheet[row_num, diffusivity_err_col] = "D (cm^2/s)", "D err (cm^2/s)"
        sheet[row_num, kinetic_factor_col], sheet[row_num, kinetic_factor_err_col] = "L (cm^2/s)", "L err (cm^2/s)"
        sheet[row_num, alpha_col], sheet[row_num, alpha_err_col] = "α", "α err"    
    end  
    function _write_deconvolved_diffusivity_table(sheet, deconvolution_object, isotherm, diffusivities, row)
        sheet[row, step_start_col, dim=1] = collect(1:length(activities(isotherm; component=1))) 
        sheet[row, activity_col, dim=1] = [a.val for a in activities(isotherm; component=1)]
        sheet[row, activity_err_col, dim=1] = [a.err for a in activities(isotherm; component=1)]
        sheet[row, conc_col, dim=1] = [a.val for a in concentration(isotherm; component=1)]
        sheet[row, conc_err_col, dim=1] = [a.err for a in concentration(isotherm; component=1)]
        sheet[row, conc_g_col, dim=1] = [a.val for a in concentration(isotherm; component=1, pol_units=:g, gas_units=:g)]
        sheet[row, conc_g_err_col, dim=1] = [a.err for a in concentration(isotherm; component=1, pol_units=:g, gas_units=:g)]
        sheet[row, mass_frac_col, dim=1] = [w.val for w in penetrant_mass_fractions(isotherm; component=1)]
        sheet[row, mass_frac_err_col, dim=1] = [w.err for w in penetrant_mass_fractions(isotherm; component=1)]
    
        sheet[row, ln_activity_col, dim=1] = [a.val for a in deconvolution_object.lna]
        sheet[row, ln_activity_err_col, dim=1] = [a.err for a in deconvolution_object.lna]
        sheet[row, ln_mass_frac_col, dim=1] = [w.val for w in deconvolution_object.lnw]
        sheet[row, ln_mass_frac_err_col, dim=1] = [w.err for w in deconvolution_object.lnw]
        sheet[row, slope_col, dim=1] = [item.val for item in deconvolution_object.slopes]
        sheet[row, slope_err_col, dim=1] = [item.err for item in deconvolution_object.slopes]
        
        sheet[row, diffusivity_col, dim=1] = [d.val for d in diffusivities]
        sheet[row, diffusivity_err_col, dim=1] = [d.err for d in diffusivities]
        sheet[row, kinetic_factor_col, dim=1] = [d.val for d in deconvolution_object.kinetic_factors]
        sheet[row, kinetic_factor_err_col, dim=1] = [d.err for d in deconvolution_object.kinetic_factors]
        sheet[row, alpha_col, dim=1] = [d.val for d in deconvolution_object.thermodynamic_factors]
        sheet[row, alpha_err_col, dim=1] = [d.err for d in deconvolution_object.thermodynamic_factors]
    end

end
KiAnCo = KineticAnalysisHelper


"""
    write_kinetic_analysis(excel_file, iso::IsothermData, transient_system::TransientSorptionSystem)
    function write_kinetic_analysis(filepath, iso::IsothermData, transient_system)

Take an isotherm and corresponding `TransientSorptionSystem` (see [`readtemplate`](@ref) for getting the system) and write out the deconvoluted kinetic and thermodynamic and kinetic contributions for all fitted models.

This writer function mainly outputs and formats the results of [`deconvolute_diffusivity`](@ref). 

Note: Generally other functions used to generate `TransientSorptionSystems` and call this method, such as the vapor sorption apparatus. Using this function manually is generally discouraged.
"""
function write_kinetic_analysis(excel_file, iso::IsothermData, transient_system::Union{TransientSorptionSystem, Nothing})
    sheet = XLSX.addsheet!(excel_file)
    XLSX.rename!(sheet, KiAnCo.default_kinetic_sheet_name)
    
    if isnothing(transient_system)
        throw(MissingException("No transient system")) 
    end

    semi_thickness_cm = transient_system.semi_thickness_cm
    num_steps = transient_system.setup.num_steps
    
    # Write out data using the fickian model
    fickian_deconvolution = deconvolute_diffusivity(iso, transient_system.fickian_models, semi_thickness_cm)
    current_row = KiAnCo.step_start_row
    KiAnCo._write_kinetic_analysis_header_row(sheet, current_row - 1, "W/ fickian model")
    KiAnCo._write_deconvolved_diffusivity_table(sheet, fickian_deconvolution, iso, transient_system.fickian_diffusivities, current_row)
    
    bh_deconvolution = deconvolute_diffusivity(iso, transient_system.bh_models, semi_thickness_cm)
    current_row += num_steps + 2
    KiAnCo._write_kinetic_analysis_header_row(sheet, current_row - 1, "W/ BH model")
    KiAnCo._write_deconvolved_diffusivity_table(sheet, bh_deconvolution, iso, transient_system.bh_diffusivities, current_row)
    
    mbh_deconvolution = deconvolute_diffusivity(iso, transient_system.mbh_models, semi_thickness_cm)
    current_row += num_steps + 2
    KiAnCo._write_kinetic_analysis_header_row(sheet, current_row - 1, "W/ MBH model")
    KiAnCo._write_deconvolved_diffusivity_table(sheet, mbh_deconvolution, iso, transient_system.mbh_diffusivities, current_row)
end

function write_kinetic_analysis(filepath::AbstractString, iso::IsothermData, transient_system::Union{TransientSorptionSystem, Nothing})  # todo untested
    XLSX.openxlsx(filepath, mode="rw") do xf  
        return write_kinetic_analysis(xf, iso, transient_system)
    end
end


module IsostericHeatAnalysisHelper
    const default_isosteric_heat_sheet_name = "Isosteric Heat Analysis"
end
IHAH = IsostericHeatAnalysisHelper

"""
    write_isosteric_heat_analysis(workbook, analysis::IsostericHeatAnalysis; name)
    write_isosteric_heat_analysis(filepath::AbstractString, analysis::IsostericHeatAnalysis; name)

Write an isosteric heat analysis to an Excel workbook. If a worksheet name isn't provided, a default one is used.
See [`isosteric_heat_of_sorption`](@ref)
"""
function write_isosteric_heat_analysis(workbook::XLSX.XLSXFile, analysis::IsostericHeatAnalysis; name=IHAH.default_isosteric_heat_sheet_name)
    sheet = XLSX.addsheet!(workbook, name)
    
    num_isotherms = length(analysis.isotherms)
    writer_row_position = 1
    writer_col_position = 1
    
    # first show the fitted, interpolated data
    sheet[writer_row_position, writer_col_position] = "Fitted, interpolated data" 
    writer_row_position += 1
    sheet[writer_row_position, writer_col_position] = "Sampled CC/CC"
    sheet[writer_row_position, writer_col_position + 1, dim=2] = [conc.val for conc in analysis.sampled_concentrations]
    writer_row_position +=1
    sheet[writer_row_position, writer_col_position] = "T"
    sheet[writer_row_position + 1, writer_col_position, dim=1] = [temp.val for temp in analysis.temperatures]
    for (idx, pressure_vector) in enumerate(analysis.pressure_vectors)
        sheet[writer_row_position, writer_col_position + idx] = "P[MPa]"
        sheet[writer_row_position + 1, writer_col_position + idx, dim=1] = [meas.val for meas in pressure_vector]
    end

    # next show the resulting ln(P) vs 1/T curves 
    writer_row_position += num_isotherms + 2

    sheet[writer_row_position, writer_col_position] = "Manipulated data for getting ∂ln(P)/∂(1/T)" 
    writer_row_position += 1
    sheet[writer_row_position, writer_col_position] = "Sampled CC/CC"
    sheet[writer_row_position, writer_col_position + 1, dim=2] = [conc.val for conc in analysis.sampled_concentrations]
    writer_row_position +=1
    sheet[writer_row_position, writer_col_position] = "1/T (K^-1)"
    sheet[writer_row_position + 1, writer_col_position, dim=1] = [temp.val for temp in analysis.inverse_temperature_vector]
    for (idx, vec) in enumerate(analysis.ln_pressure_vectors)
        sheet[writer_row_position, writer_col_position + idx] = "ln(P[MPa])"
        sheet[writer_row_position + 1, writer_col_position + idx, dim=1] = [vec_meas.val for vec_meas in vec] 
    end

    # finally, show the found isosteric heats as a function of concentration
    writer_row_position += num_isotherms + 2
    sheet[writer_row_position, writer_col_position] = "Isosteric heats as a function of concentration" 
    writer_row_position += 1
    sheet[writer_row_position, writer_col_position] = "Conc (CC/CC)" 
    sheet[writer_row_position, writer_col_position + 1] = "ΔH_s (J/mol)" 
    sheet[writer_row_position, writer_col_position + 2] = "ΔH_s err (J/mol)" 
    writer_row_position += 1
    sheet[writer_row_position, writer_col_position, dim=1] = [conc.val for conc in analysis.sampled_concentrations]
    sheet[writer_row_position, writer_col_position + 1, dim=1] = [heat.val for heat in analysis.isosteric_heat_at_conc]
    sheet[writer_row_position, writer_col_position + 2, dim=1] = [heat.err for heat in analysis.isosteric_heat_at_conc]
 
end

function write_isosteric_heat_analysis(filepath::AbstractString, analysis::IsostericHeatAnalysis; name=IHAH.default_isosteric_heat_sheet_name)
    XLSX.openxlsx(filepath, mode="rw") do xf  
        return write_isosteric_heat_analysis(xf, analysis; name=name)
    end
end

function write_vant_hoff_dual_mode_analysis(workbook::XLSX.XLSXFile, analysis::VantHoffDualModeAnalysis; name="Van't Hoff Dual Mode Analysis")
    sheet = XLSX.addsheet!(workbook, name)
    wrp = 1  # writer row position
    wcp = 1  # writer col position

    function write_vector_of_maybe_measurements(row, col, vec; write_errors=true)
        sheet[row, col, dim=1] = strip_measurement_to_value(vec)
        if write_errors
            if eltype(vec) <: Measurement
                sheet[row, col+1, dim=1] = [vecval.err for vecval in vec]
            else
                sheet[row, col+1, dim=1] = ["No uncertainty" for _ in vec]
            end
        end
    end
    function write_value_of_maybe_measurements(row, col, val; write_errors=true)
        sheet[row, col] = strip_measurement_to_value(val)
        if write_errors
            if typeof(val) <: Measurement
                sheet[row, col+1] = val.err
            else
                sheet[row, col+1] = "No uncertainty"
            end
        end
    end
    function write_parameter_headers(row_start, col_start; row_col_iter=(1, 0)) 
        rowpos = 0
        colpos = 0
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "CH'_0 (CC/CC)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "b_0 (1/MPa)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "kd_0 ((CC/CC)/MPa)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "m_CH` (1/K)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "ΔHb (J/mol)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "ΔHkd (J/mol)"
    end
    function write_matrix_at_position(row, col, mat)
        for idx in CartesianIndices(mat)
            sheet[row + idx[1] - 1, col + idx[2] - 1] = mat[idx[1], idx[2]]
        end
    end

    sheet[wrp, wcp] = "Parameters"
    sheet[wrp, wcp + 1] = "Initial"
    sheet[wrp, wcp + 2] = "Final"
    sheet[wrp, wcp + 3] = "Final err"

    write_parameter_headers(wrp, wcp)
    
    write_value_of_maybe_measurements(wrp+1, wcp+1, analysis.ch0_initial; write_errors=false); write_value_of_maybe_measurements(wrp+1, wcp+2, analysis.ch0_final)
    write_value_of_maybe_measurements(wrp+2, wcp+1, analysis.b0_initial; write_errors=false); write_value_of_maybe_measurements(wrp+2, wcp+2, analysis.b0_final)
    write_value_of_maybe_measurements(wrp+3, wcp+1, analysis.kd0_initial; write_errors=false); write_value_of_maybe_measurements(wrp+3, wcp+2, analysis.kd0_final)
    write_value_of_maybe_measurements(wrp+4, wcp+1, analysis.mch_initial; write_errors=false); write_value_of_maybe_measurements(wrp+4, wcp+2, analysis.mch_final)
    write_value_of_maybe_measurements(wrp+5, wcp+1, analysis.ΔHb_initial; write_errors=false); write_value_of_maybe_measurements(wrp+5, wcp+2, analysis.ΔHb_final)
    write_value_of_maybe_measurements(wrp+6, wcp+1, analysis.ΔHkd_initial; write_errors=false); write_value_of_maybe_measurements(wrp+6, wcp+2, analysis.ΔHkd_final)
    
    sheet[wrp+7, wcp] = "RSS"
    write_value_of_maybe_measurements(wrp+7, wcp+1, analysis.initial_rss; write_errors=false); write_value_of_maybe_measurements(wrp+7, wcp+2, analysis.final_rss)

    wcp+=5
    sheet[wrp, wcp] = "Covariance Matrix"
    write_parameter_headers(wrp, wcp)
    write_parameter_headers(wrp, wcp; row_col_iter=(0, 1))
    write_matrix_at_position(wrp + 1, wcp + 1, analysis.covariance_matrix)

    wcp+=9
    sheet[wrp, wcp] = "1/T"
    write_vector_of_maybe_measurements(wrp + 1, wcp, 1 ./ [temperature(iso) for iso in analysis.given_isotherms]; write_errors=false)
    sheet[wrp, wcp + 1] = "ch (standalone)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 1, [model.ch for model in analysis.standalone_models]; write_errors=false)
    sheet[wrp, wcp + 2] = "ln(b) (standalone)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 2, log.([model.b for model in analysis.standalone_models]); write_errors=false)
    sheet[wrp, wcp + 3] = "ln(kd) (standalone)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 3, log.([model.kd for model in analysis.standalone_models]); write_errors=false)
    sheet[wrp, wcp + 4] = "ch (initial)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 4, [model.ch for model in analysis.initial_models]; write_errors=false)
    sheet[wrp, wcp + 5] = "ln(b) (initial)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 5, log.([model.b for model in analysis.initial_models]); write_errors=false)
    sheet[wrp, wcp + 6] = "ln(kd) (initial)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 6, log.([model.kd for model in analysis.initial_models]); write_errors=false)
    sheet[wrp, wcp + 7] = "ch (final)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 7, [model.ch for model in analysis.final_models]; write_errors=false)
    sheet[wrp, wcp + 8] = "ln(b) (final)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 8, log.([model.b for model in analysis.final_models]); write_errors=false)
    sheet[wrp, wcp + 9] = "ln(kd) (final)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 9, log.([model.kd for model in analysis.final_models]); write_errors=false)
    wcp=1

    wrp+=7

    wrp+= max(2, length(analysis.given_isotherms) - 8)  # move the writer by the longest vector
    sheet[wrp, wcp] = "Isotherms"
    wrp+=1


    for idx in eachindex(analysis.final_models, analysis.initial_models, analysis.final_predicted_isotherms, analysis.initial_predicted_isotherms)
        sheet[wrp, wcp] = "Temp (K):"
        write_value_of_maybe_measurements(wrp, wcp+1, temperature(analysis.given_isotherms[idx]))
        
        # isotherm stuff
        wrp+=1
        if analysis.use_fugacity sheet[wrp, wcp] = "Fugacity (MPa)"
        else sheet[wrp, wcp] = "Pressure (MPa)" end
        sheet[wrp, wcp + 1], sheet[wrp, wcp + 2], sheet[wrp, wcp + 3], sheet[wrp, wcp + 4] = "Err (MPa)", "CC/CC", "CC/CC err", "DM standalone (CC/CC)"
        sheet[wrp, wcp + 5], sheet[wrp, wcp + 6], sheet[wrp, wcp + 7], sheet[wrp, wcp + 8], sheet[wrp, wcp + 9]  = "DM init (CC/CC)", "DM final (CC/CC)", "err (CC/CC)", "DM Henry Contrib (CC/CC)", "err (CC/CC)"
        sheet[wrp, wcp + 10], sheet[wrp, wcp + 11] = "DM Langmuir Contrib (CC/CC)", "err (CC/CC)"
        
        
        wrp+=1
        
        if analysis.use_fugacity
            used_pressures = fugacities(analysis.given_isotherms[idx]; component=1)
        else
            used_pressures = partial_pressures(analysis.given_isotherms[idx]; component=1)
        end
        write_vector_of_maybe_measurements(wrp, wcp, fugacities(analysis.given_isotherms[idx]; component=1))
        write_vector_of_maybe_measurements(wrp, wcp + 2, concentration(analysis.given_isotherms[idx]; component=1))
        write_vector_of_maybe_measurements(wrp, wcp + 4, strip_measurement_to_value(predict_concentration(analysis.standalone_models[idx], used_pressures)))
        write_vector_of_maybe_measurements(wrp, wcp + 5, strip_measurement_to_value(predict_concentration(analysis.initial_models[idx], used_pressures)))
        write_vector_of_maybe_measurements(wrp, wcp + 6, predict_concentration(analysis.final_models[idx], used_pressures))
        write_vector_of_maybe_measurements(wrp, wcp + 8, henry_mode_concentration(analysis.final_models[idx], used_pressures))
        write_vector_of_maybe_measurements(wrp, wcp + 10, langmuir_mode_concentration(analysis.final_models[idx], used_pressures))
        
        # dual mode fittings
        wrp+=0
        wcp+=14
        sheet[wrp - 1, wcp] = "Dual Mode Parameters"
        sheet[wrp + 1, wcp], sheet[wrp + 2, wcp], sheet[wrp + 3, wcp] = "CH' (CC/CC)", "b (1/MPa)", "kd ((CC/CC)/MPa)"
        sheet[wrp, wcp + 1], sheet[wrp, wcp + 2], sheet[wrp, wcp + 3], sheet[wrp, wcp + 4], sheet[wrp, wcp + 5] = "standalone", "standalone err", "initial", "final", "final err"
        wrp+=1
        write_vector_of_maybe_measurements(wrp, wcp + 1, [analysis.standalone_models[idx].ch, analysis.standalone_models[idx].b, analysis.standalone_models[idx].kd])
        write_vector_of_maybe_measurements(wrp, wcp + 3, [analysis.initial_models[idx].ch, analysis.initial_models[idx].b, analysis.initial_models[idx].kd]; write_errors=false)
        write_vector_of_maybe_measurements(wrp, wcp + 4, [analysis.final_models[idx].ch, analysis.final_models[idx].b, analysis.final_models[idx].kd])

        wrp+= 1 + num_steps(analysis.given_isotherms[idx])
        wcp = 1

    end

end

function write_vant_hoff_dual_mode_analysis(filepath::AbstractString, analysis::VantHoffDualModeAnalysis; name="Van't Hoff Dual Mode Analysis")
    XLSX.openxlsx(filepath, mode="rw") do xf  
        return write_vant_hoff_dual_mode_analysis(xf, analysis; name=name)
    end
end


"""
    write_dual_mode_diffusivity_deconvolution(filepath, analysis; [analysis_name])
    write_dual_mode_diffusivity_deconvolution(filepath, analyses; [analysis_names])

Write a Dual Mode diffusivity deconvolution `analysis` out to a .csv titles by `filepath`. If analysis names are not provided, default names are created.
Methods are provided for a single analysis, or vectors of analyses (all written to one file).

See [`dual_mode_diffusivity_deconvolution`](@ref)
"""
function write_dual_mode_diffusivity_deconvolution(
    filepath::AbstractString, analysis::DualModeDiffusivityDeconvolution; 
    analysis_name = missing)
    write_dual_mode_diffusivity_deconvolution(filepath, [analysis]; analyses_names = [analysis_name])
end
function write_dual_mode_diffusivity_deconvolution(
    filepath::AbstractString,
    analyses::AbstractVector{DualModeDiffusivityDeconvolution};
    analyses_names::AbstractVector = missing)
    permeability_deconvolution_results = []

    if ismissing(analyses_names)
        analyses_names = "Analysis " .* string.(1:length(analyses))
    end

    for (idx, analysis) in enumerate(analyses)
        langmuir_diffusivity, henry_diffusivity = analysis.langmuir_mode_diffusivity, analysis.henry_mode_diffusivity
        f = langmuir_diffusivity / henry_diffusivity
        x_values = [xval.val for xval in analysis.regression_x_values]

        # ch_b_over_1_plus_fp = model.ch * model.b / (1 .+ model.b .* corresponding_fugacities)
        push!(permeability_deconvolution_results, [analyses_names[idx]])
        push!(permeability_deconvolution_results, ["Temperature (K): ", strip_measurement_to_value(analysis.temperature_k)])
        if analysis.model.use_fugacity == false
            push!(permeability_deconvolution_results, ["Pressures (MPa): ", analysis.pressures_mpa...])
        else
            push!(permeability_deconvolution_results, ["Fugacities (MPa): ", analysis.pressures_mpa...])
        end
        push!(permeability_deconvolution_results, ["Permeabilities (Barrer): ", analysis.permeabilities_barrer...])
        push!(permeability_deconvolution_results, ["Langmuir Diffusivity (cm^2/s): ", langmuir_diffusivity.val, "+/-", langmuir_diffusivity.err])
        push!(permeability_deconvolution_results, ["Henry Diffusivity (cm^2/s): ", henry_diffusivity.val, "+/-", henry_diffusivity.err])
        push!(permeability_deconvolution_results, ["F (Dh/Dd): ", f.val, "+/-", f.err])
        push!(permeability_deconvolution_results, ["ch'*b/(1+bf)", x_values...])
    end
    writedlm(
        joinpath(
            @__DIR__, 
            filepath, 
        ), permeability_deconvolution_results, ",")

end