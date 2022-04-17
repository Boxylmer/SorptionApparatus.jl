"""
Holds methods to write out data that can only be found by combining two or more apparatus systems. 

For example, transient and equilibrium data are required to get kinetic and thermodynamic factors of diffusion
"""
module KineticAnalysisHelper
    using MembraneBase
    
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
    fickian_deconvolution = MobilityFactorAnalysis(iso, transient_system.fickian_models, semi_thickness_cm)
    current_row = KiAnCo.step_start_row
    KiAnCo._write_kinetic_analysis_header_row(sheet, current_row - 1, "W/ fickian model")
    KiAnCo._write_deconvolved_diffusivity_table(sheet, fickian_deconvolution, iso, transient_system.fickian_diffusivities, current_row)
    
    bh_deconvolution = MobilityFactorAnalysis(iso, transient_system.bh_models, semi_thickness_cm)
    current_row += num_steps + 2
    KiAnCo._write_kinetic_analysis_header_row(sheet, current_row - 1, "W/ BH model")
    KiAnCo._write_deconvolved_diffusivity_table(sheet, bh_deconvolution, iso, transient_system.bh_diffusivities, current_row)
    
    mbh_deconvolution = MobilityFactorAnalysis(iso, transient_system.mbh_models, semi_thickness_cm)
    current_row += num_steps + 2
    KiAnCo._write_kinetic_analysis_header_row(sheet, current_row - 1, "W/ MBH model")
    KiAnCo._write_deconvolved_diffusivity_table(sheet, mbh_deconvolution, iso, transient_system.mbh_diffusivities, current_row)
end

function write_kinetic_analysis(filepath::AbstractString, iso::IsothermData, transient_system::Union{TransientSorptionSystem, Nothing})  # todo untested
    XLSX.openxlsx(filepath, mode="rw") do xf  
        return write_kinetic_analysis(xf, iso, transient_system)
    end
end

