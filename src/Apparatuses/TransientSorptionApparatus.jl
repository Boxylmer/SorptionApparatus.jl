struct TransientSorptionApparatus end

struct TransientSorptionSetup{STT, NST, SMT, ICT, RSCT}
    step_data::AbstractVector{<:TransientStepData}
    semi_thickness::STT
    num_steps::NST
    interpolation_method::SMT
    interpolation_count::ICT
    resampling_count::RSCT
end

struct TransientSorptionSystem{TSS, FSPT, BHSPT, MBHSPT, FFET, BFET, MFET, STCT}
    setup::TSS  # need to fix with actual VaporSorptionSystem type
    fickian_models::Vector{FickianSorptionModel}
    bh_models::Vector{BerensHopfenbergSorptionModel}
    mbh_models::Vector{ModifiedBerensHopfenbergSorptionModel}
    resampled_steps::Vector{TransientStepData}
    fickian_step_predictions::FSPT
    bh_step_predictions::BHSPT
    mbh_step_predictions::MBHSPT
    fickian_fit_errors::FFET
    bh_fit_errors::BFET
    mbh_fit_errors::MFET

    fickian_diffusivities::Vector{Any}
    bh_diffusivities::Vector{Any}
    mbh_diffusivities::Vector{Any}

    semi_thickness_cm::STCT
end 

function TransientSorptionSystem(template_filepath::String; kwargs...)
    tss = TransientSorptionSetup(template_filepath)
    return TransientSorptionSystem(tss; kwargs...)
end

TransientSorptionSetup(path::String; kwargs...) = readtemplate(TransientSorptionApparatus(), path; kwargs...)

function TransientSorptionSystem(tss::TransientSorptionSetup; verbose=false)
    fickian_models = Vector{FickianSorptionModel}()
    bh_models = Vector{BerensHopfenbergSorptionModel}()
    mbh_models = Vector{ModifiedBerensHopfenbergSorptionModel}()
    resampled_steps = Vector{TransientStepData}()

    fickian_step_predictions = Vector{Vector{Number}}()
    bh_step_predictions = Vector{Vector{Number}}()
    mbh_step_predictions = Vector{Vector{Number}}()

    fickian_fit_errors = Vector{Number}()
    bh_fit_errors = Vector{Number}()
    mbh_fit_errors = Vector{Number}()

    fickian_diffusivities = Vector{Any}()
    bh_diffusivities = Vector{Any}()
    mbh_diffusivities = Vector{Any}()
    semi_thickness_cm = tss.semi_thickness
    if verbose
        n = length(tss.step_data) * 4
        p = Progress(n, dt=0.5,
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=10)
    end
    for step_data in tss.step_data
        push!(
            fickian_models, 
            fit_transient_sorption_model(
                step_data, :FickianSorptionModel; interpolation_method=tss.interpolation_method, interpolation_datapoints=tss.interpolation_count,
                uncertainty_method=:Bootstrap)
        )
        if verbose; next!(p); end
        push!(
            bh_models,
            fit_transient_sorption_model(
                step_data, :BerensHopfenbergSorptionModel; interpolation_method=tss.interpolation_method, interpolation_datapoints=tss.interpolation_count,
                uncertainty_method=:Bootstrap, num_uncertainty_resamples=tss.resampling_count)
        )
        if verbose; next!(p); end
        push!(
            mbh_models,
            fit_transient_sorption_model(
                step_data, :ModifiedBerensHopfenbergSorptionModel; interpolation_method=tss.interpolation_method, interpolation_datapoints=tss.interpolation_count,
                uncertainty_method=:Bootstrap)
        )
        if verbose; next!(p); end
        push!(
            resampled_steps,
            resample(step_data, tss.interpolation_count, tss.interpolation_method)
        )
        push!(
            fickian_step_predictions,
            predict_sorption(fickian_models[end], step_data.time)
        )
        push!(
            bh_step_predictions,
            predict_sorption(bh_models[end], step_data.time)
        )
        push!(
            mbh_step_predictions,
            predict_sorption(mbh_models[end], step_data.time)
        )
        push!(
            fickian_fit_errors,
            rss(fickian_models[end], step_data)
        )
        push!(
            bh_fit_errors,
            rss(bh_models[end], step_data)
        )
        push!(
            mbh_fit_errors,
            rss(mbh_models[end], step_data)
        )
        push!(
            fickian_diffusivities, 
            get_diffusivity(fickian_models[end], semi_thickness_cm)
        )
        push!(
            bh_diffusivities, 
            get_diffusivity(bh_models[end], semi_thickness_cm)
        )
        push!(
            mbh_diffusivities, 
            get_diffusivity(mbh_models[end], semi_thickness_cm)
        )
        if verbose; next!(p); end
    end

    semi_thickess = tss.semi_thickness
    # construct the system struct
    system = TransientSorptionSystem(
        tss, fickian_models, bh_models, mbh_models, resampled_steps, 
        fickian_step_predictions, bh_step_predictions, mbh_step_predictions,
        fickian_fit_errors, bh_fit_errors, mbh_fit_errors,
        fickian_diffusivities, bh_diffusivities, mbh_diffusivities, 
        semi_thickess
    )
    return system

end

module TSAHelper
    ###### define templating system ######
    const default_sheet_name = "Transient Sorption Input"
    const default_file_name = "Transient Sorption Template (rename this!).xlsx"

    const overall_header, overall_val, overall_err = "A1", "B1", "C1"
    const semi_thickness_header, semi_thickness_value, semi_thickness_error = "A3", "B3", "C3"
    const num_steps_header, num_steps_value, num_steps_error = "A4", "B4", "C4"
    const interpolation_explanation = "B5"
    const interpolation_method_header, interpolation_method_value, interpolation_method_err = "A6", "B6", "C6"
    const interpolation_count_header, interpolation_count_value, interpolation_count_err = "A7", "B7", "C7"
    const resampling_explanation = "B8"
    const resampling_count_header, resampling_count_value, resampling_count_err = "A9", "B9", "C9"

    const transient_pressure_mode_header, transient_pressure_mode_value, transient_pressure_mode_explanation = "A11","B11", "C11"

    # sorption steps (data positions relative to starting points)
    const step_start_col, step_start_row = 5, 1
    const step_col_spacing = 10  # rows between step data tables

    # add these to the starting point to get their positions
    const relative_mf_row_pos, relative_kf_row_pos, relative_mr_row_pos, relative_kr_row_pos, relative_beta_row_pos = 1, 2, 3, 4, 5
    const relative_d_row_pos, relative_fit_error_row_pos = 6, 7
    const relative_fick_col_pos,     relative_bh_col_pos,      relative_mbh_col_pos = 2, 4, 6
    const relative_fick_err_col_pos, relative_bh_err_col_pos,  relative_mbh_err_col_pos =  3, 5, 7

    const relative_data_row_start = 8  # assumed that we start at relative col 0
    const relative_input_time_col = 0
    const relative_input_dimless_sorption_col = 1
    const relative_output_resampled_time_col = 8
    const relative_output_resampled_sorption_col = 9

    const default_num_steps = 9
end
function generatetemplate(::TransientSorptionApparatus; folder = "", filename=TSAHelper.default_file_name, standalone=true) # todo check if the file exists instead of standalone
    writemode="w"
    if !standalone writemode = "rw" end
    XLSX.openxlsx(joinpath(folder, filename), mode=writemode) do xf
        if !standalone
            sheet = XLSX.addsheet!(xf, TSAHelper.default_sheet_name)
        else
            sheet = xf[1]
        end
        
        XLSX.rename!(sheet, TSAHelper.default_sheet_name)
        sheet[TSAHelper.overall_header] = "Overall Properties"; sheet[TSAHelper.overall_val] = "Value"; sheet[TSAHelper.overall_err] = "Uncertainty"
        
        sheet[TSAHelper.semi_thickness_header] = "semi-thickness (cm)"; sheet[TSAHelper.semi_thickness_value] = ""; sheet[TSAHelper.semi_thickness_error] ="";
        sheet[TSAHelper.num_steps_header], sheet[TSAHelper.num_steps_value], sheet[TSAHelper.num_steps_error] = "# steps", TSAHelper.default_num_steps, "---"
        
        sheet[TSAHelper.interpolation_method_header], sheet[TSAHelper.interpolation_method_value], sheet[TSAHelper.interpolation_method_err] = "Interpolation method", "root", ">Options<\n
        \"linear\", \"root\", \"log\", \"none\""
        sheet[TSAHelper.interpolation_count_header], sheet[TSAHelper.interpolation_count_value], sheet[TSAHelper.interpolation_count_err] = "Interpolation size", 1000, "---"
        
        sheet[TSAHelper.resampling_count_header], sheet[TSAHelper.resampling_count_value], sheet[TSAHelper.resampling_count_err] = "Resampling count", 50, "---"

        sheet[TSAHelper.transient_pressure_mode_header], sheet[TSAHelper.transient_pressure_mode_value] = "transient pressure mode", "FALSE", "<- TRUE if using transient pressure instead of dimensionless sorption."

        # nice giblets
        sheet[TSAHelper.step_start_row, TSAHelper.step_start_col + TSAHelper.default_num_steps * TSAHelper.step_col_spacing] = "-> Copy the step template pattern and adjust the \"# steps\" value accordingly to add steps. ->"
        sheet[TSAHelper.interpolation_explanation] = ">Interpolation<\n
        Interpolate the input data to some size of points according to some interpolation function\n
        making the data easier to fit while losing little information in the process.\n 
        It can also help the parameter optimization understand where the \"valuable\" information is.\n
        You almost always want 'root', as ideal fickian sorption is mostly linear with it.\n 
        'linear' is even interpolation (it will just reduce the data size)\n 
        and 'log' will very heavily cut down long sorption steps (think on the order of months or years, not days)."

        sheet[TSAHelper.TSAHelper.resampling_explanation] = "Resampling\n
        This is the number of times the data is randomly resampled and fit.\n This is done to determine a histogram "
        # add step tables
        for step_idx in 1:TSAHelper.default_num_steps
            position_idx = step_idx - 1  # off by 1 errors are the worst
            start_col = TSAHelper.step_start_col + position_idx * TSAHelper.step_col_spacing
            start_row = TSAHelper.step_start_row # right now we don't move the row

            # label corner
            sheet[start_row, start_col] = "Step " * string(step_idx)

            # row headers
            sheet[start_row + TSAHelper.relative_mf_row_pos, start_col] = "M (fick)"
            sheet[start_row + TSAHelper.relative_kf_row_pos, start_col] = "kf (1/s)"
            sheet[start_row + TSAHelper.relative_mr_row_pos, start_col] = "M (relax)"
            sheet[start_row + TSAHelper.relative_kr_row_pos, start_col] = "kr (1/s)"
            sheet[start_row + TSAHelper.relative_beta_row_pos, start_col] = "beta"
            sheet[start_row + TSAHelper.relative_d_row_pos, start_col] = "D (cm^2/s)"
            sheet[start_row + TSAHelper.relative_fit_error_row_pos, start_col] = "Fit error"

            # column headers
            sheet[start_row, start_col + TSAHelper.relative_fick_col_pos] = "Fick. Model"
            sheet[start_row, start_col + TSAHelper.relative_fick_err_col_pos] = "Fick. err"
            sheet[start_row, start_col + TSAHelper.relative_bh_col_pos] = "BH Model"
            sheet[start_row, start_col + TSAHelper.relative_bh_err_col_pos] = "BH err"
            sheet[start_row, start_col + TSAHelper.relative_mbh_col_pos] = "MBH Model"
            sheet[start_row, start_col + TSAHelper.relative_mbh_err_col_pos] = "MBH err"

            # data placeholders
            sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_fick_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_bh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_fick_col_pos] = "TBD"
            
            sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_bh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = "TBD"
            
            sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_fick_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_bh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = "TBD"
        
            sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_fick_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_bh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = "TBD"

            sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_fick_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_bh_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = "---"
            sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = "TBD"

            sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_fick_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_bh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = "TBD"
            sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = "TBD"

            # data headers and placement info 
            sheet[start_row + TSAHelper.relative_data_row_start, start_col + TSAHelper.relative_input_time_col] = "t - t0"
            sheet[start_row + TSAHelper.relative_data_row_start, start_col + TSAHelper.relative_input_dimless_sorption_col] = "Dim.less sorption or transient pressure (Pa)"
            sheet[start_row + TSAHelper.relative_data_row_start, start_col + TSAHelper.relative_output_resampled_time_col] = "Resampled time"
            sheet[start_row + TSAHelper.relative_data_row_start, start_col + TSAHelper.relative_output_resampled_sorption_col] = "Resampled dim.less sorption"
            
            # nice (repeating) giblets
            sheet[start_row + TSAHelper.relative_data_row_start + 1, start_col + TSAHelper.relative_input_time_col] = "Paste your (from sorption start) time data here!"
            sheet[start_row + TSAHelper.relative_data_row_start + 2, start_col + TSAHelper.relative_input_dimless_sorption_col] = "Paste your (dimensionless) sorption data here!"
            
        end
        # add future results headers
        
    end
    return nothing
end

function readtemplate(::TransientSorptionApparatus, path::String; sheet_name=TSAHelper.default_sheet_name, apparatus_setup=nothing)  # read a template into a sorption setup struct
    xf = XLSX.readxlsx(path)
    sheet = xf[sheet_name]
    
    num_steps = sheet[TSAHelper.num_steps_value] 
    semi_thickness = sheet[TSAHelper.semi_thickness_value] ± (ismissing(sheet[TSAHelper.semi_thickness_error]) ? 0 : sheet[TSAHelper.semi_thickness_error])
    _interpolation_method = sheet[TSAHelper.interpolation_method_value]
    resampling_count = sheet[TSAHelper.resampling_count_value]

    if _interpolation_method == "none"
        interpolation_method = nothing
    elseif _interpolation_method == "linear"
        interpolation_method = :Linear
    elseif _interpolation_method == "root"
        interpolation_method = :Root
    elseif _interpolation_method == "log"
        interpolation_method = :Log
    else
        interpolation_method = nothing
    end

    if !isnothing(interpolation_method)
        interpolation_count = sheet[TSAHelper.interpolation_count_value]
    else
        interpolation_count = nothing
    end
   
    transient_step_vector = Vector{TransientStepData}()
    is_in_transient_pressure_mode = safe_parse_bool(sheet[TSAHelper.transient_pressure_mode_value])
    for step_idx in 1:num_steps
        position_idx = step_idx - 1  # off by 1 errors are the worst
        start_col = TSAHelper.step_start_col + position_idx * TSAHelper.step_col_spacing
        start_row = TSAHelper.step_start_row # right now we don't move the row

        time_colidx = start_col + TSAHelper.relative_input_time_col
        sorption_colidx = start_col + TSAHelper.relative_input_dimless_sorption_col
        
        time_vector = convert(Array{Float64}, XLSX.gettable(
            sheet, XLSX.ColumnRange(time_colidx, time_colidx); first_row=start_row + TSAHelper.relative_data_row_start, infer_eltypes=false
            )[1][1])

        dimensionless_sorption_or_pressure_vector = convert(Array{Float64},XLSX.gettable(
            sheet, XLSX.ColumnRange(sorption_colidx, sorption_colidx); first_row=start_row + TSAHelper.relative_data_row_start, infer_eltypes=false
            )[1][1])
            
        if is_in_transient_pressure_mode
            pressure_vector = dimensionless_sorption_or_pressure_vector
            equilibrium_moles_sorbed = moles_sorbed_during_step(apparatus_setup, step_idx)

            sorption_vector = [dimensionless_mass_sorbed_during_step(apparatus_setup, step_idx, equilibrium_moles_sorbed, pres) for pres in pressure_vector]
        else
            sorption_vector = dimensionless_sorption_or_pressure_vector
        end
        transient_data = TransientStepData(time_vector, sorption_vector)
        push!(transient_step_vector, transient_data)
    end
    setup = TransientSorptionSetup(transient_step_vector, semi_thickness, num_steps, interpolation_method, interpolation_count, resampling_count)
    return setup
end

function write_transient_sorption_system_to_sheet(system::TransientSorptionSystem, sheet)
    # set up the for loop for each step
    for step_idx in 1:system.setup.num_steps
        position_idx = step_idx - 1  # off by 1 errors are the worst
        start_col = TSAHelper.step_start_col + position_idx * TSAHelper.step_col_spacing
        start_row = TSAHelper.step_start_row # right now we don't move the row
    
        # start placing the fittings
        sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_fick_col_pos] = system.fickian_models[step_idx].m_f.val
        sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = system.fickian_models[step_idx].m_f.err
        sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_bh_col_pos] = system.bh_models[step_idx].m_f.val
        sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = system.bh_models[step_idx].m_f.err
        sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = system.mbh_models[step_idx].m_f.val
        sheet[start_row + TSAHelper.relative_mf_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = system.mbh_models[step_idx].m_f.err
        
        sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_fick_col_pos] = system.fickian_models[step_idx].k_f.val
        sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = system.fickian_models[step_idx].k_f.err
        sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_bh_col_pos] = system.bh_models[step_idx].k_f.val
        sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = system.bh_models[step_idx].k_f.err
        sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = system.mbh_models[step_idx].k_f.val
        sheet[start_row + TSAHelper.relative_kf_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = system.mbh_models[step_idx].k_f.err
        
        sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_bh_col_pos] = system.bh_models[step_idx].m_r.val
        sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = system.bh_models[step_idx].m_r.err
        sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = system.mbh_models[step_idx].m_r.val
        sheet[start_row + TSAHelper.relative_mr_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = system.mbh_models[step_idx].m_r.err
    
        sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_bh_col_pos] = system.bh_models[step_idx].k_r.val
        sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = system.bh_models[step_idx].k_r.err
        sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = system.mbh_models[step_idx].k_r.val
        sheet[start_row + TSAHelper.relative_kr_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = system.mbh_models[step_idx].k_r.err

        sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = system.mbh_models[step_idx].beta.val
        sheet[start_row + TSAHelper.relative_beta_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = system.mbh_models[step_idx].beta.err
  
        sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_fick_col_pos] = ismissing(system.fickian_diffusivities[step_idx]) ? missing : system.fickian_diffusivities[step_idx].val

        sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = ismissing(system.fickian_diffusivities[step_idx]) ? missing : system.fickian_diffusivities[step_idx].err
        sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_bh_col_pos] = ismissing(system.bh_diffusivities[step_idx]) ? missing :  system.bh_diffusivities[step_idx].val
        sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = ismissing(system.bh_diffusivities[step_idx]) ? missing : system.bh_diffusivities[step_idx].err
        sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = ismissing(system.mbh_diffusivities[step_idx]) ? missing : system.mbh_diffusivities[step_idx].val
        sheet[start_row + TSAHelper.relative_d_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = ismissing(system.mbh_diffusivities[step_idx]) ? missing : system.mbh_diffusivities[step_idx].err
        
        sheet[start_row + TSAHelper.relative_fit_error_row_pos, start_col + TSAHelper.relative_fick_col_pos] = system.fickian_fit_errors[step_idx].val
        sheet[start_row + TSAHelper.relative_fit_error_row_pos, start_col + TSAHelper.relative_fick_err_col_pos] = system.fickian_fit_errors[step_idx].err
        sheet[start_row + TSAHelper.relative_fit_error_row_pos, start_col + TSAHelper.relative_bh_col_pos] = system.bh_fit_errors[step_idx].val
        sheet[start_row + TSAHelper.relative_fit_error_row_pos, start_col + TSAHelper.relative_bh_err_col_pos] = system.bh_fit_errors[step_idx].err
        sheet[start_row + TSAHelper.relative_fit_error_row_pos, start_col + TSAHelper.relative_mbh_col_pos] = system.mbh_fit_errors[step_idx].val
        sheet[start_row + TSAHelper.relative_fit_error_row_pos, start_col + TSAHelper.relative_mbh_err_col_pos] = system.mbh_fit_errors[step_idx].err

        # place resampled sorption
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_output_resampled_time_col, dim=1] = system.resampled_steps[step_idx].time
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_output_resampled_sorption_col, dim=1] = system.resampled_steps[step_idx].dimensionlesssorption
         
        # place predictions
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_fick_col_pos, dim=1] = [item.val for item in system.fickian_step_predictions[step_idx]]
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_fick_err_col_pos, dim=1] = [item.err for item in system.fickian_step_predictions[step_idx]]
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_bh_col_pos, dim=1] = [item.val for item in system.bh_step_predictions[step_idx]]
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_bh_err_col_pos, dim=1] = [item.err for item in system.bh_step_predictions[step_idx]]
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_mbh_col_pos, dim=1] = [item.val for item in system.mbh_step_predictions[step_idx]]
        sheet[start_row + TSAHelper.relative_data_row_start + 1,  start_col + TSAHelper.relative_mbh_err_col_pos, dim=1] = [item.err for item in system.mbh_step_predictions[step_idx]]
    end
    return nothing
end

function processtemplate(::TransientSorptionApparatus, template_path::String, results_path::Union{String, Nothing}; 
    sheet_name=TSAHelper.default_sheet_name, overwrite=false, verbose=false)  # read a template and write it out with options
    # load and calculate the system
    system = TransientSorptionSystem(template_path; verbose)
    if isnothing(results_path) return system end
    # copy the template (this is the only time we will not force by default)
    if template_path == results_path
        @warn("The transient sorption template and results file paths were the same")
        return # not really dealing with this behavior yet, it would directly modify the file
    end
    cp(template_path, results_path; force=overwrite)
    
    # open the results file and start adding in the calculations done
    XLSX.openxlsx(results_path, mode="rw") do xf
        sheet = xf[sheet_name]
            write_transient_sorption_system_to_sheet(system, sheet)
    end
    return "Transient sorption template processed sucessfully."
end
function processtemplate(::TransientSorptionApparatus, template_path::String; kwargs...)
    return processtemplate(TransientSorptionApparatus(), template_path, nothing; kwargs...)
end

function is_template_valid(::TransientSorptionApparatus, args...; kwargs...)  # super duper inefficient, but should be "foolproof" in that it just reads the template into a setup
    try 
        readtemplate(TransientSorptionApparatus(), args...; kwargs...)
        return true
    catch
        return false   
    end
end
