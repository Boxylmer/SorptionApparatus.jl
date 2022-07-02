using SorptionApparatus
using MembraneBase
using Test
using Revise

@testset "SorptionApparatus.jl" begin
    # Vapor Sorption Apparatus
    generatetemplate(VaporSorptionApparatus(), joinpath(@__DIR__, "template_results", "generated vapor template (empty).xlsx"))
    vapor_system = processtemplate(
        VaporSorptionApparatus(), 
        joinpath(@__DIR__, "test_templates", "methanol_template_with_transients.xlsx"), 
        joinpath(@__DIR__, "template_results", "methanol_results_with_transients.xlsx"); 
        overwrite=true
    ) 
    
    @test concentration(vapor_system.isotherm)[end].val ≈ 120.99026680342382

    # Gas Sorption Apparatus
    generated_gas_template_path = joinpath(@__DIR__, "template_results", "generated gas template (empty).xlsx")
    tpbo_75_co2_27c_template_path = joinpath(@__DIR__, "test_templates", "TPBO-0.75_CO2_27C.xlsx")
    tpbo_75_co2_27c_resaved_template_path = joinpath(@__DIR__, "template_results", "TPBO-0.75_CO2_27C_saved_template.xlsx")
    tpbo_75_co2_27c_result_path = joinpath(@__DIR__, "template_results", "TPBO-0.75_CO2_27C_results.xlsx")

    generatetemplate(GasSorptionApparatus(), generated_gas_template_path)
    gas_system = processtemplate(
        GasSorptionApparatus(), 
        tpbo_75_co2_27c_template_path, 
        tpbo_75_co2_27c_result_path
        ; 
        overwrite=true
    ) 
    @test concentration(gas_system.isotherm)[end].val ≈ 135.74857876430664

    # do saved templates avoid making changes to the template values?
    gas_setup = readtemplate(GasSorptionApparatus(), tpbo_75_co2_27c_template_path)
    savetemplate(gas_setup, tpbo_75_co2_27c_resaved_template_path)
    gas_system_reread = processtemplate(GasSorptionApparatus(), tpbo_75_co2_27c_resaved_template_path)
    @test concentration(gas_system.isotherm)[end].val == concentration(gas_system_reread.isotherm)[end].val
    @test concentration(gas_system.isotherm)[end].err == concentration(gas_system_reread.isotherm)[end].err
    @test concentration(gas_system.isotherm)[1].val == concentration(gas_system_reread.isotherm)[1].val
    @test concentration(gas_system.isotherm)[1].err == concentration(gas_system_reread.isotherm)[1].err

    # Transient Sorption Apparatus
    generatetemplate(TransientSorptionApparatus(), joinpath(@__DIR__, "template_results", "generated transient template (empty).xlsx"); standalone=true)
end

nothing