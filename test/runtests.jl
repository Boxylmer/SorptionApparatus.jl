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
    generatetemplate(GasSorptionApparatus(), joinpath(@__DIR__, "template_results", "generated gas template (empty).xlsx"))
    gas_system = processtemplate(
        GasSorptionApparatus(), 
        joinpath(@__DIR__, "test_templates", "TPBO-0.75_CO2_27C.xlsx"), 
        joinpath(@__DIR__, "template_results", "TPBO-0.75_CO2_27C_results.xlsx"); 
        overwrite=true
    ) 
    @test concentration(gas_system.isotherm)[end].val ≈ 135.74857876430664

    # Transient Sorption Apparatus
    generatetemplate(TransientSorptionApparatus(), joinpath(@__DIR__, "template_results", "generated transient template (empty).xlsx"); standalone=true)
end

nothing