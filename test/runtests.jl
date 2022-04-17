using SorptionApparatus
using Test
using Revise

@testset "SorptionApparatus.jl" begin
    # Vapor Sorption Apparatus
    # processtemplate(
    #     VaporSorptionApparatus(), 
    #     joinpath(@__DIR__, "test_templates", "methanol_template_with_transients.xlsx"), 
    #     joinpath(@__DIR__, "template_results", "methanol_results_with_transients.xlsx"); 
    #     overwrite=true
    # ) 

    # Gas Sorption Apparatus
    processtemplate(
        GasSorptionApparatus(), 
        joinpath(@__DIR__, "test_templates", "TPBO-0.75_CO2_27C.xlsx"), 
        joinpath(@__DIR__, "template_results", "TPBO-0.75_CO2_27C_results.xlsx"); 
        overwrite=true
    ) 

end

nothing