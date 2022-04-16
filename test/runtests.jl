using SorptionApparatus
using Test
using Revise

@testset "SorptionApparatus.jl" begin
    # Vapor Sorption Apparatus
    processtemplate(
        VaporSorptionApparatus(), 
        joinpath(@__DIR__, "test_templates", "methanol_template_with_transients.xlsx"), 
        joinpath(@__DIR__, "template_results", "methanol_results_with_transients.xlsx"); 
        overwrite=true
    ) 

end
