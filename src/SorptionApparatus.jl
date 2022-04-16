module SorptionApparatus
    using MembraneBase
    using MembraneEOS
    using SorptionModels
    using XLSX

    include(joinpath("Apparatuses", "TransientSorptionApparatus.jl"))
    include(joinpath("Apparatuses", "VaporSorptionApparatus.jl"))
    # include(joinpath("Apparatuses", "GasSorptionApparatus.jl"))

    include(joinpath("ApparatusAnalyses.jl"))
end
