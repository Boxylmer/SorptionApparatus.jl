module SorptionApparatus
    using MembraneBase
    using MembraneEOS
    using SorptionModels
    using XLSX
    using Measurements

    include(joinpath("Apparatuses", "ApparatusFuncs.jl"))
    export generatetemplate
    export readtemplate
    export processtemplate

    include(joinpath("Apparatuses", "TransientSorptionApparatus.jl"))
    export TransientSorptionApparatus
    include(joinpath("Apparatuses", "VaporSorptionApparatus.jl"))
    export VaporSorptionApparatus
    include(joinpath("Apparatuses", "GasSorptionApparatus.jl"))
    export GasSorptionApparatus

    include(joinpath("ApparatusAnalyses.jl"))
end
