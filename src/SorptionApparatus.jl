module SorptionApparatus
    using MembraneBase
    using MembraneEOS
    using SorptionModels
    using XLSX
    using Measurements
    using ProgressMeter

    include(joinpath("Apparatuses", "ApparatusFuncs.jl"))
    export generatetemplate
    export readtemplate
    export processtemplate
    export savetemplate

    include(joinpath("Apparatuses", "TransientSorptionApparatus.jl"))
    export TransientSorptionApparatus
    include(joinpath("Apparatuses", "VaporSorptionApparatus.jl"))
    export VaporSorptionApparatus
    include(joinpath("Apparatuses", "GasSorptionApparatus.jl"))
    export GasSorptionApparatus

    include(joinpath("ApparatusAnalyses.jl"))
end
