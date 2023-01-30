struct SimpleParams{FT<:Real} <: YSWParams
     c :: FT  # Gravity Wave speed
     f :: FT  # Coriolis parameter
    τd :: FT  # Linear Damping Timescale
     ν :: FT  # Hyperviscosity coefficient
    nν :: Int # Order of the hyperviscous operato
end

struct ForcingParams1D{FT<:Real} <: ForcingParams
     c :: FT  # Gravity Wave speed
     f :: FT  # Coriolis parameter
    τd :: FT  # Linear Damping Timescale
     ν :: FT  # Hyperviscosity coefficient
    nν :: Int # Order of the hyperviscous operator
    convection :: Forcing1D
    ϕforcing   :: Forcing1D
end

struct ForcingParams2D{FT<:Real} <: ForcingParams
     c :: FT  # Gravity Wave speed
     f :: FT  # Coriolis parameter
    τd :: FT  # Linear Damping Timescale
     ν :: FT  # Hyperviscosity coefficient
    nν :: Int # Order of the hyperviscous operator
    convection :: Forcing2D
    ϕforcing   :: Forcing2D
end

function DefineParams(
    G :: AbstractGrid;
    c :: Real = 20,      # units in m s⁻¹
    f :: Real = 0,       # units in s⁻¹
   τd :: Real = 1,       # units in days
    ν :: Real = 0,
   nν :: Int  = 1,
   τc :: Real = 1,       # units in hours
   ϕc :: Real = c^2,     # units in m² s²
   rc :: Real = 10,      # units in km
   Sc :: Real = 5.e-14, # units in m**-2 s**-1
   τl :: Real = 1,       # units in days
   ϕ0 :: Real = 0,       # units in m² s²
   Fl :: Real = 0,
   convection :: Bool = false,
   wtg        :: Bool = false,
   ϕforcing   :: Bool = false,
   convectionfunction! :: Function = calcYangConvection!,
   ϕfunction!          :: Function = calcYangLargeScale,
)

   if convection || ϕforcing
       P = DefineParams(
           G,c,f,τd*86400,ν,nν,
           τc*3600,ϕc,rc*1000,Sc,
           τl*86400,ϕ0,Fl,wtg,
           convectionfunction!,ϕfunction!
       )
   else; P = SimpleParams{eltype(G)}(c,f,τd*86400,ν,nν)
   end

end

DefineParams(
     G :: OneDGrid,
     c :: Real,
     f :: Real,
    τd :: Real,
     ν :: Real,
    nν :: Int ,
    τc :: Real,
    ϕc :: Real,
    rc :: Real,
    Sc :: Real,
    τl :: Real,
    ϕ0 :: Real,
    Fl :: Real,
    wtg                :: Bool,
    convectionfunction :: Function,
    ϕfunction          :: Function
) = ForcingParams1D{eltype(G)}(
    c,f,τd,ν,nν,
    GenerateConvection(τc,ϕc,rc,Sc,convectionfunction,G),
    GenerateϕForcing(τl,ϕ0,Fl,wtg,ϕfunction,G)
)

DefineParams(
     G :: TwoDGrid,
     c :: Real,
     f :: Real,
    τd :: Real,
     ν :: Real,
    nν :: Int ,
    τc :: Real,
    ϕc :: Real,
    rc :: Real,
    Sc :: Real,
    τl :: Real,
    ϕ0 :: Real,
    Fl :: Real,
    wtg                :: Bool,
    convectionfunction :: Function,
    ϕfunction          :: Function
) = ForcingParams2D{eltype(G)}(
    c,f,τd,ν,nν,
    GenerateConvection(τc,ϕc,rc,Sc,convectionfunction,G),
    GenerateϕForcing(τl,ϕ0,Fl,wtg,ϕfunction,G)
)