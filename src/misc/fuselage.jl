#Fuselage class
@kwdef mutable struct Fuselage
    # General Fuselage Properties
    weight::Float64 = 0.0
    sigma_factor::Float64 = 1.0
    moment::Float64 = 0.0

    # Layout
    layout::FuselageLayout = FuselageLayout()

    # Structures
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    #to do: Add material-> Structural members
    skin::StructuralMember = StructuralMember(σ=1.0342136834273669e8, ρ=2700)
    shell::StructuralMember = StructuralMember(σ=2.0684273668547338e8, ρ=2700) # IS just Skin + Additional
    cone::StructuralMember = StructuralMember()
    floor::StructuralMember = StructuralMember(thickness=0.127)
    insulation::StructuralMember = StructuralMember()
    window::StructuralMember = StructuralMember()
    bending_h::StructuralMember = StructuralMember(σ=2.0684273668547338e8, ρ=2700)
    bending_v::StructuralMember = StructuralMember(σ=2.0684273668547338e8, ρ=2700)

    # Loads
    
    # Misc properties
    Nland::Float64 = 6.0
    n_decks::Float64 = 1
    fstring::Float64 = 0.35
    fframe::Float64 = 0.25
    ffadd::Float64 = 0.2
    nftanks::Int64 = 1
    rEshell::Float64 = 1.0
    # tank_placement
end


fuselage = Fuselage()
Eskin = fuselage.material.E #parg[igEcap]
Ebend = Eskin * fuselage.rEshell
Gskin = Eskin * 0.5 / (1.0 + 0.3)

fusew!(fuselage,13344.666000000001, 219964.5779, 76987.602265, 21996.45779, 7698.7602265000005, 0.0, 0.0, 0.0, 0.0, 0.0, "", 54911.323281976234, 435.0, 22.0, 60.0, 8607.309570000001, 8607.309570000001, 0.4, 0.7, 439929.1558, 439929.15580000007, 5.1591026057637395, 0.3, 1.0, 34.8996, 33.528, 20.165065369407014, 17.3736, 0.0, 2.1336, 36.576, 15.8496, 0.0,Eskin,Ebend,Gskin)