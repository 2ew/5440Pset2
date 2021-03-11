function calculate_transcriptional_control_array(t::Float64,x::Array{Float64,1}, problem::Dict{String,Any})::Float64

    # initialize -
    u_variable = 0.0
    
    # alias elements of the species vector -
    mRNA = x[1]
    G = x[2]
    σ70 = x[3]

    # get stuff from the problem dictionary -
    E1 = problem["E1"] # Units: kJ/mol
    E2 = problem["E2"] # Units: kJ/mol
    R = problem["ideal_gas_constant_R"] # Units: kJ/mol-K
    T_K = problem["temperature_K"] # Units: K
    KD = problem["inducer_dissociation_constant"]*1000 # Units: nM
    n = problem["inducer_cooperativity_parameter"] # Units: dimensionless

    # TODO: write u-varible function here 
    I = σ70 #concentration of Inducer (nM)
    fI = I^n/(KD^n+I^n) # Units: dimensionless
    dG1 = 41.94 # Units: kj/mol
    W1 = exp(-dG1/R/T_K)
    dG2 = -27.67 # Units: kJ/mol
    W2 = exp(-dG2/R/T_K) # Units: dimensionless
    u_variable = (W1 + W2*fI)/(1 + W1 + W2*fI) # Units: dimensionless

    # return -
    return u_variable
end