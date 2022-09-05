export add_Hopf!
export Hopf_variables
export replace_Hopf_variable!

using HarmonicBalance: is_rearranged, rearrange_standard
using HarmonicBalance: LinearResponse.get_implicit_Jacobian
#import HarmonicBalance.get_steady_states; export get_steady_states

"Add a Hopf bifurcation to the system."
function add_Hopf!(eom::DifferentialEquation, Δω::Num)
    for var in get_variables(eom), ω in eom.harmonics[var]
        add_harmonic!(eom, var, ω + Δω)
        add_harmonic!(eom, var, ω - Δω)
    end
end


"Successively adds n Hopf bifurcations at the same frequency Δω"
add_Hopf!(eom::DifferentialEquation; frequency::Num, multiplicity::Int) = [add_Hopf!(eom, frequency) for k in 1:multiplicity]


"""
    Hopf_variables(eom, freq)
"""
Hopf_variables(eom::HarmonicEquation, freq::Num) = filter(x -> any(isequal.(freq, get_all_terms(x.ω))), eom.variables)


"""
    Obtain the Jacobian of `eom` with a free (Hopf) variable `fixed_var`.
    `fixed_var` marks the variable which is free due to U(1) symmetry. Its entry in the Jacobian matrix must be discarded
    (the discarded eigenvalue is 0, corresponding to a free phase.)
"""
function _Hopf_Jacobian(eom::HarmonicEquation, fixed_var::HarmonicVariable)
    eom_Jac = rearrange_standard(eom)
    free_idx = findall(x -> isequal(x, fixed_var), eom_Jac.variables)  
    deleteat!(eom_Jac.equations, free_idx)
    deleteat!(eom_Jac.variables, free_idx)

    # the free variable can be fixed to zero, this is also done in the corresponding HarmonicEquation later
    eom_Jac = substitute_all(eom_Jac, fixed_var => 0)
    J = get_implicit_Jacobian(eom_Jac)
end


"""
    Construct a `Problem` in the case where U(1) symmetry is present
    due to having added a Hopf variable with frequency `ω_limit`.
"""
function _Hopf_Problem(eom::HarmonicEquation, ω_limit::Num)

    eom = deepcopy(eom) # do not mutate eom
    isempty(Hopf_variables(eom, ω_limit)) ? error("No Hopf variables found!") : nothing
    !any(isequal.(eom.parameters, ω_limit)) ? error(ω_limit, " is not a parameter of the harmonic equation!") : nothing 

    fixed_var = Hopf_variables(eom, ω_limit)[end] # eliminate one of the Cartesian variables, it does not matter which
    
    # get the Hopf Jacobian before altering anything - this is the usual Jacobian but the entries corresponding
    # to the free variable are removed
    J = _Hopf_Jacobian(eom, fixed_var)

    _fix_gauge!(eom, ω_limit, fixed_var)

    # define Problem as usual but with the Hopf Jacobian (always computed implicitly)
    p = Problem(eom; Jacobian=J)
    return p
end


function get_steady_states(eom::HarmonicEquation, swept, fixed; ω_limit, kwargs...)   
    prob = _Hopf_Problem(eom, ω_limit);
    HarmonicBalance.get_steady_states(prob, swept, fixed; random_warmup=true, threading=true, kwargs...)
end


function _fix_gauge!(eom::HarmonicEquation, ω_limit::Num, fixed_var::HarmonicVariable)

    new_symbol = HarmonicBalance.declare_variable(string(ω_limit), first(get_independent_variables(eom)))
    rules = Dict(ω_limit => new_symbol, fixed_var.symbol => Num(0))
    eom.equations = expand_derivatives.(substitute_all(eom.equations, rules))
    eom.parameters = setdiff(eom.parameters, [ω_limit]) # ω_limit is now NOT a parameter anymore
    
    fixed_var.type = "Hopf"
    fixed_var.ω = Num(0)
    fixed_var.symbol = new_symbol
    fixed_var.name = var_name(new_symbol)
end

