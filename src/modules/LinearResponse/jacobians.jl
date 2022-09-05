export get_Jacobian

"""
Here stability and linear response is treated with the slow-flow approximation (SFA), see Chapter 5 of JK's thesis.
Linear response always appears as a sum of Lorentzians, but is inaccurate where these are peaked far from the drive frequency.
The Jacobian is stored in the Problem object as a function that takes a solution dictionary to give the numerical Jacobian. 
"""

"""
$(SIGNATURES)

Obtain the symbolic Jacobian matrix of `eom` (either a `HarmonicEquation` or a `DifferentialEquation`).
This is the linearised left-hand side of F(u) = du/dT.

"""
function get_Jacobian(eom::HarmonicEquation)

    rearr = !HarmonicBalance.is_rearranged(eom) ? HarmonicBalance.rearrange_standard(eom) : eom
    lhs = _remove_brackets(rearr)
    vars = _remove_brackets.(eom.variables)
    
    get_Jacobian(lhs, vars)

end

" Obtain a Jacobian from a `DifferentialEquation` by first converting it into a `HarmonicEquation`. "
function get_Jacobian(diff_eom::DifferentialEquation)
    @variables T
    harmonic_eq = get_harmonic_equations(diff_eom, slow_time=T, fast_time=first(get_independent_variables(diff_eom)))
    get_Jacobian(harmonic_eq)
end


" Get the Jacobian of a set of equations `eqs` with respect to the variables `vars`. "
function get_Jacobian(eqs::Vector{Num}, vars::Vector{Num})
    length(eqs) == length(vars) || error("Jacobians are only defined for square systems!")
    M = Matrix{Num}(undef, length(vars), length(vars))

    for idx in CartesianIndices(M)
        M[idx] = expand_derivatives(d(eqs[idx[1]], vars[idx[2]]))
    end
    M
end

get_Jacobian(eqs::Vector{Equation}, vars::Vector{Num}) = get_Jacobian(Num.(getfield.(eqs, :lhs) .- getfield.(eqs, :rhs)), vars)


"""
Code folllows for an implicit treatment of the Jacobian. Usually we rearrange the linear response equations to have time-derivatives on one side.
    This may be extremely costly. Implicit evaluation means only solving the equations AFTER numerical values have been plugged in, giving
    a constant time cost per run.
"""

# for implicit evaluation, the numerical values precede the rearrangement
# for limit cycles, the zero eigenvalue causes the rearrangement to fail -> filter it out
function _get_J_matrix(eom::HarmonicEquation; order=0)

    order > 1 && error("Cannot get a J matrix of order > 1 from the harmonic equations.\nThese are by definition missing higher derivatives")
    
    vars_simp = Dict([var => HarmonicBalance.declare_variable(var_name(var)) for var in get_variables(eom)])
    T = get_independent_variables(eom)[1]
    J = get_Jacobian(eom.equations, d(get_variables(eom), T, order))
    
    expand_derivatives.(HarmonicBalance.substitute_all(J, vars_simp))
end


# COMPILE THIS?
# return a function to find the Jacobian implicitly
function get_implicit_Jacobian(eom::HarmonicEquation)::Function
    J0 = _get_J_matrix(eom, order=0)
    J1 = _get_J_matrix(eom, order=1)

    function J(soln::OrderedDict)
        -inv(ComplexF64.(substitute_all(J1, soln))) * ComplexF64.(substitute_all(J0, soln))
    end
    J
end

