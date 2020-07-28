module SimplePolynomials

using MultivariatePolynomials
const MP = MultivariatePolynomials

# Taken from DynamicPolynomials/src/var.jl
function polyvecvar(prefix, idxset)
    [Variable("$(prefix * string(i))") for i in idxset]
end

function buildpolyvar(var)
    if isa(var, Symbol)
        :($(esc(var)) = Variable($"$var"))
    else
        isa(var, Expr) || error("Expected $var to be a variable name")
        Base.Meta.isexpr(var, :ref) || error("Expected $var to be of the form varname[idxset]")
        length(var.args) == 2 || error("Expected $var to have one index set")
        varname = var.args[1]
        prefix = string(var.args[1])
        idxset = esc(var.args[2])
        :($(esc(varname)) = polyvecvar($prefix, $idxset))
    end
end

macro polyvar(args...)
    reduce((x,y) -> :($x; $y), :(), buildpolyvar.(args))
end


struct Variable <: AbstractVariable
    name::String
end
MP.name(v::Variable) = v.name
Base.:(==)(v1::Variable, v2::Variable) = v1.name == v2.name
struct Monomial <: AbstractMonomial
    vars::Vector{Variable}
    exps::Vector{Int}
end
const MonomialLike = Union{Variable, Monomial}
MP.variables(m::Monomial) = m.vars
MP.exponents(m::Monomial) = m.exps
MP.monomial(v::Variable) = Monomial([v], [1])
MP.termtype(::Union{MonomialLike, Type{<:MonomialLike}}) = Term{Int}
struct Term{T} <: AbstractTerm{T}
    α::T
    m::Monomial
end
MP.coefficient(t::Term) = t.α
MP.monomial(t::Term) = t.m
struct Polynomial{T} <: AbstractPolynomial{T}
    terms::Vector{Term{T}}
end
MP.terms(p::Polynomial) = p.terms

end
