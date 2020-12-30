# This adds mutators to the Vimes.jl library.
# It uses this as a guide:
# http://pitest.org/quickstart/mutators/

using SourceWalk
using Vimes
egtext = """
function aa(x)
    if x > 3
        7
    end
end
"""
egi = quote
    if x > 3
        7
    end
    54
end
# Checks whether a library replacement works on a given expression.
function show_replacements(ex, fs)
    match_any = false
    Vimes.pathwalk(ex) do p, x
      for f in fs
        if matches(f, x)
          println("starts as $(x)")
          println("matches to $(f(x))")
          match_any = true
        end
      end
      return x
    end
    match_any
end

egi = SourceWalk.SourceFile("nowhere", egtext)
SourceWalk.replacement(egi, flipcond)
SourceWalk.textmap(_ ->Expr(:file, :(b*2)), "a * 2")
SourceWalk.textmap(nn_cond(flipcond), egtext)
show_replacements(egi, [flipcond])

# Turns > into >=, >= into >. The same for less-than.
function conditionals_boundary(x)
    if isexpr(x, :call)
        if x.args[1] == :>
            Expr(x.head, :>=, x.args[2:end]...)
        elseif x.args[1] == :>=
            Expr(x.head, :>, x.args[2:end]...)
        elseif x.args[1] == :<
            Expr(x.head, :<=, x.args[2:end]...)
        elseif x.args[1] == :<=
            Expr(x.head, :<, x.args[2:end]...)
        else
            nothing
        end
    else
        nothing
    end
end
show_replacements(egi, [conditionals_boundary])

neg_example = quote
    a = 3
    b = -a
    c = 7 - b
end

# This takes the unary - to a positive.
function invert_negatives(x)
    isexpr(x) && x.args[1] == :- && length(x.args) == 2 || return
    x.args[2]
end
show_replacements(neg_example, [invert_negatives])

# Turns + -> -, - -> +, etc, for integer and float.
function math_mutator(x)
    isexpr(x) && length(x.args) == 3 || return
    a, b = x.args[2:3]
    if x.args[1] == :-
        :(isa($a, Real) && isa($b, Real) ? $a + $b : $x)
    elseif x.args[1] == :+
        :(isa($a, Real) && isa($b, Real) ? $a - $b : $x)
    elseif x.args[1] == :*
        :(isa($a, Real) && isa($b, Real) ? $a / $b : $x)
    elseif x.args[1] == :/
        :(isa($a, Real) && isa($b, Real) ? $a * $b : $x)
    elseif x.args[1] == :%
        :(isa($a, Real) && isa($b, Real) ? $a * $b : $x)
    elseif x.args[1] == :*
        :(isa($a, Real) && isa($b, Real) ? $a % $b : $x)
    elseif x.args[1] == :^
        :(isa($a, InteRealger) && isa($b, Real) ? $a * $b : $x)
    elseif x.args[1] == :&
        :(isa($a, Integer) && isa($b, Integer) ? $a | $b : $x)
    elseif x.args[1] == :|
        :(isa($a, Integer) && isa($b, Integer) ? $a & $b : $x)
    elseif x.args[1] == :<<
        :(isa($a, Integer) && isa($b, Integer) ? $a >> $b : $x)
    elseif x.args[1] == :>>
        :(isa($a, Integer) && isa($b, Integer) ? $a << $b : $x)
    elseif x.args[1] == :>>>
        :(isa($a, Integer) && isa($b, Integer) ? $a << $b : $x)
    else
        nothing
    end
end

show_replacements(neg_example, [math_mutator])
function incode_math_mutator(a, b)
    isa(a, Real) && isa(b, Real) ? a - b : a + b
end

# This looks bad, but it totally works.
function stuff(a, b)
a + if 7 isa Real && b isa Real
    7 + b
else
    7 - b
end + 5
end
stuff(4, 2)

function has_a_return(a)
    x = 2a
    return x
end
function lacks_a_return(a)
    2a
end
function nothing_return(a)
    return
end
aa = Base.uncompressed_ast(first(methods(has_a_return)))
bb = @code_lowered(lacks_a_return(7))
bb.code[2].head == :return

return_eg = quote
    x = 2
    return x
end

calleg = :(empty_return(xy))

empty_return(x::AbstractString) = ""
empty_return(x::Real) = zero(x)
empty_return(x::Set{T}) where {T} = Set{T}()
empty_return(x::Array{T,1}) where {T} = Array{T,1}()
empty_return(x::Dict{K,V}) where {K,V} = Dict{K,V}()
empty_return(x::Char) = '\0'
empty_return(x) = x

function empty_returns(x)
    isexpr(x, :return) || return
    Expr(:return, Expr(:call, :empty_return, x.args[1]))
end
show_replacements(return_eg, [empty_returns])

flip_return(x::Bool) = !x
flip_return(x) = x

function flip_returns(x)
    isexpr(x, :return) || return
    Expr(:return, Expr(:call, :flip_return, x.args[1]))
end
show_replacements(return_eg, [flip_returns])


function return_nothing(x)
    isexpr(x, :return) && x.args[1] !== nothing || return
    Expr(:return, nothing)
end


function increment_integer(x)
    x isa Integer && !(x isa Bool) && return x + 1
end

function decrement_integer(x)
    x isa Integer && !(x isa Bool) && return x - 1
end


function scale_float(x)
    x isa AbstractFloat && return 0.9 * x
end

condeg = quote
    a && b
end

function flip_and_or(x)
    isexpr(x, :&&) && return Expr(:||, x.args...)
    isexpr(x, :||) && return Expr(:&&, x.args...)
end
