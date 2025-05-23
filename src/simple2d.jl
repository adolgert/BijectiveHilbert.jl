module Bijective
import ..BijectiveHilbert: encode_hilbert_zero, decode_hilbert_zero!, HilbertAlgorithm, axis_type, log_base2


struct Simple2D{T} <: HilbertAlgorithm{T}
end


Simple2D(::Type{T}) where {T} = Simple2D{T}()


axis_type(::Simple2D{T}) where {T} = Int


# The paper used the floor(log2()) code, so this replaces it.
rmin_orig(x, y) = convert(Int, floor(log2(max(x, y))) + 1)
rmin_replace(x, y) = log_base2(x | y) + one(Int)


function encode_hilbert_zero(::Simple2D{T}, X::AbstractVector{A})::T where {A, T}
    x = X[1]
    y = X[2]
    z = zero(T)
    if x == zero(A) && y == zero(A)
        return z
    end
    rmin = rmin_replace(x, y)
    w = one(A) << (rmin - 1)
    while rmin > 0
        z <<= 2
        if rmin&1 == 1  # odd
            if x < w
                if y >= w
                    # 1
                    x, y = (y - w, x)
                    z += one(T)
                end  # else x, y remain the same.
            else
                if y < w
                    x, y = ((w<<1) - x - one(w), w - y - one(w))
                    # 3
                    z += T(3)
                else
                    x, y = (y - w, x - w)
                    z += T(2)
                    # 2
                end
            end
        else  # even
            if x < w
                if y >= w
                    # Quadrant 3
                    x, y = (w - x - one(w), (w << 1) - y - one(w))
                    z += T(3)
                end  # else do nothing for quadrant 0.
            else
                if y < w
                    # 1
                    x, y = (y, x-w)
                    z += one(T)
                else
                    # 2
                    x, y = (y-w, x-w)
                    z += T(2)
                end
            end
        end
        rmin -= 1
        w >>= 1
    end
    z
end


function decode_hilbert_zero!(::Simple2D{T}, X::AbstractVector{A}, z::T) where {A,T}
    r = z & T(3)
    x = (r & 2) >> 1
    y = ((r + 1) & 2) >> 1
    z >>= 2
    rmin = 2
    w = one(A) << 1
    while z > zero(T)
        r = z & T(3)
        parity = 2 - rmin&1
        if rmin & 1 != 0
            # Nothing to do for quadrant 0.
            if r == 1
                x, y = (y, x + w)
            elseif r == 2
                x, y = (y + w, x + w)
            elseif r == 3
                x, y = ((w << 1) - x - one(A), w - y - one(A))
            end
        else
            if r == 1
                x, y = (y + w, x)
            elseif r == 2
                x, y = (y + w, x + w)
            elseif r == 3
                x, y = (w - x - one(A), (w << 1) - y - one(A))
            end
        end
        z >>= 2
        rmin += 1
        w <<= 1
    end
    X[1] = x
    X[2] = y
end

end
