function TransposetoAxes!(X::Vector{UInt}, b::Int, n::Int)
    ccall(
        (:TransposetoAxes, "test/libskilling"),
        Cvoid,
        (Ref{UInt}, Int, Int),
        X, b, n
        )
end

function AxestoTranspose!(X::Vector{UInt}, b::Int, n::Int)
    ccall(
        (:AxestoTranspose, "test/libskilling"),
        Cvoid,
        (Ref{UInt}, Int, Int),
        X, b, n
        )
end

function InverseUndo!(X::Vector{UInt32}, b::Int32, n::Int32)
    ccall(
        (:InverseUndo, "test/libskilling"),
        Cvoid,
        (Ref{UInt32}, Int32, Int32),
        X, b, n
        )
end
n = 3
X = UInt32[5, 10, 20]
InverseUndo!(X, Int32(5), Int32(3))
@show X
AxestoTranspose!(X, 5, 3)
@assert X == UInt64[0xe, 0xf, 0x14]
TransposetoAxes!(X, 5, 3)
@assert X == UInt[5, 10, 20]


t = UInt8[0b1100, 0b0110, 0b0011]
b = 4
n = length(t)
h = interleave_transpose(t, b, n)
h == 0b001011110100

X = UInt32[5, 10, 20]
X0 = copy(X)
axes_to_transpose!(X, 5, 3)
h = interleave_transpose(X, 5, 3)
transpose_to_axes!(X, 5, 3)
@assert X == X0
