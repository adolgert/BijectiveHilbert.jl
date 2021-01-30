
using CSV


function filename_to_hilbert_dims(fn)
   digits = parse.(Int, split(splitext(basename(fn))[1], "_"))
   popfirst!(digits), digits
end


function hamilton_example(fn)
    n, ms = filename_to_hilbert_dims(fn)
    cnt = prod(1 << x for x in ms)
    h = zeros(UInt32, cnt)
    X = zeros(UInt8, n, cnt)
    for (i, row) in enumerate(CSV.File(fn))
        h[i] = convert(UInt32, row[1])
        for j in 1:length(ms)
            X[j, i] = convert(UInt8, row[1 + j])
        end
    end
    hsort = sortperm(h)
    h[hsort], X[:, hsort], ms
end


function steps_are_one_away(h, X)
    n = size(X, 1)
    Y = copy(X)
    Y .= zero(UInt8)
    for i in eachindex(h)
        tdiff = 0
        for j in 1:n
            if X[j, i] > Y[j, i]
                tdiff += X[j, i] - Y[j, i]
            else Y[j] > X[j]
                tdiff += Y[j, i] - X[j, i]
            end
        end
        if tdiff > 1
            @show tdiff, i, X[:, i], Y[:, i]
        end
        Y .= X
    end
end
