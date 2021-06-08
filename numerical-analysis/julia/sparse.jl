function convertToMapRepresentation(A)
    map = Dict()
    # t = size(A)
    # n = t[1]
    # m = t[2]
    for (i, row) in enumerate(A)
        for (j, val) in enumerate(row)
            if A[i][j] != 0
                map[(i, j)] = val
            end
        end
    end
    return map
end

function sparseMultiplication(A, B)
    DictA = convertToMapRepresentation(A) 
    DictB = convertToMapRepresentation(B)
    DictC = Dict()
    for (key, val) in DictA
        for (key2, val2) in DictB
            if key[2] == key2[1]
                DictC[(key[1], key2[2])] = val * val2
            end
        end
    end
    return DictC    
end