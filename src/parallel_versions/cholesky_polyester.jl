module CholeskyPolyester

using Polyester

function cholesky_polyester_outer_safe(A)
    n = size(A, 1)
    L = zeros(n, n)
    
    @batch for i in 1:n
        for j in 1:i

            while j < i && L[j, j] == 0
                yield()
            end
            
            s = 0.0
            for k in 1:j-1
                while L[i, k] == 0 || L[j, k] == 0
                    yield()
                end
                s += L[i, k] * L[j, k]
            end
            
            if i == j
                L[i, j] = sqrt(A[i, i] - s)
            else
                L[i, j] = (A[i, j] - s) / L[j, j]
            end
            
        end
    end
    
    return L
end

function cholesky_polyester_outer(A)
    n = size(A, 1)
    L = zeros(n, n)

    @batch for i in 1:n
        for j in 1:i
            
            s = 0.0
            for k in 1:j-1
                s += L[i,k] * L[j,k]
            end

            if i == j
                L[i,j] = sqrt(A[i,i] - s)
            else
                L[i,j] = (A[i,j] - s) / L[j,j]
            end

        end
    end

    return L
end

function cholesky_polyester_inner_safe(A)
    n = size(A, 1)
    L = zeros(n, n)
    
    for j in 1:n
        s_diag = 0.0
        for k in 1:j-1
            s_diag += L[j, k]^2
        end
        L[j, j] = sqrt(A[j, j] - s_diag)
        
        for i in (j+1):n
            s = 0.0
            
            if j > 1
                chunk_size = max(1, ceil(Int, (j-1) / 100))
                n_chunks = ceil(Int, (j-1) / chunk_size)
                
                partial_sums = zeros(n_chunks)
                
                @batch for chunk in 1:n_chunks
                    start_k = (chunk-1) * chunk_size + 1
                    end_k = min(chunk * chunk_size, j-1)
                    
                    local_sum = 0.0
                    for k in start_k:end_k
                        local_sum += L[i, k] * L[j, k]
                    end
                    partial_sums[chunk] = local_sum
                end
                
                s = sum(partial_sums)
            end
            
            L[i, j] = (A[i, j] - s) / L[j, j]
        end
    end
    
    return L
end

function cholesky_polyester_inner(A)
    n = size(A, 1)
    L = zeros(n, n)

    for i in 1:n
        @batch for j in 1:i

            s = 0.0
            for k in 1:j-1
                s += L[i,k] * L[j,k]
            end

            if i == j
                L[i,j] = sqrt(A[i,i] - s)
            else
                L[i,j] = (A[i,j] - s) / L[j,j]
            end

        end
    end
    
    return L
end

end
