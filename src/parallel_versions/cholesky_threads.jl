module CholeskyThreads

using Base.Threads


function cholesky_threads_outer_safe(A)
    n = size(A, 1)
    L = zeros(n, n)
    
    for j in 1:n
        s = 0.0
        for k in 1:j-1
            s += L[j, k]^2
        end
        L[j, j] = sqrt(A[j, j] - s)
        
        @threads for i in (j+1):n
            s = 0.0
            for k in 1:j-1
                s += L[i, k] * L[j, k]
            end
            L[i, j] = (A[i, j] - s) / L[j, j]
        end
    end
    
    return L
end

function cholesky_threads_outer(A)
    n = size(A, 1)
    L = zeros(n, n)

    @threads for i in 1:n
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

function cholesky_threads_inner_safe(A)
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
                sum_parts = zeros(nthreads())
                @threads for tid in 1:nthreads()
                    local_sum = 0.0
                    chunk_size = ceil(Int, (j-1) / nthreads())
                    start_k = (tid-1) * chunk_size + 1
                    end_k = min(tid * chunk_size, j-1)
                    
                    for k in start_k:end_k
                        local_sum += L[i, k] * L[j, k]
                    end
                    sum_parts[tid] = local_sum
                end
                s = sum(sum_parts)
            end
            
            L[i, j] = (A[i, j] - s) / L[j, j]
        end
    end
    
    return L
end

function cholesky_threads_inner(A)
    n = size(A, 1)
    L = zeros(n, n)

    for i in 1:n
        @threads for j in 1:i

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
