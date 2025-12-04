module CholeskySerial
export cholesky_serial

function cholesky_serial(A)
    n = size(A, 1)
    L = zeros(n, n)

    for i in 1:n
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

end
