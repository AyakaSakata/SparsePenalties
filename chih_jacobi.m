function [ jacobi ] = chih_jacobi( chih, param )

    d = 1e-3*chih;
    
    jacobi = (chih_solve(chih+d, param)-chih_solve(chih, param))/d;

end

