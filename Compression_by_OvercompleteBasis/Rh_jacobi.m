function [ Rh_der ] = Rh_jacobi( Rh, param )

    d = Rh*1e-3;

    Rh_der = (Rh_solve(Rh+d,param) - Rh_solve(Rh, param))/d;

end

