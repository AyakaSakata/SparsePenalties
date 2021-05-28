function [ diff_R ] = Rh_solve( Rh, param )

    int_max = param.int_max;
    d_int = param.d_int;
    
    int_norm = d_int/sqrt(2*pi);

    chih = param.chih;
    q0h = param.q0h;
    Qh = param.Qh;
    R = param.R;
    
    alpha = param.alpha;
    
    d2 = (-int_max:d_int:int_max).^2;
    ker = exp(-0.5*d2);

    W = 1+sqrt(Qh/(Qh-chih))*exp(-Rh+0.5*q0h*d2/(Qh-chih));

    diff_R = R - ((W-1)./W)*ker'*int_norm/alpha;

end

