function [ diff_chih ] = chih_solve( chih, param )

    alpha = param.alpha;
    rho = param.rho;
    erfc_n = erfc(1/sqrt(2*chih));
    
    chih_n = (chih+1)*erfc_n-2*sqrt(chih)/sqrt(2*pi)*exp(-0.5/chih);
    diff_chih = chih-((1-rho)*chih_n+rho*(chih+1))/alpha;
end

