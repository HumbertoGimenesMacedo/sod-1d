function E2 = getE2(alpha_m, gam_m, alpha, gam, lambda, epsValue)
    
    C_minus = lambda * (getPsi(alpha + gam, epsValue) - alpha - gam) / 2;
    C_plus = lambda * (getPsi(alpha_m + gam_m, epsValue) + alpha_m + gam_m) / 2;
    E2 = 1 + lambda * (C_minus + C_plus); 
    
end