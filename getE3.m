function E3 = getE3(alpha, gam, lambda, epsValue)
    E3 = -lambda * (getPsi(alpha + gam, epsValue) - alpha - gam) / 2;
end