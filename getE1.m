function E1 = getE1(alpha, gam, lambda, epsValue)
    E1 = -lambda * (getPsi(alpha + gam, epsValue) + alpha + gam) / 2;
end