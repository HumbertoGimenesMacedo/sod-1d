function sigma = getSigma(z, epsValue, lamb)
    sigma = (getPsi(z, epsValue) + lamb * z^2) / 2;
end