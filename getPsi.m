function psiValue = getPsi(z, eps_psi)
    if (abs(z) < eps_psi)
        psiValue = (z^2 / eps_psi + eps_psi) / 2;
    else
        psiValue = abs(z);
    end
end