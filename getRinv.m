function Rinv = getRinv(u, a, H)

    %b2 = (gamma - 1)/a^2;
    %b1 = b2 * u^2/2;


   % Rinv = [(b1 + u/a)/2, (-b2*u-1/a)/2, b2/2; 1-b1, b2*u, -b2; (b1-u/a)/2, (-b2*u + 1/a)/2, b2/2];
   % 
   Rinv = [(u*(- u^2 + a*u + 2*H))/(2*a*(- u^2 + 2*H)), -(- u^2 + 2*a*u + 2*H)/(2*a*(- u^2 + 2*H)),  1/(- u^2 + 2*H);
          (2*(- u^2 + H))/(- u^2 + 2*H),                        (2*u)/(- u^2 + 2*H), -2/(- u^2 + 2*H);
            (u*(u^2 + a*u - 2*H))/(2*a*(- u^2 + 2*H)),   -(u^2 + 2*a*u - 2*H)/(2*a*(- u^2 + 2*H)),  1/(- u^2 + 2*H)];

end