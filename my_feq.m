function feq = my_feq(i, rho, ux, uy)

    global w;
    global cs2;
    global cx;
    global cy;

    cdotu = (cx(i)*ux + cy(i)*uy)/cs2;
    feq = rho.*w(i).*(1 + cdotu + 0.5*cdotu.*cdotu - 0.5.*(ux.^2 + uy.^2)/cs2);
