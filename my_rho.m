function [rho, ux, uy] = my_rho(rho, ux, uy)
    global fin;
    global cx;
    global cy;
    global nw;

    rho = sum(fin);
    ux(:, :, :) = 0;
    uy(:, :, :) = 0;
    for i = 1:nw
        ux = ux + fin(i,:,:) * cx(i);
        uy = uy + fin(i,:,:) * cy(i);
    end
    ux = ux./rho;
    uy = uy./rho;

endfunction
