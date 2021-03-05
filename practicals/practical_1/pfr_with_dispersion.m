function yy = pfr_with_dispersion(y, L, v, Pe, k, Np, NEQ, y0, delta_z_star, nu)

for j = 1:NEQ
    for i = 1:Np
        
        c(i, j) = y((i -1)*NEQ + j);
        
    end
end


for j = 1:NEQ
    for i = 2:Np-1
        
        c_p(i, j) = (c(i, j) - c(i-1, j))/delta_z_star;
        c_pp(i, j) = (c(i+1, j) - 2*c(i,j) + c(i-1, j))/delta_z_star^2;
               
    end
    
    i = 1;
    c_p(i, j) = (c(i+1, j) - c(i, j))/delta_z_star;
    
    i = Np;
    c_p(i, j) = (c(i, j) - c(i-1, j))/delta_z_star;
    
    
end

for j = 1:NEQ
    i = 1;
    yy((i - 1)*NEQ + j) = y0(j) - (c(i,j) - (1/Pe(j))*c_p(i, j));
    
    for i = 2:Np -1
        yy((i-1)*NEQ + j) = -c_p(i, j) + 1/Pe(j)*c_pp(i,j) + nu(j)*L/v*k*c(i, 1);
    end
    
    i = Np;
    yy((i-1)*NEQ + j) = c_p(i, j);

end