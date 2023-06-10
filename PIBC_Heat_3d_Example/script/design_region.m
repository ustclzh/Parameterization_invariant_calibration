function range = design_region(data)
d = data.d;
U = data.eig_U;
L = data.eig_L;
para = data.parameterization;
N = d*(d+1)/2;
range = ones(2,N);
range(2,:) = range(2,:) * U^(1/2);
range(1,:) = -range(1,:) * U^(1/2);
for i = 1:d
    range(1,i*(i+1)/2) = (L/U)^(1/2);
end
end