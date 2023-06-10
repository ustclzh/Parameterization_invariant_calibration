function y = mse_pde(L,data)
D = invparametrization(L,data.parameterization);
y = simulator(data,D);
y = log(mean((y-data.obs).^2)^(1/2));
end