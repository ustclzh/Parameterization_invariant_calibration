function y=log_posterior_simulator(data,D,V,v)
if min(eig(D)<=0)
    y = -inf;
    return
end
y_s=simulator(data,D);
Rss=sum((y_s-data.obs).^2)/data.sigmasq;
prior=Prior_Density_PDT(D,V,v);
y=-Rss/2+log(prior);
end