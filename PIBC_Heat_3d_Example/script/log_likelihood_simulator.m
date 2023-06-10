function y=log_likelihood_simulator(data,D)
    if min(eig(D)<=0)
    y = -inf;
    return
end
y_s=simulator(data,D);
y=sum((y_s-data.obs).^2)/data.sigmasq;
end