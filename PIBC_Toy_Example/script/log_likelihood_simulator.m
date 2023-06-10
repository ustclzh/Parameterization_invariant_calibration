function y=log_likelihood_simulator()
    y_s=simulator(data,D);
    Rss=(y_s-data.obs)*((y_s-data.obs)')/data.sigmasq;
    y=-Rss/2;
end