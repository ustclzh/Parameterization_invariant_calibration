function y=log_posterior_simulator(data,D,V,v)
d=size(D,1);
y_s=simulator(data,D);
Rss=(y_s-data.obs)*((y_s-data.obs)')/data.sigmasq;
if nargin==2
    prior=wishpdf(D,eye(d),d);
else
    prior=wishpdf(D,V,v);
end
y=-Rss/2+log(prior);
end