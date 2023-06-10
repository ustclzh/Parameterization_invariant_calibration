function y=simulator(data,D)
%% Coupled moisture and heat transfer equation
c=D(:);
model=data.model;
specifyCoefficients(model, 'm', 0,'d',0,'c', c, 'a', 0, 'f', 0);
applyBoundaryCondition(model,'neumann','Edge',[1,3,4],'q',0,'g',10);
applyBoundaryCondition(model,'dirichlet','Edge',[2],'h',1,'r',0);
result = solvepde(model);
y=interpolateSolution(result,data.x_obs(:),data.y_obs(:));
y=y(:)';
end
function fmatrix=fcoef(location, state)
n1 = 1;
nr = numel(location.x);
fmatrix = zeros(n1,nr);
fmatrix(1,:) = 1000*((location.x.^2+location.y.^2)<0.1).*ones(1,nr);
end
