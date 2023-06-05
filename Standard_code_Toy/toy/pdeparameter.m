function data=pdeparameter(Sigma)
gd=[1;0;0;1];
model = createpde(1);
g = decsg(gd, 'R', 'R');
geometryFromEdges(model,g);
%hmax = 0.1;
mesh = generateMesh(model,'GeometricOrder','quadratic','MesherVersion','R2013a');
data.model=model;
% figure; 
% pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
% xlabel('X-coordinate, meters')
% ylabel('Y-coordinate, meters')
% axis square
for i=1:5
    for j=1:10
        x(i,j)=i*sin(j*2*pi/10)/6;
        y(i,j)=i*cos(j*2*pi/10)/6;
    end
end
D_in=[x(:),y(:)];

data.x_obs=x(:);
data.y_obs=y(:);

obs=simulator(data,Sigma);
sigma_e=0.05;

data.D_in=D_in;
data.obs=obs+sigma_e*randn(1,length(obs));
data.obs_no_err=obs;
data.sigmasq=sigma_e^2;
end
