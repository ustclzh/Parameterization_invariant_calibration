function [y,obs_ini]=simulator(data,D)
% all simulators are here. 
% 1. toy example
% 2. coupled PDE 2d
% 3. 3-d heat equation
% 4. coupled PDE 3-d
% 5. 2-d heat equation
% 6. coupled PDE 3-d second
% 7. coupled PDE 3-d third 
% 8. coupled PDE 2-d second
% 9. coupled PDE 3-d fourth

if length(size(D)) == 1
    D = invparametrization(D,data.parameterization);
end
switch data.eg
    case 0
        y = toy_0(data,D);
        obs_ini = [];
    case 1
        y = toy(data,D);
        obs_ini = [];
    case 2
        [y,obs_ini] = coupled_2d(data,D);
    case 3
        y = eg_heat(data,D);
        obs_ini = [];
    case 4
        [y,obs_ini] = coupled_3d(data,D);
    case 5 
        y = eg_heat_2d(data,D);
        obs_ini = [];
    case 6
        [y,obs_ini] = coupled_3d_2(data,D);   
    case 7
        [y,obs_ini] = coupled_3d_3(data,D);  
    case 8 
        [y,obs_ini] = coupled_2d_2(data,D);
    case 9 
        [y,obs_ini] = coupled_3d_2_2(data,D);   
end
end

%%
function y = toy_0(data,D)
x = [data.x_obs(:),data.y_obs(:)];
y = diag(x * D * x');
y=y(:)';
end

%%
function y = toy(data,D)

c=D(:);
model=settingforpde();%construct the PDE geometry
specifyCoefficients(model, 'm', 0,'d',0,'c', c, 'a', 0, 'f', 0);
applyBoundaryCondition(model,'neumann','Edge',[1,3,4],'q',0,'g',10);
applyBoundaryCondition(model,'dirichlet','Edge',[2],'h',1,'r',0);
result = solvepde(model);
y=interpolateSolution(result,data.x_obs(:),data.y_obs(:));
y=y(:)';
end
function fmatrix=fcoef_toy(location, state)
n1 = 1;
nr = numel(location.x);
fmatrix = zeros(n1,nr);
fmatrix(1,:) = 1000*((location.x.^2+location.y.^2)<0.1).*ones(1,nr);

end
function model = settingforpde()
gd=[1;0;0;1];
model = createpde(1);
g = decsg(gd, 'R', 'R');
geometryFromEdges(model,g);
%hmax = 0.1;
generateMesh(model,'GeometricOrder','quadratic','MesherVersion','R2013a');
end
%%
function [y,obs_ini]=coupled_2d(data,in)
% Coupled moisture and heat transfer equation 2d
model=data.model;
% parameters 
C=2;
M0=0.8;
phi=0.3;%relative air humidity
rho=0.4*100^3;% wood density
T_air=350;% tempreture of air
T0= 300;
A=7.731706-0.014348*T_air;
B=0.008746+0.000567*T_air;
RT=8.31*T_air;% R is gas constant
E_b=1;%water diffusion 
alpha_M=8*10^(-7);%moisture transfer coefficient
alpha_T=1;%heat transfer coefficient
s=E_b/(B*log(1/phi)*T_air*RT);%Soret effect 
EMC=1/B*log(A/log(1/phi))/100;%equilibrium moisture content from paper
T_F_air=(T_air-273.15)*1.8+32;
W=330+0.452*T_F_air+0.00415*T_F_air^2;
k=0.791+4.63*10^(-4)*T_F_air-8.44*10^(-7)*T_F_air^2;
k1=6.34+7.75*10^(-4)*T_F_air-9.35*10^(-5)*T_F_air^2;
k2=1.09+2.84*10^(-2)*T_F_air-9.04*10^(-5)*T_F_air^2;
EMC=1800/W*(k*phi/(1-k*phi)+(k1*k*phi+2*k1*k2*k^2*phi^2)/(1+k1*k*phi+k1*k2*k^2*phi^2))/100;%equilibrium moisture content from wikipedia
D_M=in*10^(-9);%[55*10^(-9),0;0,5*10^(-9)];%matrix of moisture diffusion coefficients of wood
D_T=eye(2)*10^(-2);%matrix of thermal conductivity coefficients of wood 
c = [D_M(:);zeros(4,1);s*D_M(:);D_T(:)];
d=[1,0;-E_b/(1.8*C),rho*C];
q=diag([alpha_M,alpha_T]');
g=[alpha_M*EMC,alpha_T*T_air]';
u0=[M0;T0];
specifyCoefficients(model, 'm', 0,'d',d(:),'c', c, 'a', 0, 'f',[0;0]);
applyBoundaryCondition(model,'neumann','Edge',[1,2,3,4],'q',q,'g',g);
setInitialConditions(model,u0);
t_list=data.t_list;
%tic
result = solvepde(model,t_list);
%toc
obs_ini1=100*interpolateSolution(result,data.x_obs(:),data.y_obs(:),1,data.t_obs);
y=obs_ini1(:)';
obs_ini2=interpolateSolution(result,data.x_obs(:),data.y_obs(:),2,data.t_obs2);
y=[y,obs_ini2(:)'];
obs_ini = [obs_ini1,obs_ini2];
end
function fmatrix=fcoef_coupled_2d(location, state)
n1 = 2;
nr = numel(location.x);
fmatrix = zeros(n1,nr);
fmatrix(2,:) = (location.x>=0).*(location.x<0.3).*(location.y>0.025).*(location.y<0.035)*0.1.*ones(1,nr);
end

%% eg 3
function y=eg_heat(data,D)
% heat equation
c=D;
c=c(:);
model=data.model;
specifyCoefficients(model, 'm', 0,'d',10,'c', c, 'a', 0, 'f', 50);
%applyBoundaryCondition(model,'neumann','Face',[1,2],'q',0,'g',1);
applyBoundaryCondition(model,'dirichlet','Face',[1,2,3,5,4,6],'h',1,'r',0);
setInitialConditions(model,0);
results = solvepde(model,data.t_list);
obs_ini=interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),data.t_obs);
y=obs_ini(:)';
end

%%
function [y,obs_ini]=coupled_3d(data,D)
% Coupled moisture and heat transfer equation
model=data.model;
delta=2;
epsilon=0.3;
rho=856;
lambda=2.5*10^(6);
kq=0.577*eye(3);
km=2.2*10^(-8)*D;
cm=0.01;
cq=4201.4;
hq=25;
hm=0.0001;
T0=10;
Ta=110;
U0=87;
Ua=12;
q=[hq,lambda*hm;0,hm];%[hq,(1-epsilon)*lambda*hm;0,hm];%
g=q*[Ta;Ua];%[hq*Ta+lambda*hm*Ua;hm*Ua];%+lambda*hm*Ua
u0=[T0;U0];
dcoef = @(location, state) dcoef_younsi(location, state, cq, cm, rho, epsilon,lambda);
ccoef = @(location, state) ccoef_younsi(location, state, kq, km ,cm, lambda, epsilon, delta);
specifyCoefficients(model, 'm', 0,'d',dcoef,'c', ccoef, 'a', 0, 'f', [0;0]);
applyBoundaryCondition(model,'neumann','Face',[1:6],'q',-q,'g',-g);
setInitialConditions(model,u0);
%tic
results = solvepde(model,data.t_list);
%toc
obs_ini=interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),1,data.t_obs);
obs_ini_t=interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),2,data.t_obs);
y=[obs_ini(:)', obs_ini_t(:)'];%,obs_ini_t(:)'];
end

function dmatrix=dcoef_younsi(location, state, cq, cm, rho, epsilon,lambda)
n1 = 4;
nr = numel(location.x);
dmatrix = zeros(n1,nr);
dmatrix(1,:) = rho*cq*ones(1,nr);
%dmatrix(3,:) = -rho*cm*epsilon*lambda*ones(1,nr);
dmatrix(4,:) = rho*cm*ones(1,nr);
end

function cmatrix=ccoef_younsi(location, state, kq, km ,cm, lambda, epsilon, delta)
% this is conductivity/permeability/diffucivity term
% for compatibility, we always use this functional form instead of matrix
% form.\
%c = [D_M(:);zeros(4,1);s*D_M(:);D_T(:)];
nr = numel(location.x);
D_q=kq+epsilon*lambda*delta*km/cm;
D_q=D_q(:);
D_qm=delta*km/cm;
D_qm=D_qm(:);
D_mq=epsilon*lambda*km;
D_mq=D_mq(:);
D_m=km;
D_m = D_m(:);
cmatrix= [D_q*ones(1,nr);D_qm*ones(1,nr);D_mq*ones(1,nr);D_m*ones(1,nr)];
end
%% 5 eg

function y=eg_heat_2d(data,D)
% heat equation
c=D;
c=c(:);
model=data.model;
specifyCoefficients(model, 'm', 0,'d',1465550,'c', c, 'a', 0, 'f', 0);
applyBoundaryCondition(model,'neumann','Edge',[1,2],'q',0,'g',396.99);
applyBoundaryCondition(model,'dirichlet','Edge',[3,4],'h',1,'r',0);
setInitialConditions(model,20);
results = solvepde(model,data.t_list);
obs_ini=interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.t_obs);
y=obs_ini(:)';
end

%% 6 eg

function [y,obs_ini]=coupled_3d_2(data,D)
% Coupled moisture and heat transfer equation

model=data.model;
C=2;

phi=0.5;%relative air humidity
rho=4000;% wood density
M0=0.25;
T0 = 300;
T_air=353;% tempreture of air

A=7.731706-0.014348*T_air;
B=0.008746+0.000567*T_air;
RT=8.31*T_air;% R is gas constant
E_b=1;%water diffusion 
s=E_b/(B*log(1/phi)*T_air*RT);%Soret effect 

EMC=1/B*log(A/log(1/phi))/100;%equilibrium moisture content from paper
T_F_air=(T_air-273.15)*1.8+32;
W=330+0.452*T_F_air+0.00415*T_F_air^2;
k=0.791+4.63*10^(-4)*T_F_air-8.44*10^(-7)*T_F_air^2;
k1=6.34+7.75*10^(-4)*T_F_air-9.35*10^(-5)*T_F_air^2;
k2=1.09+2.84*10^(-2)*T_F_air-9.04*10^(-5)*T_F_air^2;
%EMC=1800/W*(k*phi/(1-k*phi)+(k1*k*phi+2*k1*k2*k^2*phi^2)/(1+k1*k*phi+k1*k2*k^2*phi^2))/100;%equilibrium moisture content from wikipedia

D_M=10^(-8)*D;

D_T=diag([0.52,0.2,0.31]);%matrix of thermal conductivity coefficients of wood 

alpha_M=5*10^(-7);%moisture transfer coefficient
alpha_T=20;%heat transfer coefficient
q=diag([alpha_M,alpha_T]');
g=[alpha_M*EMC,alpha_T*T_air]';
u0=[M0;T0];
%H=1;
%fcoef = @(location, state) fcoef_coupled_3d(location, state, H); 
dcoef = @(location, state) dcoef_coupled_3d(location, state, C, E_b, rho);
ccoef = @(location, state) ccoef_coupled_3d(location, state, D_M, D_T, s);
specifyCoefficients(model, 'm', 0,'d',dcoef,'c', ccoef, 'a', 0, 'f', [0;0]);
applyBoundaryCondition(model,'neumann','Face',[1:6],'q',q,'g',g);
setInitialConditions(model,u0);
%tic
results = solvepde(model,data.t_list);
%toc
obs_ini=interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),1,data.t_obs);
obs_ini_t=(interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),2,data.t_obs)-300);
y=[100*obs_ini(:)', obs_ini_t(:)'];%,obs_ini_t(:)'];
end

function fmatrix=fcoef_coupled_3d(location, state,H)%source term
    % this is the source term
    n1 = 2;
    nr = numel(location.x);
    fmatrix = zeros(n1,nr);
    fmatrix(2,:) = (location.x>=0).*(location.x<0.3).*(location.y>(H/2-0.005)).*(location.y<(H/2+0.005))*0.1.*ones(1,nr);
end
function dmatrix=dcoef_coupled_3d(location, state, C, E_b, rho)
    %d=[1,0;-E_b/(1.8*C),rho*C];
    n1 = 4;
    nr = numel(location.x);
    dmatrix = zeros(n1,nr);
    dmatrix(1,:) = ones(1,nr);
    dmatrix(2,:) = -E_b/(1.8*C)*ones(1,nr);
    dmatrix(4,:) = rho*C*ones(1,nr);
end
function cmatrix=ccoef_coupled_3d(location, state, D_M, D_T, s)%tensor permeability
    % this is conductivity/permeability/diffucivity term
    % for compatibility, we always use this functional form instead of matrix
    % form.\
    %c = [D_M(:);zeros(4,1);s*D_M(:);D_T(:)];
    nr = numel(location.x);
    cmatrix=[D_M(:)*ones(1,nr);zeros(9,1)*ones(1,nr);s*D_M(:)*ones(1,nr);D_T(:)*ones(1,nr);];
end


%% 7 eg
function [y,obs_ini]=coupled_3d_3(data,D)
% Coupled moisture and heat transfer equation
model=data.model;
rho=420;
alpha_M=5*10^(-7);
alpha_T=20;
C=1535;
s=0.00222;
T_air=353;
EMC=0.06;
M0=0.25;
T0=293;
E_b=1;%water diffusion 
%EMC=1800/W*(k*phi/(1-k*phi)+(k1*k*phi+2*k1*k2*k^2*phi^2)/(1+k1*k*phi+k1*k2*k^2*phi^2))/100;%equilibrium moisture content from wikipedia
D_M=0^(-8)*D;
D_T=diag([0.52,0.2,0.31]);%matrix of thermal conductivity coefficients of wood 
q=diag([alpha_M,alpha_T]');
g=[alpha_M*EMC,alpha_T*T_air]';
u0=[M0;T0];
%H=1;
%fcoef = @(location, state) fcoef_coupled_3d(location, state, H); 
dcoef = @(location, state) dcoef_coupled_3d_3(location, state, C, E_b, rho);
ccoef = @(location, state) ccoef_coupled_3d_3(location, state, D_M, D_T, s);
specifyCoefficients(model, 'm', 0,'d',dcoef,'c', ccoef, 'a', 0, 'f', [0;0]);
applyBoundaryCondition(model,'neumann','Face',[1:6],'q',q,'g',g);
setInitialConditions(model,u0);
%tic
results = solvepde(model,data.t_list);
%toc
obs_ini=interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),1,data.t_obs);
obs_ini_t=(interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),2,data.t_obs)-300);
y=[100*obs_ini(:)', obs_ini_t(:)'];%,obs_ini_t(:)'];
end

function fmatrix=fcoef_coupled_3d_3(location, state,H)%source term
    % this is the source term
    n1 = 2;
    nr = numel(location.x);
    fmatrix = zeros(n1,nr);
    fmatrix(2,:) = (location.x>=0).*(location.x<0.3).*(location.y>(H/2-0.005)).*(location.y<(H/2+0.005))*0.1.*ones(1,nr);
end
function dmatrix=dcoef_coupled_3d_3(location, state, C, E_b, rho)
    %d=[1,0;-E_b/(1.8*C),rho*C];
    n1 = 4;
    nr = numel(location.x);
    dmatrix = zeros(n1,nr);
    dmatrix(1,:) = ones(1,nr);
    dmatrix(2,:) = -E_b/(1.8*C)*ones(1,nr);
    dmatrix(4,:) = rho*C*ones(1,nr);
end
function cmatrix=ccoef_coupled_3d_3(location, state, D_M, D_T, s)%tensor permeability
    % this is conductivity/permeability/diffucivity term
    % for compatibility, we always use this functional form instead of matrix
    % form.\
    %c = [D_M(:);zeros(4,1);s*D_M(:);D_T(:)];
    nr = numel(location.x);
    cmatrix=[D_M(:)*ones(1,nr);zeros(9,1)*ones(1,nr);s*D_M(:)*ones(1,nr);D_T(:)*ones(1,nr);];
end
%%
function [y,obs_ini]=coupled_2d_2(data,in)
% Coupled moisture and heat transfer equation 2d
model=data.model;
% parameters 
C=2;
M0=0.8;
phi=0.3;%relative air humidity
rho=0.4*100^3;% wood density
T_air=350;% tempreture of air
T0= 300;
A=7.731706-0.014348*T_air;
B=0.008746+0.000567*T_air;
RT=8.31*T_air;% R is gas constant
E_b=1;%water diffusion 
alpha_M=8*10^(-7);%moisture transfer coefficient
alpha_T=1;%heat transfer coefficient
s=E_b/(B*log(1/phi)*T_air*RT);%Soret effect 
EMC=1/B*log(A/log(1/phi))/100;%equilibrium moisture content from paper
T_F_air=(T_air-273.15)*1.8+32;
W=330+0.452*T_F_air+0.00415*T_F_air^2;
k=0.791+4.63*10^(-4)*T_F_air-8.44*10^(-7)*T_F_air^2;
k1=6.34+7.75*10^(-4)*T_F_air-9.35*10^(-5)*T_F_air^2;
k2=1.09+2.84*10^(-2)*T_F_air-9.04*10^(-5)*T_F_air^2;
EMC=1800/W*(k*phi/(1-k*phi)+(k1*k*phi+2*k1*k2*k^2*phi^2)/(1+k1*k*phi+k1*k2*k^2*phi^2))/100;%equilibrium moisture content from wikipedia
D_M=in*10^(-9);%[55*10^(-9),0;0,5*10^(-9)];%matrix of moisture diffusion coefficients of wood
D_T=eye(2)*10^(-2);%matrix of thermal conductivity coefficients of wood 
c = [D_M(:);zeros(4,1);s*D_M(:);D_T(:)];
d=[1,0;-E_b/(1.8*C),rho*C];
q=diag([alpha_M,alpha_T]');
g=[alpha_M*EMC,alpha_T*T_air]';
u0=[M0;T0];
specifyCoefficients(model, 'm', 0,'d',d(:),'c', c, 'a', 0, 'f',[0;0]);
applyBoundaryCondition(model,'neumann','Edge',[1,2,3,4],'q',q,'g',g);
setInitialConditions(model,u0);
t_list=data.t_list;
%tic
result = solvepde(model,t_list);
%toc
obs_ini1=100*interpolateSolution(result,data.x_obs(:),data.y_obs(:),1,data.t_obs);
y=obs_ini1(:)';
obs_ini2=interpolateSolution(result,data.x_obs(:),data.y_obs(:),2,data.t_obs2);
obs_ini = [obs_ini1,obs_ini2];
end
%% 9 eg
function [y,obs_ini]=coupled_3d_2_2(data,D)
% Coupled moisture and heat transfer equation
model=data.model;
C=2;
phi=0.5;%relative air humidity
rho=4000;% wood density
M0=0.25;
T0 = 300;
T_air=353;% tempreture of air
A=7.731706-0.014348*T_air;
B=0.008746+0.000567*T_air;
RT=8.31*T_air;% R is gas constant
E_b=1;%water diffusion 
s=E_b/(B*log(1/phi)*T_air*RT);%Soret effect 
EMC=1/B*log(A/log(1/phi))/100;%equilibrium moisture content from paper
T_F_air=(T_air-273.15)*1.8+32;
W=330+0.452*T_F_air+0.00415*T_F_air^2;
k=0.791+4.63*10^(-4)*T_F_air-8.44*10^(-7)*T_F_air^2;
k1=6.34+7.75*10^(-4)*T_F_air-9.35*10^(-5)*T_F_air^2;
k2=1.09+2.84*10^(-2)*T_F_air-9.04*10^(-5)*T_F_air^2;
%EMC=1800/W*(k*phi/(1-k*phi)+(k1*k*phi+2*k1*k2*k^2*phi^2)/(1+k1*k*phi+k1*k2*k^2*phi^2))/100;%equilibrium moisture content from wikipedia
D_M=10^(-8)*D;
D_T=diag([0.52,0.2,0.31]);%matrix of thermal conductivity coefficients of wood 
alpha_M=5*10^(-7);%moisture transfer coefficient
alpha_T=20;%heat transfer coefficient
q=diag([alpha_M,alpha_T]');
g=[alpha_M*EMC,alpha_T*T_air]';
u0=[M0;T0];
%H=1;
%fcoef = @(location, state) fcoef_coupled_3d(location, state, H); 
dcoef = @(location, state) dcoef_coupled_3d_2(location, state, C, E_b, rho);
ccoef = @(location, state) ccoef_coupled_3d_2(location, state, D_M, D_T, s);
specifyCoefficients(model, 'm', 0,'d',dcoef,'c', ccoef, 'a', 0, 'f', [0;0]);
applyBoundaryCondition(model,'neumann','Face',[1:6],'q',q,'g',g);
setInitialConditions(model,u0);
%tic
results = solvepde(model,data.t_list);
%toc
obs_ini=interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),1,data.t_obs);
obs_ini_t=(interpolateSolution(results,data.x_obs(:),data.y_obs(:),data.z_obs(:),2,data.t_obs)-300);
y=100*obs_ini(:)';%,obs_ini_t(:)'];
end

function fmatrix=fcoef_coupled_3d_2(location, state,H)%source term
    % this is the source term
    n1 = 2;
    nr = numel(location.x);
    fmatrix = zeros(n1,nr);
    fmatrix(2,:) = (location.x>=0).*(location.x<0.3).*(location.y>(H/2-0.005)).*(location.y<(H/2+0.005))*0.1.*ones(1,nr);
end
function dmatrix=dcoef_coupled_3d_2(location, state, C, E_b, rho)
    %d=[1,0;-E_b/(1.8*C),rho*C];
    n1 = 4;
    nr = numel(location.x);
    dmatrix = zeros(n1,nr);
    dmatrix(1,:) = ones(1,nr);
    dmatrix(2,:) = -E_b/(1.8*C)*ones(1,nr);
    dmatrix(4,:) = rho*C*ones(1,nr);
end
function cmatrix=ccoef_coupled_3d_2(location, state, D_M, D_T, s)%tensor permeability
    % this is conductivity/permeability/diffucivity term
    % for compatibility, we always use this functional form instead of matrix
    % form.\
    %c = [D_M(:);zeros(4,1);s*D_M(:);D_T(:)];
    nr = numel(location.x);
    cmatrix=[D_M(:)*ones(1,nr);zeros(9,1)*ones(1,nr);s*D_M(:)*ones(1,nr);D_T(:)*ones(1,nr);];
end


