function Model_to_update = Model_Generator(PDT,Model_to_update,Model_updated)

if nargin == 2
    eg = Model_to_update.eg;
    switch eg 
        case 0
            Model_to_update = PDE_set_up_toy_0(PDT,Model_to_update);
        case 1
            Model_to_update = PDE_set_up_toy(PDT,Model_to_update);
        case 2
            Model_to_update = PDE_set_up_coupled_2d(PDT,Model_to_update);
        case 3
            Model_to_update = PDE_set_up_heat(PDT,Model_to_update);
        case 4
            Model_to_update = PDE_set_up_coupled_3d(PDT,Model_to_update);
        case 5
            Model_to_update = PDE_set_up_heat_2d(PDT,Model_to_update);
        case 6
            Model_to_update = PDE_set_up_coupled_3d_2(PDT,Model_to_update);
        case 7
            Model_to_update = PDE_set_up_coupled_3d_3(PDT,Model_to_update);
        case 8
            Model_to_update = PDE_set_up_coupled_2d_2(PDT,Model_to_update);
        case 9
            Model_to_update = PDE_set_up_coupled_3d_2_2(PDT,Model_to_update);
    end
else
    Model_to_update.obs = Model_updated.obs;
    Model_to_update.D_true = Model_updated.D_true;
    Model_to_update.sigmasq = Model_updated.sigmasq;
    Model_to_update.obs_no_err = Model_updated.obs_no_err;
end

end
%%
function pde_model = PDE_set_up_toy_0(PDT,pde_model)
    pde_model.eg = 0;
    gd=[1;0;0;1];
    model = createpde(1);
    g = decsg(gd, 'R', 'R');
    geometryFromEdges(model,g);
    %hmax = 0.1;
    mesh = generateMesh(model,'GeometricOrder','quadratic','MesherVersion','R2013a');
    pde_model.model=model;
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
    pde_model.D_in=D_in;
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.coeff_ind = [8,3,2,4,1,6,5,7,9];
    obs=simulator(pde_model,PDT);
    sigma_e=0.1;
    pde_model.D_true = PDT;
    pde_model.data_structure = [size(D_in,1)]; % a vector, each element is the number of observations of each simulator.
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end

%%
function pde_model = PDE_set_up_toy(PDT,pde_model)
    pde_model.eg = 1;
    gd=[1;0;0;1];
    model = createpde(1);
    g = decsg(gd, 'R', 'R');
    geometryFromEdges(model,g);
    %hmax = 0.1;
    mesh = generateMesh(model,'GeometricOrder','quadratic','MesherVersion','R2013a');
    pde_model.model=model;
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
    pde_model.D_in=D_in;
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.coeff_ind = [8,3,2,4,1,6,5,7,9];
    obs=simulator(pde_model,PDT);
    sigma_e=0.05;
    pde_model.D_true = PDT;
    pde_model.data_structure = [size(D_in,1)]; % a vector, each element is the number of observations of each simulator.
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end
%%
function pde_model=PDE_set_up_coupled_2d(PDT,pde_model)
    pde_model.eg = 2;
    % consider a simple example, which is a linear PDE system 2d
    N = 2;
    model = createpde(N);
    L = 0.3; % beam length in meters
    H = 0.05; % overall height of the beam
    R = [3 4 0 L L 0 0 0 H H];
    gdm = R';
    g = decsg(gdm, 'R', 'R');
    geometryFromEdges(model,g);
    % figure; 
    % pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
    % xlabel('X-coordinate, meters')
    % ylabel('Y-coordinate, meters')
    % axis([-0.1*L, 1.1*L, -H, 2*H])
    % axis square
    %hmax = 0.03;
    %mesh = generateMesh(model,'Hmax',hmax,'GeometricOrder','quadratic','MesherVersion','R2013a');
    mesh = generateMesh(model);
    %pdeplot(model) 
    pde_model.model=model;
    e=L/5;
    eh = H/5;
    [x,y]=meshgrid(e:e:L-e,eh:eh:H-eh);
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.t_list=0:100:200000;
    pde_model.t_obs=20:400:2000;
    pde_model.t_obs2=20:400:2000;
    D_in=repmat([pde_model.x_obs,pde_model.y_obs],length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs,length(pde_model.x_obs(:)),1);
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure = [size(D_in,1),size(D_in,1)];
    pde_model.D_in=[pde_model.D_in;pde_model.D_in];
    [obs,obs_ini]=simulator(pde_model,PDT);
    sigma_e=0.01;
    pde_model.D_true = PDT;
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_no_err=obs;
    pde_model.obs_ini=obs_ini;
    pde_model.sigmasq=sigma_e^2;
end

%% eg 3
function pde_model=PDE_set_up_heat(PDT,pde_model)
    pde_model.eg = 3;
    N = 1;
    model = createpde(N);
    gm = multicuboid(1,1,1);
    model.Geometry = gm;
    %figure; 
    %pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
    %axis square
    hmax = 0.2;
    %mesh = generateMesh(model,'Hmax',hmax);
    mesh = generateMesh(model,'Hmax',hmax);
    pde_model.model=model;
    x=(0.2:0.3:0.8)-0.5;
    y=(0.2:0.3:0.8)-0.5;
    z=0.2:0.3:0.8;
    [x,y,z]=meshgrid(x,y,z);
    t_list=(0:500)/500;
    t_obs=[101,501];
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.z_obs=z(:);
    pde_model.t_obs=t_obs(:);
    pde_model.t_list=t_list;
    obs=simulator(pde_model,PDT);
    sigma_e=0.001;
    pde_model.D_true = PDT;
    D_in=repmat([pde_model.x_obs,pde_model.y_obs,pde_model.z_obs],length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs/500,1,length(pde_model.x_obs(:)))';
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure=[size(D_in,1)];
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end


%%
function pde_model=PDE_set_up_coupled_3d(PDT,pde_model)
    pde_model.eg = 4;
    N = 2;
    W1=0.048;
    W2=0.048;
    W3=0.02;
    model = createpde(N);
    gm = multicuboid(W1,W2,W3);
    model.Geometry = gm;
    hmax = 0.005;
    %mesh = generateMesh(model,'Hmax',hmax);
    mesh = generateMesh(model);
    pde_model.model=model;
    x=-W1/3:W1/3:W1/3;
    y=-W2/3:W2/3:W2/3;
    z=W3/6:W3/3:5*W3/6;
    [x,y,z]=meshgrid(x,y,z);
    t_list=0:4:2000;
    t_obs=[200,500];
    D_in=[x(:),y(:),z(:)];
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.z_obs=z(:);
    pde_model.t_obs=t_obs(:);
    pde_model.t_list=t_list;
    D_in=repmat([pde_model.x_obs,pde_model.y_obs],length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs,length(pde_model.x_obs(:)),1);
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure = [size(D_in,1),size(D_in,1)];
    pde_model.D_in=[pde_model.D_in;pde_model.D_in];
    [obs,obs_ini]=simulator(pde_model,PDT);
    sigma_e=0.01;
    
    pde_model.D_true = PDT;
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_ini=obs_ini;
    
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end
%%
function pde_model=PDE_set_up_heat_2d(PDT,pde_model)
    pde_model.eg = 5;
    N = 1;
    model = createpde(N);
    L = 0.0366; % beam length in meters
    H =L; % overall height of the beam
    R = [3 4 0 L L 0 0 0 H H];
    gdm = R';
    g = decsg(gdm, 'R', 'R');
    geometryFromEdges(model,g);
    mesh = generateMesh(model);
    %figure; 
    %pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
    %axis square
    pde_model.model=model;
    x=(0.2:0.3:0.8);
    y=(0.2:0.3:0.8);
    [x,y]=meshgrid(x*L,y*L);
    t_list=0:100;
    t_obs=[50,100];
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.t_obs=t_obs(:);
    pde_model.t_list=t_list;
    obs=simulator(pde_model,PDT);
    sigma_e=0.1;
    pde_model.D_true = PDT;
    D_in=repmat([pde_model.x_obs,pde_model.y_obs],length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs,1,length(pde_model.x_obs(:)))';
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure=[size(D_in,1)];
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end

%%
function pde_model=PDE_set_up_coupled_3d_2(PDT,pde_model)
    pde_model.eg = 6;
    N = 2;
    W1=2;
    W2=0.2;
    W3=0.05;
    model = createpde(N);
    gm = multicuboid(W1,W2,W3);
    model.Geometry = gm;
    %hmax = 0.005;
    %mesh = generateMesh(model,'Hmax',hmax);
    mesh = generateMesh(model);
    pde_model.model=model;
    x=-W1/3:W1/3:W1/3;
    y=-W2/3:W2/3:W2/3;
    z=W3/6:W3/3:5*W3/6;
    [x,y,z]=meshgrid(x,y,z);
    t_list=0:100:20000;
    t_obs=[100,200];
    D_in=[x(:),y(:),z(:)];
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.z_obs=z(:);
    pde_model.t_obs=t_obs(:);
    pde_model.t_list=t_list;
    D_in=repmat(D_in,length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs,length(pde_model.x_obs(:)),1);
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure = [size(D_in,1),size(D_in,1)];
    pde_model.D_in=[pde_model.D_in;pde_model.D_in];
    [obs,obs_ini]=simulator(pde_model,PDT);
    sigma_e=0.01;
    
    pde_model.D_true = PDT;
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_ini=obs_ini;
    
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end
%%
function pde_model=PDE_set_up_coupled_3d_3(PDT,pde_model)
    pde_model.eg = 7;

    N = 2;
    model = createpde(N);
    gm = multicuboid(2,0.2,0.05);
    model.Geometry = gm;

    %mesh = generateMesh(model,'Hmax',hmax);
    mesh = generateMesh(model);
    W1=2;
    W2=0.2;
    W3=0.05;
    pde_model.model=model;
    x=-W1/3:W1/3:W1/3;
    y=-W2/3:W2/3:W2/3;
    z=W3/6:W3/3:5*W3/6;
    [x,y,z]=meshgrid(x,y,z);
    t_list=0:100:20000;
    t_obs=[100,200];
    D_in=[x(:),y(:),z(:)];
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.z_obs=z(:);
    pde_model.t_obs=t_obs(:);
    pde_model.t_list=t_list;
    D_in=repmat(D_in,length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs,length(pde_model.x_obs(:)),1);
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure = [size(D_in,1),size(D_in,1)];
    pde_model.D_in=[pde_model.D_in;pde_model.D_in];
    [obs,obs_ini]=simulator(pde_model,PDT);
    sigma_e=0.01;
    
    pde_model.D_true = PDT;
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_ini=obs_ini;
    
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end
%%
function pde_model=PDE_set_up_coupled_2d_2(PDT,pde_model)
    pde_model.eg = 8;
    % consider a simple example, which is a linear PDE system 2d
    N = 2;
    model = createpde(N);
    L = 0.3; % beam length in meters
    H = 0.05; % overall height of the beam
    R = [3 4 0 L L 0 0 0 H H];
    gdm = R';
    g = decsg(gdm, 'R', 'R');
    geometryFromEdges(model,g);
    % figure; 
    % pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
    % xlabel('X-coordinate, meters')
    % ylabel('Y-coordinate, meters')
    % axis([-0.1*L, 1.1*L, -H, 2*H])
    % axis square
    %hmax = 0.03;
    %mesh = generateMesh(model,'Hmax',hmax,'GeometricOrder','quadratic','MesherVersion','R2013a');
    mesh = generateMesh(model);
    %pdeplot(model) 
    pde_model.model=model;
    e=L/5;
    eh = H/5;
    [x,y]=meshgrid(e:e:L-e,eh:eh:H-eh);
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.t_list=0:100:200000;
    pde_model.t_obs=20:400:2000;
    pde_model.t_obs2=20:400:2000;
    D_in=repmat([pde_model.x_obs,pde_model.y_obs],length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs,length(pde_model.x_obs(:)),1);
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure = [size(D_in,1)];
    pde_model.D_in=pde_model.D_in;
    [obs,obs_ini]=simulator(pde_model,PDT);
    sigma_e=0.01;
    pde_model.D_true = PDT;
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_no_err=obs;
    pde_model.obs_ini=obs_ini;
    pde_model.sigmasq=sigma_e^2;
end
%%
function pde_model=PDE_set_up_coupled_3d_2_2(PDT,pde_model)
    pde_model.eg = 9;
    N = 2;
    W1=2;
    W2=0.2;
    W3=0.05;
    model = createpde(N);
    gm = multicuboid(W1,W2,W3);
    model.Geometry = gm;
    %hmax = 0.005;
    %mesh = generateMesh(model,'Hmax',hmax);
    mesh = generateMesh(model);
    pde_model.model=model;
    x=-W1/3:W1/3:W1/3;
    y=-W2/3:W2/3:W2/3;
    z=W3/6:W3/3:5*W3/6;
    [x,y,z]=meshgrid(x,y,z);
    t_list=0:100:20000;
    t_obs=[100,200];
    D_in=[x(:),y(:),z(:)];
    pde_model.x_obs=x(:);
    pde_model.y_obs=y(:);
    pde_model.z_obs=z(:);
    pde_model.t_obs=t_obs(:);
    pde_model.t_list=t_list;
    D_in=repmat(D_in,length(pde_model.t_obs),1);
    D_t=repmat(pde_model.t_obs,length(pde_model.x_obs(:)),1);
    D_t=D_t(:);
    pde_model.D_in=[D_in,D_t];
    pde_model.data_structure = [size(D_in,1)];
    [obs,obs_ini]=simulator(pde_model,PDT);
    sigma_e=0.01;
    
    pde_model.D_true = PDT;
    pde_model.obs=obs+sigma_e*randn(1,length(obs));
    pde_model.obs_ini=obs_ini;
    
    pde_model.obs_no_err=obs;
    pde_model.sigmasq=sigma_e^2;
end

