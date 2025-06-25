clear

% Specify the base folder where the mesh folders are located
baseFolder = 'C:/Users/Alessandro/Documents/GITHUB/Heat_eq_FEM_MATLAB';
saveFolder = 'C:/Users/Alessandro/Documents/GITHUB/Heat_eq_FEM_MATLAB/figure';
t_max = 10;
dt =0.02;
t1 = t_max/4;
t2 = t_max/2;
t3 = t_max*3/4;
t4 = t_max;
teta = 0.5;


%loading ref solution
Ref = load('solRef.dat');
xRef = Ref(:,1);
yRef = Ref(:,2);
uRef = Ref(:,3);
interp = scatteredInterpolant(xRef, yRef, uRef);

iter_chol_tot = zeros(1, 5); 
iter_jac_tot = zeros(1, 5); 
epsilon_Chol = zeros(1, 5); 
epsilon_Jac = zeros(1, 5);  
h = zeros(1, 5);            
r_chol_tot = zeros(1, 5);   
r_jac_tot = zeros(1, 5);

fig = 0;

% Iterate from 0 to 4
for i = 0:5 % Modificato per arrivare fino a 4 come previsto
    %--------------UPLOAD THE MASH FROM THE CORRECRT FOLDER----------------
    % Construct the folder path for the current iteration
    currentFolder = fullfile(baseFolder, ['mesh', num2str(i)]);
    
    coordFile = fullfile(currentFolder, sprintf('mesh%i.coord',  i));
    topolFile = fullfile(currentFolder, sprintf('mesh%i.topol',  i));
    boundFile = fullfile(currentFolder, sprintf('mesh%i.bound',  i));
    traceFile = fullfile(currentFolder, sprintf('mesh%i.trace',  i));
    trackFile = fullfile(currentFolder, sprintf('mesh%i.track',  i));

    % Check if the files exist before attempting to load
    if exist(coordFile, 'file') && exist(topolFile, 'file') && exist(boundFile, 'file') && exist(traceFile, 'file') && exist(trackFile, 'file')
        % Load the mesh files
        coord = load(coordFile);
        topol = load(topolFile);
        bound = load(boundFile);
        trace = load(traceFile);
        track = load(trackFile);
        
        % Display information
        disp(['Loaded mesh data for iteration ', num2str(i)]);
    else
        disp(['Mesh files not found for iteration ', num2str(i), '. Skipping to next iteration.']);
         
    end
     coord = load(coordFile);
     topol = load(topolFile);
     bound = load(boundFile);
     trace = load(traceFile);
     track = load(trackFile);
        

    x=coord(:,1);
    y=coord(:,2);    
    u0 = zeros(size(coord,1),1);
    time_1 = tic;
    %-------------------COMPUTATION OF THE SOLUTIONS-----------------------
    [u_chol, iter_chol, resvec_chol, u_jac, iter_jac, resvec_jac, delta] = my_fem(topol, coord, bound, t_max, dt, teta, u0);
    time_2 = toc(time_1);
    disp(time_2);
    iter_chol_tot(i+1) = iter_chol;
    iter_jac_tot(i+1) = iter_jac;
    
    %-------------------SCATTER PLOT OF THE SOLUTION-----------------------
    fig = fig + 1;
    figure(fig);
    scatter3(x, y, u_chol(:,end), 50, u_chol(:,end), 'filled'); % 50 è la dimensione dei marker, 'filled' li riempie
    colorbar;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('U-value');
    xlim([-1, 1]); 
    ylim([-1, 1]); 
    zlim([0, 1]);
    title('3D plot of the solution at the final time step');
    saveas(gcf, fullfile(saveFolder, sprintf('u_chol_mesh%i.png',  i)));   

    %----------------PLOT OF THE RESIDUAL NORM AT EACH ITERATION-----------
    fig = fig + 1;
    figure(fig)
    semilogy(0:iter_chol,resvec_chol)
    xlabel('iterations');
    ylabel('semilog-residual norm');
    title('plot of the residual norm at each iteration, Cholesky precondition');
    saveas(gcf, fullfile(saveFolder, sprintf('Semilog_Residual_chol_mesh%i.png',  i)));  

    fig = fig + 1;
    figure(fig)
    semilogy(0:iter_jac,resvec_jac)
    xlabel('iterations');
    ylabel('semilog-residual norm');
    title('plot of the residual norm at each iteration, Jacobi precondition');
    saveas(gcf, fullfile(saveFolder, sprintf('Semilog_Residual_jac_mesh%i.png',  i)));  
    
    %-----------------SOLUTION ALONG THE BOUNDARY OF THE MESH--------------
    fig = fig + 1;
    figure(fig)
    plot(trace(:,2),u_chol(trace(:,1),t1/dt+1))
    hold on
    plot(trace(:,2),u_chol(trace(:,1),t2/dt+1))
    plot(trace(:,2),u_chol(trace(:,1),t3/dt+1))
    plot(trace(:,2),u_chol(trace(:,1),t4/dt+1))
    hold off
    legend('2.5','5.0','7.5','10.0');
    xlabel('Arc length');
    ylabel('u');
    xlim([0, 6]); 
    ylim([0, 1]);
    title(sprintf('Solution along the boundary of mesh %i',i))
    saveas(gcf, fullfile(saveFolder, sprintf('trace_u_chol_mesh%i.png',  i)));  

    %-----------------PLOT OF THE SOLUTION ON THE TRACKING POINTS----------
    fig = fig + 1;
    figure(fig)
    plot(0:dt:t_max,u_chol(track(1),:))
    hold on
    plot(0:dt:t_max,u_chol(track(2),:))
    plot(0:dt:t_max,u_chol(track(3),:))
    hold off
    legend('P1','P2','P3');
    xlabel('time');
    ylabel('u');
    xlim([0, 10]); 
    ylim([0, 1]);
    title(sprintf('the three tracking points for all the times of mesh %i',i))
    saveas(gcf, fullfile(saveFolder, sprintf('track_chol_mesh%i.png',  i)));  
    
    %------------------VALUES ON THE TRACKING POINTS AT DIFFERENT TIME---------------------------------   
    P1(1) = u_chol(track(1),t1/dt);
    P1(2) = u_chol(track(1),t2/dt);
    P1(3) = u_chol(track(1),t3/dt);
    P1(4) = u_chol(track(1),t4/dt);

    P2(1) = u_chol(track(2),t1/dt);
    P2(2) = u_chol(track(2),t2/dt);
    P2(3) = u_chol(track(2),t3/dt);
    P2(4) = u_chol(track(2),t4/dt);

    P3(1) = u_chol(track(3),t1/dt);
    P3(2) = u_chol(track(3),t2/dt);
    P3(3) = u_chol(track(3),t3/dt);
    P3(4) = u_chol(track(3),t4/dt);
    
    tabella = table([t1,t2,t3,t4],P1, P2, P3);

    
    disp(tabella);
    %---------------------EPSILON------------------------------------------
    eps_chol = 0;
    eps_jac = 0;

    soli = interp(x, y);
    eps_chol = sqrt(sum((u_chol(:, end) - soli.^2) .* (delta/3)));
    eps_jac = sqrt(sum((u_jac(:, end) - soli.^2) .* (delta/3)));

    epsilon_Chol(i + 1) = eps_chol;
    epsilon_Jac(i + 1) = eps_jac;

    %--------------------------RATIO---------------------------------------
    % Creation of all index pairs
    combination = nchoosek(1:size(topol,2), 2);
    
    % Calculation of Euclidean distances between all pairs of points
    distance = sqrt((x(combination(:, 1)) - x(combination(:, 2))).^2 + (y(combination(:, 1)) - y(combination(:, 2))).^2);
    
    % finding the max distance
    htot = max(distance);
    h(i + 1) = htot;
    r_chol = 0;
    r_jac = 0;
    
    if i>0
        r_chol = (epsilon_Chol(i) ./ epsilon_Chol(i + 1)) .* (h(i + 1) ./ h(i)).^2;
        r_jac = (epsilon_Jac(i) ./ epsilon_Jac(i + 1)) .* (h(i + 1) ./ h(i)).^2;

        r_chol_tot(i + 1) = r_chol;
        r_jac_tot(i + 1) = r_jac;
    end
end



%Display last table
table([0,1,2,3,4], iter_chol_tot, epsilon_Chol, r_chol_tot, iter_jac_tot, epsilon_Jac, r_jac_tot)
