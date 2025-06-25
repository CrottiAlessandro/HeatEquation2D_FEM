% FEM function to solve a time-dependent linear partial differential equation
function [u_chol_tot, iter_chol, resvec_chol, u_jac_tot, iter_jac, resvec_jac, my_delta] = my_fem(topol, coord, bound, t_max, dt, teta, u0)


    % Extracting node coordinates
    x = coord(:,1);
    y = coord(:,2);

    % Sizes of topol and coord matrices
    [dim1,dim2] = size(topol);
    [dim3,dim4] = size(coord);

    % Number of time steps
    N=t_max/dt;

    % Initialization of solution matrices
    u_jac_tot = zeros(dim3, N+1);
    u_chol_tot = zeros(dim3,N+1);
    u_jac_tot(:,1) = u0;
    u_chol_tot(:,1) = u0;

    %implentation of a stiffness matrix
    [H, M, delta] = stiff_mat(topol, coord, sparse(dim3,dim3), sparse(dim3,dim3));
    my_delta = zeros(dim3,1);    

    %computation of delta
    for i=1:dim3
        
        for j=1:dim3
            if (topol(j,1)==i || topol(j,2)==i || topol(j,3)==i)
                my_delta(i)=my_delta(i)+delta(j);
            end
        end
    end

    % Coefficients and matrices for the time-stepping scheme
    coeff_N = M/dt;
    coeff_teta1 = teta*H;
    coeff_teta2 = (1-teta)*H;
    K1 = coeff_N + coeff_teta1;
    K2 = coeff_N - coeff_teta2;

    tmp_bound_k1=K1(:,bound(:,1)); 
    
    for i=1:N

        % Right-hand side vectors
        f_chol = K2 * u_chol_tot(:,i);
        f_jac = K2 * u_jac_tot(:,i);
    
    
        % Imposing boundary conditions
        K1(bound(:,1),:) = 0;
        K1(:,bound(:,1)) = tmp_bound_k1; 
        f_chol = f_chol - K1(:,bound(:,1)) * bound(:,2); 
        f_chol(bound(:,1)) = bound(:,2);
        f_jac = f_jac - K1(:,bound(:,1)) * bound(:,2); 
        f_jac(bound(:,1)) = bound(:,2);
        K1(:,bound(:,1)) = 0; 

        % Applying Dirichlet boundary conditions
        nbound=size(bound,1);
        for j = 1:nbound
            K1(bound(j,1),bound(j,1)) = 1; 
        end

       
        % Solving the linear system using Preconditioned Conjugate Gradient
        L = ichol(K1);
        tol = 1e-8;
        Jac=zeros(dim3,dim3);
        for j=1:dim3
             Jac(j,j)=1/K1(j,j);
        end

        % Solving the linear system
        [u_chol,~,~,iter_chol,resvec_chol] = pcg(K1,f_chol,tol,5000,L,L');
        [u_jac,~,~,iter_jac,resvec_jac] = pcg(K1,f_jac,tol,5000,Jac);

        % Updating the solution matrices
        u_chol_tot(:,i+1) = u_chol;
        u_jac_tot(:,i+1) = u_jac;
    end

end

    
