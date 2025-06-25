
% STIFFNESS MATRIX
function [H,M,delta] = stiff_mat(topol, coord, H, M)
    
    % Coefficient matrix for mass matrix
    coef_M = (ones(3,3)+diag(ones(1, 3)));

    % Dimensions of topol matrix
    [dim1, dim2] = size(topol);

    % Extracting node coordinates
    x = coord(:,1);
    y = coord(:,2);

    % Matrices for coefficients
    b_mat = zeros(3, 3);
    c_mat = zeros(3, 3);    

    for z=1:dim1

        i=topol(z,1);
        j=topol(z,2);
        m=topol(z,3);
        
        % Element surface measure
        delta(z) = 0.5*det([1 x(i) y(i); 1 x(j) y(j); 1 x(m) y(m)]);
        
        % Coefficients for linear shape functions
        b = [y(j) - y(m); y(m) - y(i); y(i) - y(j)];
        c = [x(m) - x(j); x(i) - x(m); x(j) - x(i)];


        % Computation of the local stiffness matrix for a triangular element
        b_mat = b * b';
        c_mat = c * c';

        % Local stiffness matrix
        Hloc = (1/(4*delta(z)))*(b_mat+c_mat);        
        Mloc = (delta(z)/12)*coef_M;

        % Assembling the global stiffness and mass matrices
        for i = 1:3
            row = topol(z,i);
            for j = 1:3
                col = topol(z,j);
                H(row,col) = H(row,col) + Hloc(i,j);
                M(row,col) = M(row,col) + Mloc(i,j);
            end 
        end


    end

