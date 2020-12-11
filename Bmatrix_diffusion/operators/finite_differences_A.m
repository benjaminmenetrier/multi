%function A = finite_differences(ell_tild,N,BC_type)
% Constructs the finites differences matrix of the implicit discretization 
% of the heat equation in dimension 1. 
% - BC_type can be equal to: 1 -> Dirichlet 
%                            2 -> Neumann
%                            3 -> Periodic
% - ell_tild is the normalized Matern length scale
% - N is the dimension of the matrix
%  The output is a square sparse matrix of size N*N. 

function A = finite_differences_A(ell_tild,N,BC_type)


%Diagonal elements
line_indices=1:N;
column_indices=1:N;
A = sparse(line_indices,column_indices,(1+2*ell_tild^2)*ones(N,1),N,N);

%Sub-diagonal elements
line_indices=1:N-1;
column_indices=2:N;
A = A + sparse(line_indices,column_indices,-ell_tild^2*ones((N-1),1),N,N);

%Sup-diagonal
line_indices=2:N;
column_indices=1:N-1;
A = A + sparse(line_indices,column_indices,-ell_tild^2*ones((N-1),1),N,N);

%Boundary conditions
switch(BC_type)
    
    case(1)% Dirichlet, i.e. value=0 at boundary
    BC=sparse(N,N); % empty matrix, A is by default in Dirichlet BCs
    
    case(2)% Neumann, i.e. derivative=0 at the boundary
    line_indices=[1,N];
    column_indices=[1,N];
    BC=sparse(line_indices,column_indices,-ell_tild^2*[1,1],N,N);
    
    case(3)% Periodic, i.e. value and derivative resp. equal on each side
    line_indices=[1,N];
    column_indices=[N,1];
    BC=sparse(line_indices,column_indices,-ell_tild^2*[1,1],N,N);
end
%Boundary conditions application
A=A+BC;


end