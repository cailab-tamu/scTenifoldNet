function [L,Lnorm] = sbe_laplacian_matrix(A)
% Get graph Laplacian matrix
%
%   L = laplacian(g)
%
% graph Laplacian matrix is defined by L = D - A, where D is vertex degree
% diagonal matrix and A is adjacency matrix.
%
% See also: adjacency
% L = diag(sum(A)) - A;

% Systems Biology & Evolution Toolbox
% Author: James Cai
% Email: jcai@tamu.edu
% Website: https://github.com/jamesjcai/SBEToolbox_lite

% https://github.com/dtuia/KEMA/blob/7378c0fce50a818c2fb59f5de9344ea5c1929fa4/general_routine/laplacian.m
% https://github.com/KavehFathian/clear/blob/e72318ad52082485442f55ad5a7ee9b989b85677/Algorithms/Helpers/NormalizeLap.m
% https://github.com/hungrydoggy/Pinocchio/blob/5664503b210005fa4f1fc053e237ba3ecf6a7945/skeletonizer/matlab/toolbox/compute_mesh_laplacian.m
% https://github.com/fljohnston1/otto-group-product/blob/02a6f35f8c144ed52a2097a1965d561172f2e701/SpectralClustering.m

%{
A=[  0     1     0     1     0     0     0     0     0
     1     0    -1     1    -1     0     0     1     0
     0    -1     0     0     1     1     0     0     0
     1     1     0     0    -1     0     1     1     0
     0    -1     1    -1     0     1     0    -1     1
     0     0     1     0     1     0     0    -1     1
     0     0     0     1     0     0     0     1     0
     0     1     0     1    -1    -1     1     0    -1
     0     0     0     0     1     1     0    -1     0];

%}

D=sum(abs(A));    %
L=diag(D)-A;
if nargout>1
    D(D~=0)=sqrt(1./D(D~=0));
    D=diag(D);
    % Lnorm=D*L*D;
    Lnorm=eye(size(A,1))-D*A*D;     % L = I-D^-1/2*W*D^-1/2
end



