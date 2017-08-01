function Xi = Xi_diff_MT(Xi_T,Xi_L,diff,kmax,N)
% Function to build Xi operator with diffusion effects for MT (4 states per
% k-value). 
% This function efficiently uses spdiags to build the large band diagonal
% matrix that results from applying Xi_L and Xi_T to each order of states.
% Diffusion changes the matrices applied to each order, which this function
% takes care of.
% Shaihan Malik July 2017
[bDL, bDT] = EPG_diffusion_weights(diff.G,diff.tau,diff.D,0:kmax);

% Xi has 3 non-zero diagonals, number -1,0,1. can construct it separately based on these
diags = zeros(N,3);

%%% Diagonal 0
XiT = diag(Xi_T,0)*bDT';
XiL = diag(Xi_L,0)*bDL';
% Combine them
XiA = cat(1,XiT,XiL);
diags(:,1) = XiA(:);

%%% Diagonal -1
XiL = Xi_L(2,1)*bDL;
diags(3:4:end,2) = XiL(:);

%%% Diagonal +1
XiL = Xi_L(1,2)*bDL;
diags(4:4:end,3) = XiL(:);

%%% Now use sparse diagonal function to define overall Xi
Xi = spdiags(diags,[0 -1 1],N,N);

end