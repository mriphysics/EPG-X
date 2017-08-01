function Xi = Xi_diff_BM(Xi_T,Xi_L,diff,kmax,N)
% Function to build Xi operator with diffusion effects for BM (6 states per
% k-value). 
% This function efficiently uses spdiags to build the large band diagonal
% matrix that results from applying Xi_L and Xi_T to each order of states.
% Diffusion changes the matrices applied to each order, which this function
% takes care of.
% Shaihan Malik July 2017
[bDL, bDT] = EPG_diffusion_weights(diff.G,diff.tau,diff.D,0:kmax);

% Xi has 3 non-zero diagonals, number -3,0,3. can construct it separately based on these
diags = zeros(N,3);

%%% Diagonal 0
XiT = diag(Xi_T,0)*bDT';
XiL = diag(Xi_L,0)*bDL';
% interleave
XiA = cat(1,XiT(1:2,:),XiL(1,:),XiT(3:4,:),XiL(2,:));
diags(:,1) = XiA(:);

%%% Diagonal -3
XiT = diag(Xi_T,-2)*bDT'; %<-- it's diagonal -2 of XiT
XiL = diag(Xi_L,-1)*bDL'; %<-- and diagonal -1 of XiL
XiA = zeros(6,N/6);
XiA(1:2,:) = XiT;
XiA(3,:) = XiL;
% rest are zero
diags(:,2) = XiA(:);

%%% Diagonal +3
XiT = diag(Xi_T,2)*bDT'; %<-- it's diagonal +2 of XiT
XiL = diag(Xi_L,1)*bDL'; %<-- and diagonal +1 of XiL
XiA = zeros(6,N/6);
XiA(4:5,:) = XiT;% need to put leading zeros in here
XiA(6,:) = XiL;
% rest are zero
diags(:,3) = XiA(:);

%%% Now use sparse diagonal function to define overall Xi
Xi = spdiags(diags,[0 -3 3],N,N);

end