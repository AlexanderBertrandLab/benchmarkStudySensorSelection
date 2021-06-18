function [J,G] = optFunSTECS(w,R1,R2,lambda,groupSelector,eps)
% OPTFUNSTECS The optimization function of the STECS optimization problem.
%
%   Input parameters:
%       w [DOUBLE]: the currect solution of the optimization problem
%       R1 [DOUBLE]: the target covariance matrix
%       R2 [DOUBLE]: the interference covariance matrix
%       groupSelector [BINARY]: a nbVariables.nbGroups x nbGroups binary
%           matrix, indicating per group (column) which variables of the
%           covariance matrices belong to that group with ones at the
%           corresponding positions.
%       lambda [DOUBLE]: the regularization parameter
%
%   Output parameters:
%       J [DOUBLE]: the evaluated objection function
%       G [DOUBLE]: the evaluated gradient


% REFERENCES TO BE ADDED
% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% evaluation objective function
l12norm = 0;
for gr = 1:size(groupSelector,2)
   l12norm = l12norm + sqrt(norm(w(groupSelector(:,gr)==1),2)^2+eps); 
end
J = w'*R2*w+1/(w'*R1*w)+lambda*l12norm;

%% evaluation gradient
gradReg = zeros(size(w));
for i = 1:length(gradReg)
    gradReg(i) = w(i)/(sqrt(norm(w(groupSelector(:,groupSelector(i,:)==1)==1),2)^2+eps)) ;
end
G = 2*R2*w - (2*R1*w)/(w'*R1*w)^2+lambda*gradReg;
end