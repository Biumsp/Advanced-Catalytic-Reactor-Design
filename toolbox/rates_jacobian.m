% Enrico Bussetti 210307

function [MC, new_RO] = rates_jacobian(RO)
    % Computes the analytical jacobian matrix of an array of reaction rates
    % when the composition dependece of the rates is expressed through 
    % a power-law (at the moment it does not handle -1 as order)
    %
    % Arguments:
    % RO: Reaction-Orders matrix (Species x Reaction)
    % 
    % Outputs:
    % - MC: (Matrix) Multiplicative Coefficients (coming from the differentiation)
    % - new_RO: (Matrix) Reaction Orders after differentiation
    
    MC = RO;
    new_RO = RO - 1;
   
end
    
    