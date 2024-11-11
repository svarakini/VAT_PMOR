function [loaddistr_lhs,loaddistr_rhs]=displacementDistribution(DisplacemntBoundaryCase, ...
                                                                    u_limit_lhs,l_limit_lhs,nodenum_lhs,...
                                                                    u_limit_rhs,l_limit_rhs,nodenum_rhs)
% svara: 02/12/2021                            
                                                                
switch DisplacemntBoundaryCase
    
    case 'uniform'                                                                
        loaddistr_lhs=ones(1,nodenum_lhs);
        loaddistr_rhs=ones(1,nodenum_rhs);
        
    case 'linear-Symm'                                                                
        loaddistr_lhs=linspace(u_limit_lhs,l_limit_lhs,nodenum_lhs);
        loaddistr_rhs=linspace(u_limit_rhs,l_limit_rhs,nodenum_rhs);
        
    case 'linear-Asymm'                                                                
        loaddistr_lhs=linspace(u_limit_lhs,l_limit_lhs,nodenum_lhs);
        loaddistr_rhs=linspace(l_limit_rhs,u_limit_rhs,nodenum_rhs);
        
    case 'linear-Symm-MidPoint-Low'                                                                
        loaddistr_lhs1=linspace(u_limit_lhs,l_limit_lhs,round(nodenum_lhs/2));
        loaddistr_lhs2=flip(loaddistr_lhs1(1,1:end-1));
        loaddistr_lhs=[loaddistr_lhs1 loaddistr_lhs2];
        
        loaddistr_rhs1=linspace(u_limit_rhs,l_limit_rhs,round(nodenum_rhs/2));
        loaddistr_rhs2=flip(loaddistr_rhs1(1,1:end-1));
        loaddistr_rhs=[loaddistr_rhs1 loaddistr_rhs2];
        
    case 'linear-Symm-MidPoint-High'                                                                
        loaddistr_lhs1=linspace(l_limit_lhs,u_limit_lhs,round(nodenum_lhs/2));
        loaddistr_lhs2=flip(loaddistr_lhs1(1,1:end-1));
        loaddistr_lhs=[loaddistr_lhs1 loaddistr_lhs2];
        
        loaddistr_rhs1=linspace(l_limit_rhs,u_limit_rhs,round(nodenum_rhs/2));
        loaddistr_rhs2=flip(loaddistr_rhs1(1,1:end-1));
        loaddistr_rhs=[loaddistr_rhs1 loaddistr_rhs2];
    
    case 'linear-Asymm-MidPoint-LHS:High-RHS:Low'                                                                
        loaddistr_lhs1=linspace(l_limit_lhs,u_limit_lhs,round(nodenum_lhs/2));
        loaddistr_lhs2=flip(loaddistr_lhs1(1,1:end-1));
        loaddistr_lhs=[loaddistr_lhs1 loaddistr_lhs2];
        
        loaddistr_rhs1=linspace(l_limit_rhs,u_limit_rhs,round(nodenum_rhs/2));
        loaddistr_rhs2=flip(loaddistr_rhs1(1,1:end-1));
        loaddistr_rhs=[loaddistr_rhs1 loaddistr_rhs2];
        
    case 'random'
        loaddistr_lhs=rand(1,nodenum_lhs);
        loaddistr_rhs=rand(1,nodenum_rhs);
    
end