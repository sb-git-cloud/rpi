function [H_rpi, h_rpi] = compute_rpi_git(A,B,K)
%COMPUTE_RPI_GIT computes the robust invariant control admissible set for
%the closed-loop dynamics e(k+1)=(A-B*K)e(k)+w(k), where K is a stabilizing 
% controller for the error dynamics, and with state and input constraints
% 
%   He*e \leq he (state),  and Hv*v \leq hv (input).
%
% The disturbances w are ssumed to be bounded by
%
%   Hw*w \leq hw.
%
% For the linear feedback gain the input constraints can be expressed as 
% additional constraints on the error e:
%
%   Hv*v \leq hv = -Hv*K*e\leq hv.
%
% This allows us to compute a robust forward invariant set for the closed
% loop without inputs. THe computations are based on 
% 
% [1]
% "Predictive Control for Linear and Hybrid Systems" by Francesco Borrelli, 
% Alberto Bemporad, Manfred Morari, Cambridge University Press, 2017, ISBN
% 1107016886, 9781107016880; and
%
% [2]
% Gilbert, E. G., & Tan, K. (1991). 
% Linear systems with state and control constraints: The theory and 
% application of maximal output admissible sets. IEEE Transactions on
% Automatic Control, 36, 1008–1020.

% State constraint He*e\leq he
He = ;
he = ;

% Define constraint set for control v s.t. Hv*v\leq hv
Hv = ;
hv = ;

% Add to state constraints
He = [He;-Hv*K];
he = [he;hv];

% Bounded set for disturbance w: Hw*w\leq b
Hw = ;
hw = ;

% Init
H_rpi_new = He;
h_rpi_new = he;
isconverged = false;
Acl = A-B*K;

options = optimoptions('linprog','Display','off');
while isconverged == false
    % Old values
    H_rpi = H_rpi_new;
    h_rpi = h_rpi_new;
    % New values
    [H_pre, h_pre] = compute_robust_pre(Acl, H_rpi, h_rpi, Hw, hw, options);
    [H_rpi_new, h_rpi_new] = compute_intersection(H_rpi, h_rpi, H_pre, ...
        h_pre, options);
    
    % Check if converged by checking if new set is subset of old set
    isconverged = is_subset(H_rpi_new,h_rpi_new,H_rpi,h_rpi,options);
end

end

%% Functions
function is_subset = is_subset(H_new,h_new,H,h,options)
    % As in Algorithm 3.2 [2]

    tol = 1e-10;
    is_subset = true;
    for i=1:size(H_new,1)
        [~, f] = linprog(-H_new(i,:)', [H; H_new(i,:)], [h; h_new(i)+1], ...
            [],[],[],[],options);
        if isempty(f) || -f > h_new(i)+tol
            is_subset = false;
            break;
        end
    end
end

function [H_new, h_new] = compute_intersection(H1, h1, H2, h2, options)
    % Compute minimal intersection

    % Possibly non-minimal intersection: combine half spaces
    H_new = [H1;H2];
    h_new = [h1;h2];

    % Minimal representation: remove redundant inequalities
    i=1;
    while i <= size(H_new,1) && size(H_new,1)>1
        % Check row i by cutting it
        Hi = H_new(i,:);
        hi = h_new(i);
        H_new = [H_new(1:i-1,:); H_new(i+1:end,:)];
        h_new = [h_new(1:i-1);h_new(i+1:end)];
        
        % maximize Hi*x under residual constraints (half planes)
        % Add Hi*x \leq hi+1 for bounded solution
        A = [H_new; Hi];
        b_loc = [h_new;hi+1];        
        [~,f] = linprog(-Hi',A,b_loc,[],[],[],[],options);
        

        if -f>hi
            % Row i is not redundant
            H_new = [Hi; H_new];
            h_new = [hi; h_new];
            i = i+1;
        end
    end
end

function [H_pre, h_pre] = compute_robust_pre(A,He,he,Hw,hw, options)
    % Worst case scenario of disturbance given state constraints
    
    H_pre = He*A;
    for i=1:size(he,1)
        [~,h_pre(i)] = linprog(He(i,:)',Hw,hw,[],[],[],[],options);
    end
    h_pre = he+h_pre';
end