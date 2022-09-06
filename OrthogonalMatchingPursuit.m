function result = OrthogonalMatchingPursuit(Dict, K,y)
            % We assume that all the columns in dictionary are normalized.
             % Initialization
            % Solves approximation problem using OMP
            SaveAllApproximations = false
        % Maximum residual norm
        MaxResNorm = 1e-4
        % Indicates if we should stop on exceeding residual norm
        StopOnResidualNorm = true
        % Indicates if we should stop when residual norm stops improving
        StopOnResNormStable = true
        % Maximum number of iterations for approximation
%         MaxIters
        Verbose = false
        % Minimum Sparsity
        MinK = 4
        % Ignored atom (which won't be considered in identification step)
        IgnoredAtom = -1
            [N, D] = size(Dict);
            if ~exist('K', 'var')
                % No sparsity level has been pre-specified.
                K = -1;
            end
            % Maximum number of iterations
            maxIter = N;
            if K > 0
                % We have to consider pre-specified sparsity level
                maxIter = K;
            end
            MaxIters = maxIter;
            if K > 0 && MinK >= K
                MinK = 0;
            end
        
        
            r=y;
      
      
            
            % Active indices 
            omega = [];
            % Estimate
            z = zeros(D, 1);
            if SaveAllApproximations
                ZZ = zeros(D, N);
            end
            if StopOnResNormStable
                oldResNorm = norm(r);
            end
            maxIter = MaxIters;
            ignored_atom =IgnoredAtom;
            for iter=1:maxIter
                % Compute inner products
                AH=Dict';
                innerProducts =AH * r;
                % Mark the inner products of already selected columns as 0.
                innerProducts(omega) = 0;
                if ignored_atom > 0
                    % forcefully ignore this atom
                    innerProducts(ignored_atom) = 0;
                end 
                innerProducts = abs(innerProducts);
                % Find the highest inner product
                [~, index] = max(innerProducts);
                % Add this index to support
                omega = [omega, index];
                % Solve least squares problem
                subdict = Dict(:, omega);
                tmp = linsolve(subdict, y);
                % Updated solution
                z(omega) = tmp;
                if SaveAllApproximations
                    % We need to keep track of all approximations
                    ZZ(omega, iter) = tmp;
                end
                % Let us update the residual.
                r = y - subdict * tmp;
                if StopOnResidualNorm || StopOnResNormStable
                    resNorm = norm(r);
                    if resNorm < MaxResNorm
                        break;
                    end
                    if StopOnResNormStable
                        change = abs(oldResNorm  - resNorm);
                        if change/oldResNorm < .01
                            % No improvement
                            break;
                        end
                    end
                end
            end
            % Solution vector
            result.z = z;
            % Residual obtained
            result.r = r;
            % Number of iterations
            result.iterations = iter;
            % Solution support
            result.support = omega;
            % Iterative solutions
            if SaveAllApproximations
                % Remove extra 0s
                result.Z = ZZ(:, 1:iter);
            end
            result = result.z;
end
            