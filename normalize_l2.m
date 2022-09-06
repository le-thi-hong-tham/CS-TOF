 function X = normalize_l2(X)
            % Normalizes all points in X by the column-wise l-2 norm
            norms =sqrt(sum(X .* conj(X), 1));
            X = bsxfun(@rdivide, X, norms + eps);
            %X = spx.norm.scale_columns(X, norms);
        end
