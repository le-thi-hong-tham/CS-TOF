function noise = createNoise(signal, snrDb)
            snr = db2pow(snrDb);
            [m, n] = size(signal);
            noise = randn(m, n);
            % we treat each column as a separate data vector
            for i=1:n
                signalNorm = norm(signal(:,i));
                noiseNorm = norm(noise(:,i));
                actualNormRatio = signalNorm / noiseNorm;
                requiredNormRatio = sqrt(snr);
                factor = actualNormRatio  / requiredNormRatio;
                noise(:,i) = factor .* noise(:,i);
            end