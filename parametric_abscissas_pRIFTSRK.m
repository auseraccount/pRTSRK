    %% To demonstrate that the parametric abscissas of pRIFTSRK are non-decreasing.
    xv = (0:0.1:200)';
    phi = zeros(length(xv), 9);
    phi(:, 1) = ones(size(xv));
    for i = 1:8
        phi(:, i+1) = phi(:, i) + xv.^i/factorial(i);
    end
    TSRK_flagv = [223 323 324 424 425 626 827 1128]; 
    for tsrk = 1:length(TSRK_flagv)
        TSRK_flag = TSRK_flagv(tsrk);  
        stage = floor(TSRK_flag/100);
        step  = mod(floor(TSRK_flag/10), 10);
        order = mod(TSRK_flag, 10);
        if ~exist('SSPIF-TSRK-methods-master')
            fprintf('Download TSRK file from https://github.com/SSPmethods/SSPIF-TSRK-methods');
            urlwrite('https://codeload.github.com/SSPmethods/SSPIF-TSRK-methods/zip/refs/heads/master', 'SSPIF-TSRK-methods-master.zip');
            unzip('SSPIF-TSRK-methods-master.zip', '.');
        end
        tsrkfilename = ['./SSPIF-TSRK-methods-master/eSSPTSRKplus methods/' ...
            num2str(stage) 's' num2str(step) 'k' num2str(order) 'pSSPTSRK+.mat'];
        load(tsrkfilename);
        fprintf('TSRK method loaded: step = %d, stage = %d, order = %d\n', step, stage, order);
        % Transform the TSRK coefficients into Butcher tableau
        % Dtheta = [D;theta(2:-1:1)]; ABhat = [Ahat; Bhat]; AB = [A; B];
        tildeD = [1 0; 0 1; D]; tildeA = [0 zeros(1,size(A,2)); Ahat, A]; tildeB = [Bhat B];
        e = ones(size(A,2)+1,1); L = [1 0]';
        tildeD = [tildeD; theta]; tildeAB = [tildeA; tildeB];
        c = tildeAB*e - tildeD*L; 
        syms x;
        hatD(:, 1, 1) = ones(size(xv)); hatD(:, 2, 2) = ones(size(xv)); orderp = order;
        hatc = zeros(length(xv),stage+2); 
        %% Recurrent approximations
        hatc(:, 1) = (-1)*ones(size(xv)); hatc(:, 2) = zeros(size(xv));
        psi(:, 1) = phi(:, 1); psi(:, 2) = phi(:, orderp);
        for i = 3:length(c)
            psi(:, i) =  tildeD(i,1)*psi(:, 1) + tildeD(i,2)*psi(:, 2);
            for j = 1:i-1
                psi(:, i) = psi(:, i) + tildeAB(i,j).* psi(:, j) .* xv;
            end
            for j = 1:2
                hatD(:, i,j) = 1./psi(:,i) .* tildeD(i,j) .* psi(:, j);
                for k = 1:i-1
                    hatD(:, i,j) = hatD(:,i,j) + 1./psi(:,i) .* xv .* tildeAB(i,k) .* psi(:,k) .* hatD(:,k,j);
                end
            end
            hatc(:,i) = -hatD(:,i,1);
            for j = 1:i-1
                hatAB(:,i,j) = 1./psi(:,i) .* tildeAB(i,j) .* psi(:,j);
                for k = 3: i-1
                    hatAB(:,i,j) = hatAB(:,i,j) + 1./psi(:,i) .* xv .* tildeAB(i,k) .* psi(:,k) .* hatAB(:,k,j);
                end
                hatc(:,i) = hatc(:,i) + hatAB(:,i,j);
            end
        end 
        indv = 1:2:length(xv);
        figure;
        for i =  1:size(hatc, 2)
            plot(xv(indv), hatc(indv, i),'markersize', 4); hold on;
        end
        axis([0 200 -0.05 1]);
        title_str = ['TSRK(' num2str(stage) ',' num2str(order) ')'];
        title(title_str);
        xlabel('$\tau \kappa$', 'interpreter', 'latex', 'fontsize', fs-4);
        ylabel('$\hat{c}_i$','interpreter', 'latex', 'fontsize', fs-4);
    end
