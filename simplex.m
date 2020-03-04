


function [stat, sol] = simplex(c,A,b,BAS)

    [m n] = size(A);

    % compute N
    N = [];
    in = 0;
    for i = 1:n
        for j = 1:m
            if BAS(j,1) == i;
                in = 1;
            end
        end
        if in == 0c
            N = [N; i];
        else
            in = 0;
        end
    end


    % c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cB = [];
    for i = 1:m
        cB(i,1) = c(BAS(i,1),1);
    end
    cN = [];
    for i = 1:n-m
        cN(i,1) = c(N(i,1),1);
    end
    % A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AB = [];
    for i = 1:m
        AB(:,i) = A(:,BAS(i,1));
    end
    AN = [];
    for i = 1:n-m
        AN(:,i) = A(:,N(i,1));
    end
    ABinv = inv(AB);
    % x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xB = mtimes(ABinv,b);
    xN = zeros(n-m,1);



    % begin loop

    % compute the reduced costs
    cBarT = cN.' - mtimes(mtimes(cB.',ABinv),AN);


    while any(cBarT > 0)
        % if reduced costs are strictly less than 0 then STOP,OPTIMAL
        testPos = cBarT > 0;
        for i=1:n-m
            if testPos(1,i) == 1 % else find positive cj
                cj = cBarT(i,1);
                dN = zeros(n-m,1); % determine dN
                dN(i,1) = 1;
                j = N(i,1); % j is next non-basic variable
                i = n-m; % BREAK LOOP
                break;
            end
        end




        % determine dB

        dB = mtimes(-ABinv,A(:,j));
        % construct d
        d = zeros(n,1);
        for i=1:n
            for k=1:m
                if BAS(k,1) == i
                    d(i,1) = dB(k,1);
                end
            end
            for k=1:n-m
                if N(k,1) == i
                    if i == j
                        d(i,1) = 1;
                    else
                        d(i,1) = 0;
                    end
                end
            end
        end


        % if d is non negative the STOP, UNBOUNDED
        if all(dB > 0)
            stat = 1;
            sol = d;
            return;
        end


        % determine Theta and new pivot variable
        ThetaVec = -xB./dB;
        % find minimium positive value
        min = Inf;
        for i=1:m
            if ThetaVec(i,1) > 0 & ThetaVec(i,1) < min
                min = ThetaVec(i,1);
                l = BAS(i,1);
            end
        end
        theta = min;


        % PIVOT to new basis
        for i=1:m
            if BAS(i,1) == l
                BAS(i,1) = j;
                i = m;
            end
        end
        for i=1:n-m
            if N(i,1) == j
                N(i,1) = l;
                i = n-m;
                break;
            end
        end
        BAS = sort(BAS);
        N = sort(N);

        % re create vectors and matrices

        % c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cB = [];
        for i = 1:m
            cB(i,1) = c(BAS(i,1),1);
        end
        cN = [];
        for i = 1:n-m
            cN(i,1) = c(N(i,1),1);
        end
        % A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AB = [];
        for i = 1:m
            AB(:,i) = A(:,BAS(i,1));
        end
        AN = [];
        for i = 1:n-m
            AN(:,i) = A(:,N(i,1));
        end
        ABinv = inv(AB);
        % x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xB = mtimes(ABinv,b);
        xN = zeros(n-m,1);

        % compute the reduced costs
        cBarT = cN.' - mtimes(mtimes(cB.',ABinv),AN);

    end

    stat = 0;

    % construct x
    x = zeros(n,1);
    for i=1:n
        for j=1:m
            if BAS(j,1) == i
                x(i,1) = xB(j,1);
            end
        end
        for j=1:n-m
            if N(j,1) == i
                x(i,1) = 0;
            end
        end
    end

    sol = x;


end
