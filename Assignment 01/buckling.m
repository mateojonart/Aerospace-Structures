function [sigCritica,Fail] = buckling(data,x,Tn,sig,m,Tm)
   
    Fail = zeros(size(sig,1),1);
    l = zeros(data.nel,1);
    for ii = 1:data.nel
        [xel] = [x(Tn(ii,:),:)];
        l(ii) = (sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2));  
    end

    
    sigCritica = zeros(size(sig,1),1);

    for n = 1:size(m,1)
       for jj = 1:size(Tm)
           if Tm(jj) == n
                sigCritica(jj) = -(pi^2 * m(n,1) * m(n,4)) / (l(jj)^2 * m(n,2)); 
           end
       end
    end

    for ii = 1:size(sig,1)
        if sig(ii) < 0
            if abs(sig(ii)) > abs(sigCritica(ii))
                Fail(ii) = 1;
            else  
                Fail(ii) = 0;
            end
        else 
            Fail(ii) = 0;
        end                                 
    end

    % Identify which element fails by buckling.
    m = 0;
    for ii = 1:size(Fail,1)
        if Fail(ii) == 1
        disp(['The element ', num2str(ii), ' fails by buckling.']);
        else 
            m = m+1;
        end
    end

    if m == data.nel
        disp('None of the elements fail by buckling.');
    end

end



