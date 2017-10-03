function V = integrate2dfunctri(func,C,qrule)
if(qrule==1)
    P1 = C{1};  P2 = C{2};  P3 = C{3};
    tricentroid = 1/3*(P1+P2+P3);
    A = [1, P1; 1, P2; 1, P3];
    area = 0.5*abs(det(A));
    fval = func(tricentroid(1), tricentroid(2));
    V = area*fval;
elseif(qrule==6)
    P1 = C{1};   P2 = C{2};   P3 = C{3};    
    area = 0.5*abs(det([1, P1; 1, P2; 1, P3]));    
    Q = quadrature_rule(6,2); % 6th order triangular quadrature;
    qw = Q(:,1);  qx = Q(:,2); qy = Q(:,3); % Obtain Quadrature points and weights    
    V = 0;
    for q = 1:qrule
        xhat = (P2(1)-P1(1))*qx(q) + (P3(1)-P1(1))*qy(q) + P1(1);
        yhat = (P2(2)-P1(2))*qx(q) + (P3(2)-P1(2))*qy(q) + P1(2);
        fval = func(xhat, yhat);
        V = V + qw(q)*2*area*fval;
    end
end
end