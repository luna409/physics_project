function output = circleFit(X, method)
Xtranspose=X';
[~,datanum]=size(X);
if method==1
    A=[2*Xtranspose(1:datanum,1) 2*Xtranspose(1:datanum,2) ones(datanum,1)];
    b=[Xtranspose(1:datanum,1).^2+Xtranspose(1:datanum,2).^2];
    theta=A\b;
    output=[theta(1);theta(2);sqrt(theta(3)+theta(1)^2+theta(2)^2)];
end
function objectivefunction=J(x)
    objectivefunction=0;
    for i=1:datanum
        objectivefunction=objectivefunction+(norm(X(:,i)-[x(1);x(2)])-x(3))^2;
    end
    end
if method==2
    A=[2*Xtranspose(1:datanum,1) 2*Xtranspose(1:datanum,2) ones(datanum,1)];
    b=[Xtranspose(1:datanum,1).^2+Xtranspose(1:datanum,2).^2];
    theta=A\b;
    tempoutput=[theta(1);theta(2);sqrt(theta(3)+theta(1)^2+theta(2)^2)];
    output=fminsearch(@(x)J(x),tempoutput);
end
end

