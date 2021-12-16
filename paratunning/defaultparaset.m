function [ sprate, lambdaA, lambdaB, deltaIvalue, itemax, Kextend, bcmax ] = defaultparaset()
% assign hyperparameters
%   sprate, lambdaA, lambdaB: 12*1 array, for hyperparameter tunning
%   deltaIvalue, itemax, Kextend, bcmax: fixed values

itemax = 5000;
Kextend = 2;
bcmax = 10;

sprateshort = [1];
lambdaAshort = [5,10,20];
lambdaBshort = [0.1,0.2];
paralist = zeros(length(sprateshort)*length(lambdaAshort)*length(lambdaBshort),3);
t = 1;
for i=1:length(sprateshort)
    for j=1:length(lambdaAshort)
        for k=1:length(lambdaBshort)
            paralist(t,1) = sprateshort(i);
            paralist(t,2) = lambdaAshort(j);
            paralist(t,3) = lambdaBshort(k);
            t = t+1;
        end
    end
end

sprate = paralist(:,1);
lambdaA = paralist(:,2);
lambdaB = paralist(:,3);
deltaIvalue = 10;

end

