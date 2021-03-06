
function [Integral] = calculateIntegral(SplineType, func, K, PC, Price, r, T, F, Klevels )

if nargin == 8
    Klevels = [0.1 6*F];                                                    % Defines integration limits
end

%% Calculate implied volatilities
O = ones(size(PC));                                                         
IV = BlackIV(PC,Price,F*O,K,r*O,T*O);                                       % Calculates the implied volatilites based on the BSM model, function in misc folder

    

indxnonNaN = ~isnan(IV);                                                    % remove NaN arguments from VIX Calculations, this defines a vector of size IV with 1 for not NaN and 0 for a Nan argument
K = K(indxnonNaN);                                                          % Filter strikes based on non-NaN IV's
PC = PC(indxnonNaN);                                                        % Filter opions based on non-NaN IV's
IV = IV(indxnonNaN);                                                        % Filter IV's based on non-NaN IV's

OTM = (K >= F).*(PC == 1) + (K < F).*(PC == 2);                             % Defines all out of the money options, as either the first or second term will be zero because PC is either 1 (call) or 2 (put)

Kdiff = sort(unique(K));                                                    % Creates a vecor of all the unique strikes, note Kdiff is not the same size and K!
IVdiff = NaN*Kdiff;                                                         % initialises a vector, the same size as Kdiff, that is the size of all the unique strikes, with value NaN


for i=1:numel(Kdiff)                                                        % This for loop iterates over all the unique strike prices
    indx = find(K==Kdiff(i));                                               % Creates an index of all the positions of each unique strike
    datatmp  =[IV(indx) OTM(indx)];                                         % Creates vector including IV's and OTM options
    datatmp = sortrows(datatmp,-2);                                         % Sorts the data of IV and OTM in terms of descending (largest first) IV's
    IVdiff(i) = datatmp(1,1);                                               % Takes the largest IV and appends it into the vecor for IVdiff
end


 clf
 plot(Kdiff, IVdiff, 'LineStyle', 'no', 'Marker', 'o')
 xlabel('Strike')
 ylabel('Implied Vol')
 title(['DTM = ' num2str(T*365)])
% error
try
    f = @(k) M(SplineType, k, Kdiff, IVdiff, r, T, F).*feval(func,k,F);
    Integral = quadgk(f,Klevels(1), Klevels(2));
catch
    disp('went into NaN in integral calc')
    Integral = NaN;
end




function Price = M(SplineType, k, K, IV, r, T, F)

switch SplineType
    case 1 % Cubic Spline
        if 0
            IVp = interp1(K,IV,k,'spline', NaN);
            indleft = find(k<min(K));
            IVp(indleft) = IV(1);
            indright = find(k>max(K));
            IVp(indright) = IV(end);
        else
            Mon = K./F;                                                    % Calculates moneyness vector
            beta = regress(IV, [ones(size(Mon)) Mon  Mon.^2 Mon.^3]);      % Fits IV with cubic function
            Monk = (k./F)';                                                         
            Datak = [ones(size(Monk)) Monk  Monk.^2 Monk.^3];              % Vector with moneyness, moneyness squared and moneyness cubed
            IVp = (Datak*beta)';                                           % Calculates the IV fit from the regressed beta
            indleft = find(k<min(K));
            minM = min(K)/F;
            IVp(indleft) = [1 minM minM^2 minM^3]*beta;
            indright = find(k>max(K));
            maxM = max(K)/F;
            IVp(indright) = [1 maxM maxM^2 maxM^3]*beta;
        end
        
    
    case 2 % Cubic Hermite polynomial
        if 0
            IVp = interp1(K,IV,k,'spline', NaN);
            indleft = find(k<min(K));
            IVp(indleft) = IV(1);
            indright = find(k>max(K));
            IVp(indright) = IV(end);
        else
            Mon = K./F;                                                    % Calculates moneyness vector
            beta = regress(IV, [-3.*Mon Mon.^3]);                          % Fits IV with Hermite cubic function
            Monk = (k./F)';
            Datak = [-3.*Monk Monk.^3];                                    % Vector with moneyness, moneyness squared and moneyness cubed
            IVp = (Datak*beta)';                                           % Calculates the IV fit from the regressed beta
            indleft = find(k<min(K));
            minM = min(K)/F;
            IVp(indleft) = [-3.*minM minM^3]*beta;
            indright = find(k>max(K));
            maxM = max(K)/F;
            IVp(indright) = [-3.*maxM maxM^3]*beta;
        end
        
        
    case 3 % Cubic Bessel polynomial
        if 0
            IVp = interp1(K,IV,k,'spline', NaN);
            indleft = find(k<min(K));
            IVp(indleft) = IV(1);
            indright = find(k>max(K));
            IVp(indright) = IV(end);
        else
            Mon = K./F;                                                    % Calculates moneyness vector
            beta = regress(IV, [ones(size(Mon)) 6.*Mon  15.*Mon.^2 15.*Mon.^3]);      % Fits IV with cubic function
            Monk = (k./F)';                                                         
            Datak = [ones(size(Monk)) 6.*Monk  15.*Monk.^2 15.*Monk.^3];              % Vector with moneyness, moneyness squared and moneyness cubed
            IVp = (Datak*beta)';                                           % Calculates the IV fit from the regressed beta
            indleft = find(k<min(K));
            minM = min(K)/F;
            IVp(indleft) = [1 6.*minM 15.*minM^2 15.*minM^3]*beta;
            indright = find(k>max(K));
            maxM = max(K)/F;
            IVp(indright) = [1 6.*maxM 15.*maxM^2 15.*maxM^3]*beta;
        end
        
end
        




O = ones(size(k));                                                          % Defines ones vector the size of the number of strikes
Put = Black(2*O,F*O,k,r*O,T*O,IVp);                                         % Calculates option prices based on BSM model
Call = Black(1*O,F*O,k,r*O,T*O,IVp);
Price = (k>=F).*Call + (k<F).*Put;                                          % makes price vector (for puts and calls)
Price(Price<0) = 0;                                                         % Filters out negative option values 





