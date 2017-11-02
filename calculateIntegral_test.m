
function [Integral, Kdiff_original, IVdiff_original, Kdiff, IVdiff, model] = calculateIntegral(SplineType, func, K, PC, Price, r, T, F, Klevels )

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
%% Strike Filtering
for Y = 1:2
    for z = 1:(length(Kdiff)-1)
        if (Kdiff(z+1) - Kdiff(z)) > 100
            Kdiff(z+1) = NaN; IVdiff(z+1) = NaN;
        else
        end
    end
    indx = find(~isnan(Kdiff) == 1);
    Kdiff = Kdiff(indx); IVdiff = IVdiff(indx);
end

%% Moneyness filter, range filter
moneyness = Kdiff./F;
counter = 0;
indx = ones(length(moneyness),1);
for z = 1:length(moneyness)
   if moneyness(z) > 1.1
       counter = counter + 1;
       indx(z) = 0;
   else
   end
end
if counter ~= 0
    index = find(indx == 1);
    Kdiff = Kdiff(index); IVdiff = IVdiff(index); moneyness = Kdiff./F;
else
end
% Step-size filter
counter = 0;
indx = ones((length(moneyness)-1),1);
for z=1:(length(moneyness)-1)
    step_size = sqrt((moneyness(z)-moneyness(z+1)).^2);
    if step_size > 0.02
        indx(z) = 0;
        counter = counter + 1;
    else
    end
end
if counter ~= 0
    index = find(indx == 1);
    Kdiff = Kdiff(index); IVdiff = IVdiff(index); moneyness = Kdiff./F;
else
end
if length(moneyness) > 0
    if ( moneyness(length(moneyness))  - moneyness(length(moneyness) - 1)  ) > 0.01
        Kdiff(length(moneyness)) = []; IVdiff(length(moneyness)) = [];
    else
    end
    Kdiff_original = Kdiff; IVdiff_original = IVdiff;
else
    Kdiff_original = Kdiff; IVdiff_original = IVdiff;
end

%% Outlier Filtering
p = polyfit(Kdiff, IVdiff, 3);
model = p(1).*Kdiff.^3 + p(2).*Kdiff.^2 + p(3).*Kdiff + p(4);


diff = sqrt((model - IVdiff).^2);
for z = 1:length(diff)
   if (100*diff(z)) > 3
       IVdiff(z) = NaN; Kdiff(z) = NaN; 
   else
   end
       
end
index = find(~isnan(IVdiff) == 1);
IVdiff = IVdiff(index); Kdiff = Kdiff(index); 

%% Monotonicity Filtering
min_indx = find(IVdiff == min(IVdiff));
left_IV = flipud(IVdiff(1:min_indx));
right_IV = IVdiff((min_indx + 1):length(IVdiff));
K_left = flipud(Kdiff(1:min_indx));
K_right = Kdiff((min_indx + 1):length(Kdiff));

% Left side
counter = 0;
indx = [];
for z = 1:(length(left_IV)-8)
   if and(and(and(and(left_IV(z) < left_IV(z+1), left_IV(z) < left_IV(z+2)),and(left_IV(z) < left_IV(z+3), left_IV(z) < left_IV(z+4))),and(left_IV(z) < left_IV(z+5), left_IV(z) < left_IV(z+6))), and(left_IV(z) < left_IV(z+7), left_IV(z) < left_IV(z+8))) == 1
      indx = [indx;1]; 
   else
       indx = [indx;0];
       counter = counter+1;
   end
end
if counter ~= 0
    index = find(indx == 1);
    left_IV = flipud(left_IV(index)); K_left = flipud(K_left(index)); 
else
    left_IV = flipud(left_IV); K_left = flipud(K_left);
end

if (left_IV(length(left_IV)) - left_IV(length(left_IV) - 1)) < 0
    left_IV(length(left_IV)) = []; K_left(length(K_left)) = [];
else
end

% Right side
counter = 0;
indx = [];
for z = 1:(length(right_IV)-6)
   if and(and(and(right_IV(z) < right_IV(z+1), right_IV(z) < right_IV(z+2)),and(right_IV(z) < right_IV(z+3), right_IV(z) < right_IV(z+4))),and(right_IV(z) < right_IV(z+5), right_IV(z) < right_IV(z+6))) == 1
      indx = [indx;1]; 
   else
       indx = [indx;0];
       counter = counter+1;
   end
end
if counter ~= 0
    index = find(indx == 1);
    right_IV = right_IV(index); K_right = K_right(index); 
else
    right_IV = right_IV; K_right = K_right;
end

if (right_IV(length(right_IV)) - right_IV(length(right_IV) - 1)) < 0
    right_IV(length(right_IV)) = []; K_right(length(K_right)) = [];
else
end

% Re-package IVs and strikes
IVdiff = [left_IV; right_IV]; Kdiff = [K_left; K_right];


%% Plot implied vol against strike
% clf;
% plot(Kdiff_original./F, model, 'LineStyle', 'no', 'Marker', 'o')
% hold on
% plot(Kdiff./F, IVdiff, 'LineStyle', 'no', 'Marker', '*')
% hold on
% plot(Kdiff_original./F, IVdiff_original, 'LineStyle', 'no', 'Marker', '+')
% xlabel('Moneyness')
% ylabel('Implied Vol')
% title(['DTM = ' num2str(T*365)])
% legend('best fit', 'filtered data', 'raw data')


% Use model implied vol and original K for integral calc. This is a TEST
% ONLY
%Kdiff = Kdiff_original; IVdiff = model;


%% Continue with integral calc
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





