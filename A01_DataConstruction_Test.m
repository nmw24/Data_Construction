function [] = DataComp(a,b)

clear all
clc
addpath('misc');
addpath('data');
addpath('outputSummary')
paths = getPaths;

%% Real VIX Data
filename = 'vixarchive.xlsx';
sheet = 1;
yRange1 = 'E6052:E6072';
yRange2 = 'E6092:E6112';
yRange3 = 'E6155:E6175';
yRange4 = 'E6219:E6239';
yRange5 = 'E6282:E6302';

VIX1 = xlsread(filename, sheet, yRange1); VIX1 = VIX1./100;
VIX2 = xlsread(filename, sheet, yRange2); VIX2 = VIX2./100;
VIX3 = xlsread(filename, sheet, yRange3); VIX3 = VIX3./100;
VIX4 = xlsread(filename, sheet, yRange4); VIX4 = VIX4./100;
VIX5 = xlsread(filename, sheet, yRange5); VIX5= VIX5./100;
% size(VIX1)
% size(VIX2)
% size(VIX3)
% size(VIX4)
% size(VIX5)

Recorded_VIX = [VIX1, VIX2, VIX3, VIX4, VIX5];
%% Data Construction
% size(busdays('1/1/14', '1/31/14'))
% size(busdays('3/1/14', '3/31/14'))
% size(busdays('6/1/14', '6/30/14'))
% size(busdays('9/1/14', '9/30/14'))
% size(busdays('12/1/14', '12/30/14'))

if nargin == 0
    %datesUse = busdays('1/1/14', '1/31/14');
    Dates = [busdays('1/1/14', '1/31/14'), busdays('3/1/14', '3/31/14'), busdays('6/1/14', '6/30/14'), busdays('9/1/14', '9/30/14'), busdays('12/1/14', '12/30/14')];
end

% startdate   = a;
% enddate     = b;



frequency = 30;                                                                                    % only works in particular cases
VIX_diff = NaN*ones(size(Dates));

for y = 1:1
    datesUse = Dates(:, y);
    Act_VIX = Recorded_VIX(:, y);
    
    %% Initialise empty variables
    QuoteTime = [];
    ExpTime = [];
    F = [];
    F = [];
    TT = [];
    straddle_VIX = [];
    %% Iterate over each business day
    for k = 1:length(datesUse)                                                                          % This for loop iterates over all the different trading days in our sample
        tic;
        current_date = datestr(datesUse(k),26);
        disp(datestr(datesUse(k), 26))                                                                  % Displays the current trading day in format 'yyyy/mm/dd'
        
        %try
        load([paths.SPXfolder '/outputSummary/SPXVol_' datestr(datesUse(k), 'yyyymmdd') '.mat']);
        
        uniqueTimes = transpose(unique(QuoteDateTime));                                                 % unique(QuoteDateTime) will find all the unique times of quotes given in the current day
        
        
        filtertime = (mod(minute(uniqueTimes), frequency)==0);                                          % minute(uniqueTimes) recieves a serial date and returns the minute of the hour.
        % The mod(minute(uniqueTimes), frequency) will return 0 if frequency divides the minute(uniqueTimes) without a remainder
        uniqueTimes = uniqueTimes(filtertime);                                                          % This establishes the unique times which have been filtered by the frequency
        
        %% Iterate over each frequency in a given day
        for i = 1:length(uniqueTimes)                                                                   % This for loop iterates over the unique times in each trading day
            indx =  find(ismember(QuoteDateTime,uniqueTimes(i)));                                      % ismember() will find a vector the same size as its first input with 1's where the two matricies are equal and 0's elsewhere
            % The find() function which returns a vector containing the index positions of all the nonzero elements of its input, in this case that find the index positions of all the nonzero elements of when quote Time equal unique time
            
            currentTime = uniqueTimes(i);
            
            
            %% initalise variables
            Ftmp = nan(length(indx),1);
            rtmp = nan(length(indx),1);
            Ttmp = nan(length(indx),1);
            Vtmp = nan(length(indx),1);
            
            mat = [];
            calc_VIX = [];
            Datatmp = [];
            for j = 1:length(indx)                                                                      % This for loop iterates over every value in the index, that is each unique time for which there is a quote
                currentExp = ExpiryDateTime{indx(j)};                                                       % In each iteration this sets the value for the current expiry date
                
                r = rDTMTime(indx(j));                                                                  % Defines the current value of the risk-free-rate
                T = (datenum(currentExp) - datenum(currentTime))/365;                                   % Calculates the time to maturity, we divide by 365 and not 252 as time to maturity counts non-trading days
                opdata = OptionMatrix{indx(j)};                                                         % Defines the option data for this iteration
                
                [~, Ftmp_] = GetFutures(opdata, r, T);                                                  % Function GetFutures() calculates the futures price, both that of interpolation for an ATM futures and the average futures price
                K = opdata(:,1);                                                                        % Defines all the strike prices
                
                %% calculate VIX index for the current time to maturity
                clf;
                InTheMoneyCall = opdata(:,1)>Ftmp_;
                InTheMoneyPut = opdata(:,1)<Ftmp_;
                
                %% Calls Monotonicity Filter
                Calls = opdata(~InTheMoneyCall,2);
                Strike_Call = opdata(~InTheMoneyCall,1);
                index = find(~isnan(Calls) == 1);
                Calls = Calls(index);
                Strike_Call = Strike_Call(index);
                
                midway_indx = floor(length(Calls)/2);
                Calls_Lower = Calls(midway_indx+1:end);
                Strike_Lower = Strike_Call(midway_indx+1:end);
                Calls_Upper = Calls(1:midway_indx);
                Strike_Upper = Strike_Call(1:midway_indx);
                
                counter = 0;
                Calls_indx = [];
                for y = 1:(length(Calls_Lower) - 2)
                    %if (Calls_Upper(y) > Calls_Upper(y+1)) == 1
                    if and(Calls_Lower(y) > Calls_Lower(y+1), Calls_Lower(y) > Calls_Lower(y+2)) == 1
                        Calls_indx = [Calls_indx; 1];
                    else
                        Calls_indx = [Calls_indx; 0];
                        counter = counter + 1;
                    end
                    
                end
                if counter ~= 0
                    disp(['Total of ' num2str(counter) ' violations of monotonicity, calls lower'])
                    index = find(Calls_indx == 1);
                    Calls_Lower = Calls_Lower(index); Strike_Lower = Strike_Lower(index);
                else
                end   
                
                counter = 0;
                Calls_indx = [];
                Calls_Upper = flipud(Calls_Upper); Strike_Upper = flipud(Strike_Upper);
                for Y = 1:(length(Calls_Upper) - 2)
                    %if (Calls_Upper(y) < Calls_Upper(y+1)) == 1
                    if and(Calls_Upper(Y) < Calls_Upper(Y+1), Calls_Upper(Y) < Calls_Upper(Y+2)) == 1
                        Calls_indx = [Calls_indx; 1];
                    else
                        Calls_indx = [Calls_indx; 0];
                        counter = counter + 1;
                    end
                    
                end
                if counter ~= 0
                    disp(['Total of ' num2str(counter) ' violations of monotonicity, calls upper'])
                    index = find(Calls_indx == 1);
                    Calls_Upper = flipud(Calls_Upper(index)); Strike_Upper = flipud(Strike_Upper(index));
                else
                    Calls_Upper = flipud(Calls_Upper); Strike_Upper = flipud(Strike_Upper);
                end

               
                Calls = [Calls_Upper; Calls_Lower];
                Strike_Call = [Strike_Upper; Strike_Lower];
                
                if length(Calls) > 1
                    if Calls(1) < Calls(2)
                        Calls(1) = []; Strike_Call(1) = [];
                    else
                    end
                    
                    if Calls(length(Calls)) > Calls(length(Calls)-1)
                        Calls(length(Calls)) = []; Strike_Call(length(Strike_Call)) = [];
                    else
                    end
                else
                end
                

                
                %% Puts Monotonicity Filter
                Puts = opdata(~InTheMoneyPut,3);
                Strike_Put = opdata(~InTheMoneyPut,1);
                index = find(~isnan(Puts) == 1);
                Puts = Puts(index);
                Strike_Put = Strike_Put(index);

                midway_indx = floor(length(Puts)/2);
                Puts_Lower = Puts(midway_indx+1:end);
                Strike_Lower = Strike_Put(midway_indx+1:end);
                Puts_Upper = Puts(1:midway_indx);
                Strike_Upper = Strike_Put(1:midway_indx);
                
                counter = 0;
                Puts_indx = [];
                for g = 1:(length(Puts_Lower) - 2)
                    %if (Puts_Lower(y) < Puts_Lower(y+1)) == 1
                    if and(Puts_Lower(g) < Puts_Lower(g+1), Puts_Lower(g) < Puts_Lower(g+2)) == 1
                        Puts_indx = [Puts_indx; 1];
                    else
                        Puts_indx = [Puts_indx; 0];
                        counter = counter + 1;
                    end
                    
                end
                if counter ~= 0
                    disp(['Total of ' num2str(counter) ' violations of monotonicity, puts lower'])
                    index = find(Puts_indx == 1);
                    Puts_Lower = Puts_Lower(index); Strike_Lower = Strike_Lower(index);
                else
                end

                
                counter = 0;
                Puts_Upper = flipud(Puts_Upper); Strike_Upper = flipud(Strike_Upper);
                Puts_indx = [];
                for J = 1:(length(Puts_Upper) - 2)
                    %if (Puts_Upper(y) > Puts_Upper(y+1)) == 1
                    if and(Puts_Upper(J) > Puts_Upper(J+1), Puts_Upper(J) > Puts_Upper(J+2)) == 1
                        Puts_indx = [Puts_indx; 1];
                    else
                        Puts_indx = [Puts_indx; 0];
                        counter = counter + 1;
                    end
                    
                end
                if counter ~= 0
                    disp(['Total of ' num2str(counter) ' violations of monotonicity, puts upper'])
                    index = find(Puts_indx == 1);
                    Puts_Upper = flipud(Puts_Upper(index)); Strike_Upper = flipud(Strike_Upper(index));
                else
                    Puts_Upper = flipud(Puts_Upper); Strike_Upper = flipud(Strike_Upper);
                end
                Puts = [Puts_Upper; Puts_Lower];
                Strike_Put = [Strike_Upper; Strike_Lower];
                
                if length(Puts) > 1
                    if Puts(1) > Puts(2)
                        Puts(1) = []; Strike_Put(1) = [];
                    else
                    end
                    
                    if Puts(length(Puts)) < Puts(length(Puts)-1)
                        Puts(length(Puts)) = []; Strike_Put(length(Strike_Put)) = [];
                    else
                    end
                else
                end
                

                
                %% Plots
                
%                 plot(Strike_Call./Ftmp_, Calls);
%                 hold on;
%                 plot(Strike_Put./Ftmp_, Puts);
%                 xlim([0 2.0])
%                 title([datestr(currentTime,30) ' DTM = ' num2str(T*365)]);
%                 shg;
%                 pause(0.1)
%                 continue;
                
                % determine OTM options
                % OTM = (K<=FATMtmp_)*2 + (K>FATMtmp_)*1;
                
                %[opdata(~InTheMoneyPut,1) opdata(~InTheMoneyPut,3)]
                
                opdatatmp = [[opdata(~InTheMoneyCall,1);opdata(~InTheMoneyPut,1)] [0*opdata(~InTheMoneyCall,1)+1;0*opdata(~InTheMoneyPut,1)+2] [opdata(~InTheMoneyCall,2); opdata(~InTheMoneyPut,3)]];                           % Defines the new option data in a matrix where the first column is the strike prices, second is calls, third is put (second and third columns will have a 1 or 0 for true or false), then column four gives the price
                opdatatmp2 = [[Strike_Call; Strike_Put] [0*Strike_Call+1;0*Strike_Put+2] [Calls; Puts]];     
                opdatatmp = opdatatmp(sum(isnan(opdatatmp),2)==0,:);                                    % remove NaN
                isequal(opdatatmp, opdatatmp2)
                %[Strike_Put Puts]
                
                TT = [TT; T*365];
                %calculateIntegral( func, K, PC, Price, r, T, F, Klevels )
                try
                    VIX = sqrt(2/T*exp(r*T)*calculateIntegral_test(1, @(K_,F_) 1./K_.^2, opdatatmp(:,1), opdatatmp(:,2), opdatatmp(:,3), r, T, Ftmp_ ));       % Calculates VIX values
                catch
                    disp('went into NaN in construction')
                    VIX = nan;
                end
                calc_VIX = [calc_VIX; VIX];
                diff = Act_VIX(j) - VIX;
                diff = abs(diff);
                VIX_diff(k, y) = diff;%[VIX_diff; diff];
                [TT(j) VIX diff]
                mat = [mat; TT(j)];
                %% data to save
                Ftmp(j) = Ftmp_;
                rtmp(j) = r;
                Ttmp(j) = T;
                Vtmp(j) = VIX;
                Datatmp{j,1} = opdata;
            end
            weight1 = (mat(2)-30)/(mat(2)-mat(1));
            weight2 = 1-weight1;
            straddle_VIX = [straddle_VIX;weight1*calc_VIX(1)+weight2*calc_VIX(2)];
            
            %
            %
            %         Data(counter).QuoteTime = currentTime;
            %         Data(counter).F = Ftmp;
            %         Data(counter).r = rtmp;
            %         Data(counter).T = Ttmp;
            %         Data(counter).VIX = Vtmp;
            %         Data(counter).OpData = Datatmp;
            %
            %         counter = counter + 1;                                                                      % Adds one to counter so each day doesn't overwrite the previous
            %         save(['data/test_data_4.mat'], 'Data', 'current_date');
        end
        
        %     catch
        %        disp(['Date: ' datestr(datesUse(k), 'yyyymmdd') ' Is not recorded.']);
        %     end
        
        
        toc
    end
    
    
end

Act_VIX
straddle_VIX

size(Act_VIX)
size(straddle_VIX)

%plot(datesUse, Act_VIX,'b--o', datesUse, straddle_VIX, 'g--*')


 
% weight_short = NaN*ones(length(mat)-1,1); weight_long = weight_short; 
% VIX =  NaN*ones(length(mat)-2,1);
% 
% t = mat(1);
% for i=2:length(mat)
%     T = mat(i);
%     weight_short(i-1) = (T-30)/(T-t);
%     weight_long(i-1) = (30-t)/(T-t);
% end
% [weight_short weight_long]
% for i=1:length(weight_short)-1
%     VIX(i) = calc_VIX(1).*weight_short(i) + calc_VIX(i+1).*weight_long(i);
% end
% VIX

% mat = unique(mat)
% %% Plots
% subplot(3,2,1)
% plot(Dates(:, 1), VIX_diff(:, 1), 'b--o')
% datetick('x',19)
% ylabel('Difference/100');
% title(['VIX difference 2007, maturity = ' num2str(mat(1))])
%
%
% subplot(3,2,2)
% plot(Dates(:, 2), VIX_diff(:, 2), 'b--o')
% datetick('x', 19)
% ylabel('Difference/100');
% title(['VIX difference 2007, maturity = ' num2str(mat(2))])
%
%
% subplot(3,2,3)
% plot(Dates(:, 3), VIX_diff(:, 3), 'b--o')
% datetick('x', 19)
% ylabel('Difference/100');
% title(['VIX difference 2007, maturity = ' num2str(mat(3))])
%
%
% subplot(3,2,4)
% plot(Dates(:, 4), VIX_diff(:, 4), 'b--o')
% datetick('x', 19)
% ylabel('Difference/100');
% title(['VIX difference 2007, maturity = ' num2str(mat(4))])
%
%
%
% subplot(3,2,5)
% plot(Dates(:, 5), VIX_diff(:, 5), 'b--o')
% datetick('x', 19)
% ylabel('Difference/100');
% title(['VIX difference 2007, maturity = ' num2str(mat(5))])





















