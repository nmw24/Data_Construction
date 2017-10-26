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

Recorded_VIX = [VIX1, VIX2, VIX3, VIX4, VIX5];

if nargin == 0
    Dates = [busdays('1/1/14', '1/31/14'), busdays('3/1/14', '3/31/14'), busdays('6/1/14', '6/30/14'), busdays('9/1/14', '9/30/14'), busdays('12/1/14', '12/30/14')];
end


frequency = 30;                                                                                    % only works in particular cases
VIX_diff = NaN*ones(size(Dates));


%% chooser
month = 1;
datesUse = Dates(:, month);
disp(['Number of days in chosen month: ' num2str(length(datesUse))]);
Act_VIX = Recorded_VIX(:, month);

%% Initialise empty variables
QuoteTime = [];
ExpTime = [];
F = [];
F = [];
TT = [];
quote_times = [];
straddle_VIX = [];
calc_VIX = NaN*ones(14,7);

%%
day_number = 1;
%% Iterate over each business day
for k = day_number:day_number                                                                        % This for loop iterates over all the different trading days in our sample
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
        quote_times = [quote_times; currentTime];
        
        
        %% initalise variables
        Ftmp = nan(length(indx),1);
        rtmp = nan(length(indx),1);
        Ttmp = nan(length(indx),1);
        Vtmp = nan(length(indx),1);
        
        mat = [];
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
            
            
            [Calls, Strike_Call, Puts, Strike_Put] = Monotonicity_Filter(opdata(~InTheMoneyCall,2), opdata(~InTheMoneyCall,1), opdata(~InTheMoneyPut,3),opdata(~InTheMoneyPut,1));
            
            
            %% Plots
            
            %                 plot(Strike_Call./Ftmp_, Calls);
            %                 hold on;
            %                 plot(Strike_Put./Ftmp_, Puts);
            %                 xlim([0 2.0])
            %                 title([datestr(currentTime,30) ' DTM = ' num2str(T*365)]);
            %                 shg;
            %                 pause(0.1)
            %                 continue;
            
           %% VIX Calculations 
            opdatatmp = [[Strike_Call; Strike_Put] [0*Strike_Call+1;0*Strike_Put+2] [Calls; Puts]];
            
            TT = [TT; T*365];
            try
                VIX = sqrt(2/T*exp(r*T)*calculateIntegral_test(1, @(K_,F_) 1./K_.^2, opdatatmp(:,1), opdatatmp(:,2), opdatatmp(:,3), r, T, Ftmp_ ));       % Calculates VIX values
            catch
                disp('went into NaN in construction')
                VIX = nan;
            end
            calc_VIX(i,j) = VIX;
            diff = Act_VIX(day_number) - VIX;
            diff = abs(diff);
            VIX_diff(k, month) = diff;
            if diff > 0.03
                [TT(j) VIX diff]
            else
            end
            
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
    end
end


Act_VIX = Act_VIX(1)*ones(length(straddle_VIX), 1);
quote_times = [8.30;9;9.30;10;10.30;11;11.30;12;12.30;13;13.30;14;14.30;15];

clf
plot(quote_times, Act_VIX, 'b--o');
hold on;
plot(quote_times, calc_VIX(:,1), 'g--*');
hold on;
plot(quote_times, calc_VIX(:,2), 'b--*');
hold on;
plot(quote_times, calc_VIX(:,3), 'r--*');
hold on;
plot(quote_times, calc_VIX(:,4), 'c--*');
hold on;
plot(quote_times, calc_VIX(:,5), 'm--*');
hold on;
plot(quote_times, calc_VIX(:,6), 'y--*');
hold on;
plot(quote_times, calc_VIX(:,7), 'k--*');
legend('VIX', 'DTM = 15', 'DTM = 50', 'DTM = 78', 'DTM = 106', 'DTM = 169', 'DTM = 260', 'DTM = 351');




















