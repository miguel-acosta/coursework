% read data
data = csvread('Data_for_HW_1.csv', 1, 0);

% rename states 1, 2, 5, 6 to 1, 2, 3, 4
data(data(:,4)==5,4)=3;
data(data(:,4)==6,4)=4;

% parse data
userID = data(:,1);
questionID = data(:,2);
chosenAction = data(:,3);
state = data(:,4);

% compute state-dependent choice probabilities
P = nan(4, 4);
for dp=8:11
    for s = 1:4
    actionChoices = chosenAction(questionID==dp & state == s, :);
    actionChoices = actionChoices - min(actionChoices);
    P(dp-7, s) = 1 - mean(actionChoices);
    end
end


