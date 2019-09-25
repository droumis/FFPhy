function [out] = JY_getrewards(index, excludetime, data)
% gets the number of rewards for each epoch
% gets the rate of rewards in rewards/second


warning('OFF','MATLAB:divideByZero');

% get the number of rewards
rewards=size(find(data{index(1)}{index(2)}.Run(:,2)>0),1);
duration=data{index(1)}{index(2)}.Stats.duration/10000;
out.rewards=rewards;
out.duration=duration;


warning('ON','MATLAB:divideByZero');