function [hpWi] = getHPWForPulse(pulse)
if ~isrow(pulse)
    pulse = pulse';
end
pulse(pulse < max(pulse)./sqrt(2)) = 0;
pulse(pulse > max(pulse)./sqrt(2)) = 1;
hpWi = diff([0 pulse]);
hpWi = diff([find(hpWi == 1,1) find(hpWi == -1,1, 'last')]);

