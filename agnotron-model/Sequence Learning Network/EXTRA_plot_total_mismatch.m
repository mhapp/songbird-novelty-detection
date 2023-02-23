bout_length = winSize/2; % we just happen to know this
num_bouts = length(total_mm) / bout_length;


boutwise_mm = zeros(num_bouts,1);

for bout_idx = 1:num_bouts
    bout_start = ((bout_idx-1)*bout_length) + 1;
    bout_end = bout_start + 19;
    bout_mm = sum(total_mm(bout_start:bout_end));
    boutwise_mm(bout_idx) = bout_mm;
end

c1_boutwise_mm = boutwise_mm(1:2:end);
c2_boutwise_mm = boutwise_mm(2:2:end);

figure(808);
subplot(211)
plot(c1_boutwise_mm)
hold on;
plot(196, c1_boutwise_mm(196), 'r*', 'MarkerSize', 12)
hold off;
title('Mismatch Response - Song 1')
ylabel('Response Strength')
makepretty;

subplot(212)
plot(c2_boutwise_mm)
hold on;
plot(196, c2_boutwise_mm(196), 'r*', 'MarkerSize', 12)

hold off;
title('Mismatch Response - Song 2')
ylabel('Response Strength')
xlabel('Bout number')
makepretty;