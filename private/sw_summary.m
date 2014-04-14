function output = sw_summary(cfg)

nsubj = 14;

n_f2b = zeros(nsubj, 1);
n_b2f = zeros(nsubj, 1);

v_f2b = cell(nsubj, 1);
v_b2f = cell(nsubj, 1);

for subj = 1:nsubj

  f2b = [];
  b2f = [];
  
  trig_dir = sprintf('%s%04d/spm/triggers/',  cfg.data, subj);
  trig_file = [trig_dir cfg.trigB];
  load(trig_file)
  
  eeg_dir = sprintf('%s%04d/eeg/', cfg.data, subj);
  chk = dir([eeg_dir 'chk*.mat']);
  
  for i = 1:numel(chk)
    chk_file = [eeg_dir chk(i).name];
    warning off
    load(chk_file, 'D')
    warning on
    
    all_sw = [D.other.CRC.SW.SW.negmax_tp] / D.Fsample;
    
    if numel(all_sw) < cfg.minW
      continue
    end
    
    [~, i_f2b] = intersect(all_sw, bSW_onset{i});
    [~, i_b2f] = intersect(all_sw, sSW_onset{i});
    
    if numel(i_f2b) ~= numel(bSW_onset{i}) ||  numel(i_b2f) ~= numel(sSW_onset{i})
      fprintf('check subj %04d, sess % 2d\n', subj, i)
    end

    f2b = [f2b D.other.CRC.SW.SW(i_f2b)];
    b2f = [b2f D.other.CRC.SW.SW(i_b2f)];
    
  end

  n_f2b(subj) = numel(f2b);  
  n_b2f(subj) = numel(b2f);
  
  % v_f2b(subj) = mean(arrayfun(@(x)(numel(x.electrodes)), f2b));
  % v_b2f(subj) = mean(arrayfun(@(x)(numel(x.electrodes)), b2f));

  % v_f2b(subj) = mean(arrayfun(@(x)(median(x.delays)), f2b));
  % v_b2f(subj) = mean(arrayfun(@(x)(median(x.delays)), b2f));
  
  v_f2b{subj} = cell2mat(arrayfun(@(x)(x.delays), f2b, 'uni', 0));
  v_b2f{subj} = cell2mat(arrayfun(@(x)(x.delays), b2f, 'uni', 0));
  
  % v_f2b(subj) = mean([f2b.amplitude]);  
  % v_b2f(subj) = mean([b2f.amplitude]);
  
  % v_f2b(subj) = max([f2b.amplitude]);  
  % v_b2f(subj) = max([b2f.amplitude]);

  % v_f2b(subj) = mean([f2b.maxslope]);  
  % v_b2f(subj) = mean([b2f.maxslope]);
  
  % v_f2b(subj) = mean([bSW_param{:}]);
  % v_b2f(subj) = mean([sSW_param{:}]);
  
end

[h, p] = ttest(n_f2b - n_b2f);
output = sprintf('N Slow Waves is %d, with p=% 6.4f, f2b % 6.2f (% 6.2f), b2f % 6.2f (% 6.2f)\n', ...
  h, p, mean(n_f2b), std(n_f2b), mean(n_b2f), std(n_b2f));

avg_f2b = cellfun(@mean, v_f2b);
avg_b2f = cellfun(@mean, v_b2f);

[h, p] = ttest(avg_f2b - avg_b2f);
output = sprintf('%sDelay is %d, with p=% 6.4f, f2b % 6.2f (% 6.2f), b2f % 6.2f (% 6.2f)\n', ...
  output, h, p, mean(avg_f2b), std(avg_f2b), mean(avg_b2f), std(avg_b2f));
   
t = 0:1:200;
k_b2f = zeros(numel(t), nsubj);
k_f2b = zeros(numel(t), nsubj);

for i = 1:nsubj
  k_f2b(:, i) = ksdensity(v_f2b{i}, t);
  k_b2f(:, i) = ksdensity(v_b2f{i}, t);
end

m_f2b = mean(k_f2b, 2);
m_b2f = mean(k_b2f, 2);

se_f2b = std(k_f2b, [], 2) / sqrt(nsubj);
se_b2f = std(k_b2f, [], 2) / sqrt(nsubj);

[h, p] = ttest(k_f2b' - k_b2f');

d_k = mean(k_f2b' - k_b2f');
h_adjusted = double(fdr(p) < 0.05) .* sign(d_k);  % fdr in private

h = figure;
hold on
plot(t(h_adjusted == 1), .04 * ones(1, numel(find(h_adjusted == 1))), '+r') 
plot(t(h_adjusted == -1), .04 * ones(1, numel(find(h_adjusted == -1))), '+b') 

errorfill(t, m_f2b', se_f2b', 'r'); % errorfill in private
errorfill(t, m_b2f', se_b2f', 'b');

saveas(h, [cfg.outp 'slow_wave_delays.pdf'])
close(h)
