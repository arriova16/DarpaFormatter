%% Load all data
% tld = 'B:\ProjectFolders\DARPA\Data\ProcessedData';
% outdir = 'B:\ProjectFolders\DARPA\Data\FinalData';
% figure_path = 'B:\ProjectFolders\DARPA\Figures';
tld = 'B:\ProjectFolders\DARPA\Data\ProcessedData';
outdir = 'C:\Users\arrio\Box\BensmaiaLab\ProjectFolders\DARPA\Figures';
figure_path = 'C:\Users\arrio\Box\BensmaiaLab\ProjectFolders\DARPA\Figures';

task_type = {'Block', 'Hybrid'};
data = struct(); ii = 1;

 monkey_list = dir(tld); monkey_list = monkey_list(3:end);
 for m = 1:length(monkey_list)
     for t = 1:length(task_type)
         dir_str = fullfile(tld, monkey_list(m).name, ['*', task_type{t}, 'Task*combtable.mat']);
         flist = dir(dir_str); % Get list of full tables of task type in folder
         for f = 1:length(flist) % Go through each file and data and parameters
             % Monkey name
             data(ii).Monkey = monkey_list(m).name;
             % Task type
             data(ii).Task = task_type{t};
             % Raw data
             temp = load(fullfile(flist(f).folder, flist(f).name));
             data(ii).ResponseTable = temp.bigtable;
             fsplit = strsplit(flist(f).name, '_');
             % Electrode numbers
             if contains(fsplit{2}, 'and') % Pairs
                 and_idx = strfind(fsplit{2}, 'and');
                 ee = [str2double(fsplit{2}(1:and_idx-1)), str2double(fsplit{2}(and_idx+3:end))];
             else % Singles
                 ee = str2double(fsplit{2});
             end
             data(ii).Electrodes = sort(ee); % Force ascending order
             ii = ii + 1;
         end
     end
 end

%% Combine tables of matching electrode
cat_data = struct(); ii = 1;
for m = 1:length(monkey_list)
    for t = 1:length(task_type)
        mt_idx = find(strcmp({data.Monkey}, monkey_list(m).name) & strcmp({data.Task}, task_type{t}));
        elec_list = {data(mt_idx).Electrodes};
        max_elecs = max(cellfun(@length, elec_list));
        elec_table = zeros(length(mt_idx), max_elecs);
        for d = 1:length(mt_idx)
            elec_table(d, 1:length(data(mt_idx(d)).Electrodes)) = data(mt_idx(d)).Electrodes;
        end

        [ue, ~, row_idx] = unique(elec_table, 'rows');
        for i = 1:size(ue,1)
            idx = mt_idx(row_idx == i);
            % Monkey name
            cat_data(ii).Monkey = monkey_list(m).name;
            % Task type
            cat_data(ii).Task = task_type{t};
            cat_data(ii).Electrodes = ue(i,ue(i,:) > 0);
            cat_data(ii).ResponseTable = cat(1, data(idx).ResponseTable);
            % Stim parameters
            cat_data(ii).ICMSAmps = unique(cat_data(ii).ResponseTable.StimAmp);
            cat_data(ii).ICMSFreqs = unique(cat_data(ii).ResponseTable.StimFreq);
            cat_data(ii).MechAmps = unique(cat_data(ii).ResponseTable.IndentorAmp);
            cat_data(ii).MechFreqs = unique(cat_data(ii).ResponseTable.IndentorFreq);
            ii = ii + 1;
        end
    end
end
            

%% Get unique ICMS parameters for each entry
dpt = 1.35;
for d = 1:length(cat_data)
    % All combinations of parameters
    param_combs = combvec(cat_data(d).ICMSAmps',...
                          cat_data(d).ICMSFreqs',...
                          cat_data(d).MechFreqs'); % Need to transpose inputs to row vectors
    summary_struct = struct(); ii = 1;
    for p = 1:size(param_combs,2)
        % Find matching idx
        idx = cat_data(d).ResponseTable.StimAmp == param_combs(1,p) & ...
              cat_data(d).ResponseTable.StimFreq == param_combs(2,p) & ...
              cat_data(d).ResponseTable.IndentorFreq == param_combs(3,p);
        if sum(idx) == 0 % Skip empty combinations

            continue
        end
        % Get unique indentor amplitudes
        uia = unique(cat_data(d).ResponseTable.IndentorAmp(idx));
        pdetected = zeros(size(uia));
        for a = 1:length(uia)
            p1u_idx = idx & cat_data(d).ResponseTable.IndentorAmp == uia(a);
            pdetected(a) = sum(strcmp(cat_data(d).ResponseTable.Response(p1u_idx), 'correct')) / sum(p1u_idx);
            if uia(a) == 0
                pdetected(a) = 1 - pdetected(a);
            end
        end
        % Compute d'
        pdetected_adj = pdetected;
        pdetected_adj(pdetected_adj == 0) = 0.01;
        pdetected_adj(pdetected_adj == 1) = 0.99;
        dprime = norminv(pdetected_adj) - norminv(pdetected_adj(1));
        % Compute threshold
        [SigmoidFun, coeffs, ~, ~, ~, warn] = FitSigmoid(uia, pdetected, 'Constraints', [0, 200; -5, 5]);
        xq = linspace(uia(1), uia(end)*2);
        y_fit = SigmoidFun(coeffs,xq);
        % Convert sigmoid to d'
        dprimeq = norminv(y_fit) - norminv(y_fit(1));
        if warn
            figure; hold on
            plot(uia, pdetected)
            plot(xq, y_fit, 'LineStyle','--')
        end
        yq_idx = find(dprimeq >= dpt, 1, 'first');
        mt = xq(yq_idx);
        if isempty(yq_idx)
            mt = NaN;
        end

        % Add to summary struct
        summary_struct(ii).ICMSAmp = param_combs(1,p);
        summary_struct(ii).ICMSFreq = param_combs(2,p);
        summary_struct(ii).IndentorFreq = param_combs(3,p);
        summary_struct(ii).IndenterAmps = uia;
        summary_struct(ii).pDetected = pdetected;
        summary_struct(ii).dPrime = dprime;
        summary_struct(ii).MechThreshold = mt;
        summary_struct(ii).SigFun = SigmoidFun;
        summary_struct(ii).SigCoeffs = coeffs;
        ii = ii + 1;
    end
    cat_data(d).Summary = summary_struct;
end

%% Bootstrap testing
num_permutations = 1e4;
[wb, pfwb_update]  = ParforWaitbar('Bootstrap testing', length(cat_data));
parfor d = 1:length(cat_data)
    % Find list of no stim conditions
    control_idx = find([cat_data(d).Summary.ICMSAmp] == 0);
    num_controls = length(control_idx);
    for c = 1:num_controls
        % Get list of control condition trials
        control_table_idx = [cat_data(d).ResponseTable.StimAmp] == 0 & ...
                            [cat_data(d).ResponseTable.StimFreq] == cat_data(d).Summary(control_idx(c)).ICMSFreq & ...
                            [cat_data(d).ResponseTable.IndentorFreq] == cat_data(d).Summary(control_idx(c)).IndentorFreq;
        % Find treatment groups
        treatment_idx = find([cat_data(d).Summary.ICMSAmp] ~= 0 & ...
                             [cat_data(d).Summary.ICMSFreq] == cat_data(d).Summary(control_idx(c)).ICMSFreq & ...
                             [cat_data(d).Summary.IndentorFreq] == cat_data(d).Summary(control_idx(c)).IndentorFreq);
        num_treatment_groups = length(treatment_idx);
        for t = 1:num_treatment_groups
            % Get observed effect size
            delta_thresholds = abs(cat_data(d).Summary(control_idx(c)).MechThreshold - ...
                                   cat_data(d).Summary(treatment_idx(t)).MechThreshold);
            % Find treatment trials
            treatment_table_idx = [cat_data(d).ResponseTable.StimAmp] == cat_data(d).Summary(treatment_idx(t)).ICMSAmp & ...
                                  [cat_data(d).ResponseTable.StimFreq] == cat_data(d).Summary(control_idx(c)).ICMSFreq & ...
                                  [cat_data(d).ResponseTable.IndentorFreq] == cat_data(d).Summary(control_idx(c)).IndentorFreq;
            % Permute the indices repeatedly and measure effect size on
            % each permutation
            combined_indices = find(control_table_idx | treatment_table_idx);
            treatment_table = cat_data(d).ResponseTable(combined_indices, ["IndentorAmp", "Response"]);
            num_indices = size(treatment_table,1);
            idx_list = 1:num_indices;
            half_length = floor(num_indices / 2);
            % Get unique amplitude
            uia = unique(cat_data(d).ResponseTable.IndentorAmp(combined_indices));
            xq = linspace(uia(1), uia(end)*2);
            % Prepare null vector
            null_delta_threshold = zeros(num_permutations,1);
            for p = 1:num_permutations
                % Get shuffled indices
                perm_idx = randperm(num_indices);
                p1_idx = idx_list(perm_idx(1:half_length));
                p2_idx = idx_list(perm_idx(half_length+1:half_length*2));

                % Sample detection probability from shuffled indices
                [p1detected, p2detected] = deal(zeros(size(uia)));
                for a = 1:length(uia)
                    p1u_idx = treatment_table.IndentorAmp(p1_idx) == uia(a);
                    p1detected(a) = sum(strcmp(treatment_table.Response(p1_idx(p1u_idx)), 'correct')) / sum(p1u_idx);
                    p2u_idx = treatment_table.IndentorAmp(p2_idx) == uia(a);
                    p2detected(a) = sum(strcmp(treatment_table.Response(p2_idx(p2u_idx)), 'correct')) / sum(p2u_idx);
                    if uia(a) == 0
                        p1detected(a) = 1 - p1detected(a);
                        p2detected(a) = 1 - p2detected(a);
                    end
                end
                
                % Fit sigmoid to each pdetected & get threshold
                [~, coeffs1, ~, ~, ~, ~] = FitSigmoid(uia, p1detected, 'Constraints', [0, 300; -5, 5]);
                mt1 = Sigmoid2MechThreshold(coeffs1, xq, dpt);
                [~, coeffs2, ~, ~, ~, ~] = FitSigmoid(uia, p2detected, 'Constraints', [0, 300; -5, 5]);
                mt2 = Sigmoid2MechThreshold(coeffs2, xq, dpt);
                null_delta_threshold(p) = abs(mt1 - mt2);
            end
            % Compute p
            cat_data(d).Summary(treatment_idx(t)).BootP = sum(null_delta_threshold > delta_thresholds) / num_permutations;
%             if cat_data(d).Summary(treatment_idx(t)).BootP < 0.05
%                 disp('Not significant')
%             end
        end
    end
    send(pfwb_update, 0);
end

%% plot
clf; hold on
group_counter = 0;
significant_counter = 0;
for d = 1:length(cat_data)
    control_idx = find([cat_data(d).Summary.ICMSAmp] == 0);
    num_controls = length(control_idx);
    for c = 1:num_controls
        treatment_idx = find([cat_data(d).Summary.ICMSAmp] ~= 0 & ...
                             [cat_data(d).Summary.ICMSFreq] == cat_data(d).Summary(control_idx(c)).ICMSFreq & ...
                             [cat_data(d).Summary.IndentorFreq] == cat_data(d).Summary(control_idx(c)).IndentorFreq);
        num_treatment_groups = length(treatment_idx);
        for t = 1:num_treatment_groups
            delta_thresholds = abs(cat_data(d).Summary(control_idx(c)).SigCoeffs(2) - ...
                                   cat_data(d).Summary(treatment_idx(t)).SigCoeffs(2));
            p_val = cat_data(d).Summary(treatment_idx(t)).BootP;
            scatter(delta_thresholds, p_val, 100,"filled")
            disp(p_val)

            group_counter = group_counter + 1;
            if p_val < 0.05
                significant_counter = significant_counter + 1;
            end
            ax = gca;
            ax.FontSize = 28;
            xlabel('Delta Threshold',"FontSize",28)
            ylabel('P-Value',"FontSize",28)
        end
    end
end


return

%% Individual electrode plots
mrkr_size = 50;
SetFont('Arial', 12)
export = true;
for d = 1%length(cat_data)
    if length(cat_data(d).Electrodes) > 1
        ee_str = sprintf(['Electrodes %d' repmat(',%d', length(cat_data(d).Electrodes) - 1)], cat_data(d).Electrodes);
    else
        ee_str = sprintf('Electrode %d', cat_data(d).Electrodes);
    end
    title_str = sprintf('%s %s-Task, %s', cat_data(d).Monkey, cat_data(d).Task, ee_str);

    amp_colors = [[.6 .6 .6]; parula(length(cat_data(d).ICMSAmps) - 1)];

    clf;
    % Probability of detection
    subplot(1,5,[1,2]); hold on; title(title_str)
    for p = 1:size(cat_data(d).Summary,2)
        % Color is determined by Amplitude
        c_idx = cat_data(d).Summary(p).ICMSAmp == cat_data(d).ICMSAmps;
        c = amp_colors(c_idx,:);
        % LineStyle is determined by Frequency
        if cat_data(d).Summary(p).ICMSFreq == 100
            mrkr_style = 'o';
        elseif cat_data(d).Summary(p).ICMSFreq == 300
            mrkr_style = '^';
        end
        % Plot
        x = cat_data(d).Summary(p).IndenterAmps;
        y = cat_data(d).Summary(p).pDetected;
        plot(x,y, 'Color', c, 'LineWidth', 1)
        scatter(x,y, mrkr_size, c, "filled", mrkr_style)
    end
    set(gca, 'YLim', [0 1], 'YTick', [0,1], 'XTick', round(linspace(0, max(cat_data(d).MechAmps), 5), 2))
    xlabel('Mech Amplitude (mm)')
    ylabel('p(Detected)')

    % Legends
    ct = ColorText(cat_data(d).ICMSAmps, amp_colors);
%     ct{1} = [ct{1}, ' ', GetUnicodeChar('mu'),'A'];
    ct(end+1)= ColorText('uA', [.6 .6 .6]); 
    text(max(cat_data(d).MechAmps) * 0.95, 0.05, ct, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    shape_text = sprintf('%s 100 Hz\n%s 300 Hz', GetUnicodeChar('EmptyCircle'), GetUnicodeChar('UpTriangle'));
    text(max(cat_data(d).MechAmps) * 0.75, 0.05, shape_text, 'Color', [.6 .6 .6], ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontName', 'JuliaMono')


    % dPrime
    subplot(1,5,[3,4]); hold on
    for p = 1:size(cat_data(d).Summary,2)
        % Color is determined by Amplitude
        c_idx = cat_data(d).Summary(p).ICMSAmp == cat_data(d).ICMSAmps;
        c = amp_colors(c_idx,:);
        % LineStyle is determined by Frequency
        if cat_data(d).Summary(p).ICMSFreq == 100
            mrkr_style = 'o';
        elseif cat_data(d).Summary(p).ICMSFreq == 300
            mrkr_style = '^';
        end
        % Plot
        x = cat_data(d).Summary(p).IndenterAmps;
        y = cat_data(d).Summary(p).dPrime;
        plot(x,y, 'Color', c, 'LineWidth', 1)
        scatter(x,y, mrkr_size, c, "filled", mrkr_style)
    end
    set(gca, 'YLim', [-1 5], 'YTick', [-1:5], 'XTick', round(linspace(0, max(cat_data(d).MechAmps), 5), 2))
    xlabel('Mech Amplitude (mm)')
    ylabel('d')

    % delta Threshold
    subplot(1,5,5); hold on
    uf = unique([cat_data(d).Summary.ICMSFreq]);
    for f = 1:length(uf)
        if uf(f) == 100
            mrkr_style = 'o';
       elseif uf(f) == 300
            mrkr_style = '^';
        end 
        uf_idx = [cat_data(d).Summary.ICMSFreq] == uf(f);
        uft = [[cat_data(d).Summary(uf_idx).ICMSAmp];...
               [cat_data(d).Summary(uf_idx).MechThreshold]];
        dt = uft(2,:) - uft(2,1);
        for i = 2:length(dt)
            amp_idx = uft(1,i) == cat_data(d).ICMSAmps;
            scatter(1, dt(i), mrkr_size, amp_colors(amp_idx,:), 'filled', mrkr_style)
        end
    end
    ymax = max(abs(gca().YLim));
    set(gca, 'XLim', [.5 1.5], 'XTick', [], 'XTickLabel', num2str(uf'), 'YLim', [-ymax, ymax], 'YTick', [-ymax, 0, ymax])
    ylabel(sprintf('%s Threshold (mm)', GetUnicodeChar('Delta')))
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 400])

    if export
        ffname = fullfile(figure_path, title_str);
        print(gcf, ffname, '-dpng', '-r300')
    end
end


%% Delta threshold plot
clf; 
subplot(1,3,1); hold on
    plot([0, 1.4], [0, 1.4], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for d = 1:length(cat_data)
        if isfield(cat_data(d).Summary, 'BootP') == 0
            continue
        end
        uf = unique([cat_data(d).Summary.ICMSFreq]);
        % amp_colors = [[.6 .6 .6]; parula(length(cat_data(d).ICMSAmps) - 1)];
        amp_colors = repmat([.2 .2 .2], [length(cat_data(d).ICMSAmps),1]);
        for f = 1:length(uf)
            if uf(f) == 100
                mrkr_style = 'o';
            elseif uf(f) == 300
                mrkr_style = '^';
            end 
            uf_idx = [cat_data(d).Summary.ICMSFreq] == uf(f) & [cat_data(d).Summary.ICMSAmp] > 0;
            
            uft = [[cat_data(d).Summary(uf_idx).ICMSAmp];...
                   [cat_data(d).Summary(uf_idx).MechThreshold];...
                   [cat_data(d).Summary(uf_idx).BootP]];
            for i = 2:size(uft,2)
                amp_idx = uft(1,i) == cat_data(d).ICMSAmps;
                if uft(3,i) < 0.05
                    fa = 0.8;
                else
                    fa = 0.1;
                end
                
                scatter(uft(2,1), uft(2,i), mrkr_size, mrkr_style,...
                    'MarkerEdgeColor', amp_colors(amp_idx,:),'MarkerFaceColor', amp_colors(amp_idx,:), 'MarkerFaceAlpha', fa)
            end
        end
    end
    xlabel('Mech Threhsold (Control)','FontSize',18)
    ylabel('Mech Threhsold (ICMS)','FontSize',18)

subplot(1,3,2); hold on
    plot([0, length(cat_data)], [0, 0], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for d = 1:length(cat_data)
        if isfield(cat_data(d).Summary, 'BootP') == 0
            continue
        end
        uf = unique([cat_data(d).Summary.ICMSFreq]);
        % amp_colors = [[.6 .6 .6]; parula(length(cat_data(d).ICMSAmps) - 1)];
        amp_colors = repmat([.2 .2 .2], [length(cat_data(d).ICMSAmps),1]);
        for f = 1:length(uf)
            if uf(f) == 100
                mrkr_style = 'o';
            elseif uf(f) == 300
                mrkr_style = '^';
            end 
            uf_idx = [cat_data(d).Summary.ICMSFreq] == uf(f);
            uft = [[cat_data(d).Summary(uf_idx).ICMSAmp];...
                   [cat_data(d).Summary(uf_idx).MechThreshold]];
            dt = uft(2,:) - uft(2,1);
            uf_idxf = find(uf_idx);
            bp = cat_data(d).Summary(uf_idxf(2:end)).BootP;
            for i = 2:length(dt)
                amp_idx = uft(1,i) == cat_data(d).ICMSAmps;
                if bp(i) < 0.05
                    fa = 0.8;
                else
                    fa = 0.1;
                end
                scatter(d, dt(i), mrkr_size, mrkr_size, mrkr_style,...
                    'MarkerEdgeColor', amp_colors(amp_idx,:),'MarkerFaceColor', amp_colors(amp_idx,:), 'MarkerFaceAlpha', fa)
            end
        end
    end
    
    ylabel(sprintf('%s Threshold (mm) (lower is better)', GetUnicodeChar('Delta')),'FontSize',18)
    xlabel('Monkey/Electrode(s)','FontSize',18)
    set(gca,'XTick', [], 'YLim', [-.4 .4])

subplot(1,3,3); hold on
    plot([0, length(cat_data)], [1, 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for d = 1:length(cat_data)
        uf = unique([cat_data(d).Summary.ICMSFreq]);
        % amp_colors = [[.6 .6 .6]; parula(length(cat_data(d).ICMSAmps) - 1)];
        amp_colors = repmat([.2 .2 .2], [length(cat_data(d).ICMSAmps),1]);
        for f = 1:length(uf)
            if uf(f) == 100
                mrkr_style = 'o';
           elseif uf(f) == 300
                mrkr_style = '^';
            end 
            uf_idx = [cat_data(d).Summary.ICMSFreq] == uf(f);
            uft = [[cat_data(d).Summary(uf_idx).ICMSAmp];...
                   [cat_data(d).Summary(uf_idx).MechThreshold]];
            dt = uft(2,:) / uft(2,1);
            for i = 2:length(dt)
                amp_idx = uft(1,i) == cat_data(d).ICMSAmps;
                scatter(d, dt(i), mrkr_size, amp_colors(amp_idx,:), 'filled', mrkr_style)
            end
        end
    end
    
    ylabel('Threshold ICMS / Threshold Control', 'FontSize',18)
    xlabel('Monkey/Electrode(s)','FontSize',18)
    set(gca,'XTick', [], 'YScale', 'log', 'YLim', [.2 5], 'YTick', [.2 5])

set(gcf, 'Position', [2600 100 2000 500])



%% Summary plot w/ significance
SetFont('Arial', 12)
clearvars ax
clf; 
 ax = gca;
      ax.FontSize = 30;
% Create axes and format
ax(1) = subplot(1,3,1); hold on
    plot([0, 1.4], [0, 1.4], 'Color', [.6 .6 .6], 'LineStyle', '--')
    xlabel('Mech Threhsold (Control)',"FontSize",24)
    ylabel('Mech Threhsold (ICMS)',"FontSize",24)
ax(2) = subplot(1,3,2); hold on
    plot([0, length(cat_data)], [0, 0], 'Color', [.6 .6 .6], 'LineStyle', '--')
    ylabel(sprintf('%s Threshold (mm) (lower is better)', GetUnicodeChar('Delta')),'FontSize',24)
    xlabel('Monkey/Electrode(s)',"FontSize",24)
    set(gca,'XTick', [], 'YLim', [-.4 .4])
ax(3) = subplot(1,3,3); hold on
    plot([0, length(cat_data)], [1, 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    ylabel('Threshold ICMS / Threshold Control',"FontSize",24)
    xlabel('Monkey/Electrode(s)',"FontSize",24)
    set(gca,'XTick', [], 'YScale', 'log', 'YLim', [.2 5], 'YTick', [.2 5])

for i = 1:length(cat_data)
    % Iterate through unique icms:indentor frequencies
    [unique_conditions, ~, cond_idx] = unique([[cat_data(i).Summary.ICMSFreq];...
                                               [cat_data(i).Summary.IndentorFreq]]', 'rows');
    for c = 1:size(unique_conditions, 1)
        % Set marker style based on frequency
        if unique_conditions(c,1) == 100
            mrkr_style = 'o';
        elseif unique_conditions(c,1) == 300
            mrkr_style = '^';
        end

        % Get all amplitudes > 0
        c_idx = find(cond_idx == c);
        icms_amps = [cat_data(i).Summary(c_idx).ICMSAmp];
        if sum(icms_amps == 0) ~= 1 || length(c_idx) < 2 % Skip no controls or only one group
            continue
        end
        % Get control threshold
        control_threshold = cat_data(i).Summary(c_idx(icms_amps == 0)).MechThreshold;
        for j = 1:length(c_idx)
            if cat_data(i).Summary(c_idx(j)).ICMSAmp == 0 % Skip control
                continue
            end
            % Use p-value for alpha
            if cat_data(i).Summary(c_idx(j)).BootP < 0.05
                sa = 0.8;
            else
                sa = 0.1;
            end
                
            % Plot individual values
            treatment_threhsold = cat_data(i).Summary(c_idx(j)).MechThreshold;
            % Ax1 - control vs treatment
            scatter(control_threshold, treatment_threhsold, 80, mrkr_style, 'MarkerEdgeColor', [.6 .6 .6],...
                'MarkerFaceColor', [.6 .6 .6], 'MarkerFaceAlpha', sa, 'Parent', ax(1))
            % Ax2 - treatment minus control
            scatter(i, treatment_threhsold - control_threshold, 80, mrkr_style, 'MarkerEdgeColor', [.6 .6 .6],...
                'MarkerFaceColor', [.6 .6 .6], 'MarkerFaceAlpha', sa, 'Parent', ax(2))
            % Ax3 - treatment / control
            scatter(i, treatment_threhsold / control_threshold, 80, mrkr_style, 'MarkerEdgeColor', [.6 .6 .6], ...
                'MarkerFaceColor', [.6 .6 .6], 'MarkerFaceAlpha', sa, 'Parent', ax(3))
        end
        
    end
end
%%

function mt = Sigmoid2MechThreshold(coeffs, xq, threshold)
    SigmoidFun = GetSigmoid(length(coeffs));
    y_fit = SigmoidFun(coeffs, xq);
    % Convert sigmoid to d'
    dprimeq = norminv(y_fit) - norminv(y_fit(1));
    yq_idx = find(dprimeq >= threshold, 1, 'first');
    mt = xq(yq_idx);
    if isempty(yq_idx)
        mt = NaN;
    end
end



