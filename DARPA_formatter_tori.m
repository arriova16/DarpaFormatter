%%%% New DARPA formatter
% 1. Convert .RPS files to .Mat files with tables
% 2. Combine tables within electrode across session
% 3. Analyze
    % a. include 1.35 threshold analysis
    % b. include ability to work from darpa and double darpa
    % c. include ability to work with whistlepig and pinot

%% Formatter for Darpa data
% going through threshold, hybrid, and block. Format and get ready for
% processing

tld = 'C:\Users\Somlab\Box\BensmaiaLab\ProjectFolders\DARPA\Data\RawData';
outdir = 'B:\ProjectFolders\DARPA\Data\ProcessedData';

folder_path = dir(tld);
monkey_list = folder_path(3:end);

%% loading rsp files
monkey_struct = struct();

for m = 3:size(monkey_list,1)
    



end
    
    
    %% Loading in proccessed data
tld = 'B:\ProjectFolders\DARPA\Data\ProcessedData';
overwrite = false;

monkey_list = dir(tld);

 

 clearvars -except monkey_list tld

 % Make structure containing all data from that monkey/electrode

 MEStruct = struct('Monkey', [], ...
                'Electrode', [], ...
                'BlockTaskTable', []); mei = 1;

for m = 3:size(monkey_list, 1)

      electrode_list = dir(fullfile(tld, monkey_list(m).name, 'Electrode*'));
      
    for e = 1:size(electrode_list,1)
        %need to grab from list, looks like only from one
        electrode_name_layout = split(electrode_list(e).name, "_", 1);
        electrode_number = (electrode_name_layout(2,1));

        subf_block = fullfile(tld, monkey_list(m).name);
        blocktask_file_list = dir(fullfile(subf_block, '*mat'));
     
        block_task_tables = cell(size(blocktask_file_list, 1),1);

            for t = 1:size(blocktask_file_list,1)
                temp = load(fullfile(blocktask_file_list(t).folder, blocktask_file_list(t).name));
                fname_split = strsplit(blocktask_file_list(t).name, '_');
                % date_col = cell2table(repmat(fname_split(2), [size(temp.block_task_tables, 1), 1]), 'VariableNames', {'Date'});
                % block_task_tables{t} = [temp.block_task_tables]; 
            end
    end
%Assign to struct
    
    MEStruct(mei).Monkey = monkey_list(m).name;
    MEStruct(mei).Electrode = electrode_number;
    MEStruct(mei).BlockTaskTable = temp;
    % mei = mei + 1;

 end

%% Analysis

for i = 1:size(MEStruct, 2)
    %Analyze the block task
    if ~isempty(MEStruct(i).BlockTaskTable)
        u_dates = unique(MEStruct(i).BlockTaskTable.Date);
        [date_dt_tables, date_dp_tables, date_coeffs] = deal(cell(size(u_dates)));
        for d = 1:length(u_dates)
            d_idx = strcmp(MEStruct(i).BlockTaskTable.Date, u_dates{d});
            [date_dt_tables{d}, date_dp_tables{d}, date_coeffs{d}] = AnalyzeHybridTable(MEStruct(i).BlockTaskTable(d_idx,:));
        end
        %All days
        [across_day_dt_table, across_day_dp_table, across_day_coeffs] = AnalyzeHybridTable(MEStruct(i).BlockTaskTable);
        %Allocate
        MEStruct(i).BlockRates = across_day_dt_table;
        MEStruct(i).BlockDPrime = across_day_dp_table;
        MEStruct(i).BlockCoeffs = across_day_coeffs;
        MEStruct(i).DailyBlockRates = date_dt_tables;
        MEStruct(i).DailyBlockDPrimes = date_dp_tables;
        MEStruct(i).DailyBlockCoeffs = date_coeffs;
    end

end
%% trying to analyze again

 bigtable = MEStruct(1).BlockTaskTable;

for i = 1:size(bigtable)
     u_icms_big = unique(bigtable.StimAmp);
    [u_test_amps_big, ~, ia] = unique(bigtable.IndentorAmp);
    [pd_strings_big, dp_strings_big] = deal(cell(1, length(u_icms_big)));
     p_detect_big = zeros([length(u_test_amps_big),length(u_icms_big)]);
     dprime_big = NaN([length(u_test_amps_big),length(u_icms_big)]);
     for u = 1:length(u_icms_big)
         % Initalize arrays
         p_detect = ones([length(u_test_amps_big),1]) * 1e-3;
         dprime = NaN([length(u_test_amps_big),1]);
        for j = 1:length(u_test_amps_big)
             trial_idx_big = ia == j & [bigtable.StimAmp] == u_icms_big(u);
             correct_idx_big = strcmp(bigtable.Response(trial_idx_big), 'correct');
             if u_test_amps_big(j) == 0
                 p_detect(j) = 1 - (sum(correct_idx_big) / sum(trial_idx_big));
             else
                 p_detect(j) = sum(correct_idx_big) / sum(trial_idx_big);
             end
        end
         p_detect_big(:,u) = p_detect; 
       % Compute d'
        pmiss_big = max([p_detect(1), 1e-3]);
         for j = 1:length(dprime)-1
            phit_big = p_detect(j+1);
            if phit_big == 1 % Correct for infinite hit rate
                 phit_big = .999;
             elseif phit_big == 0
                phit_big = 1e-3;
             end
             dprime(j+1) = norminv(phit_big) - norminv(pmiss_big);        
         end
         dprime_big(:,u) = dprime;
% 
%         % Make strings
        pd_strings_big{u} = sprintf('pDetect_%d', u_icms_big(u));
         dp_strings_big{u} = sprintf('dPrime_%d', u_icms_big(u));
     end
   bigtable_DetectionRates_1 = array2table([u_test_amps_big, p_detect_big, dprime_big], 'VariableNames', ['TestAmps', pd_strings_big, dp_strings_big]);

end

%% Plotting
% dprime_threshold = 1.35;
% sigfun = @(c,x) (c(3) .* (1./1 + exp(-c(1).*(x -c(2))))) + c(4);
% export_directory = 'C:\Users\arrio\Box\BensmaiaLab\ProjectFolders\DARPA\Practice\ExportedImages';
% for i = 1:size(MEStruct,2)
%     
%   x_data = MEStruct.BlockDPrime(i)
%   y_data = 
%   z_data = 
% 
% 
% 
% end
