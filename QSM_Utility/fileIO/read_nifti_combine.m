function read_nifti_combine(nifti_dirs, cleanup, parrecflag, dcm_name)
% function read_nifti_combine(nifti_dirs, cleanup, parrecflag, dcm_name)
% 
% combine multi-echo GRE data (NIFTI from dcm2niix) and read Params from json file and save to header.mat 
% 
% nifti_dirs: folders with nifti output from dcm2niix
% cleanup   : cleanup the original NIFTI and json files from dcm2niix
% parrecflag: flag for par/rec data
% dcm_name  : original dcm folder or eDICOM file; for reading eDICOM Params
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu

% 2024-03-21, X.L., for SIEMENS data with diff series number for mag vs. phase
% 2024-04-23, X.L., bug fix for naming convention prefix

if nargin < 2
    cleanup = 0;        % default no cleanup
    parrecflag = 0;
    dcm_name = [];
elseif nargin < 3
    parrecflag = 0;
    dcm_name = [];
elseif nargin < 4
    dcm_name = [];
end

% if multiple folders
if ~iscell(nifti_dirs)
    nifti_dirs = {nifti_dirs};
end

%% echo combinatio
for nifti_ii = 1:length(nifti_dirs)

    nifti_dir = nifti_dirs{nifti_ii};

    disp(nifti_dir);

    % combine nifti file from each echo
    json_test = dir(fullfile(nifti_dir, '*_ph.json'));

    if length(json_test) > 1
        % if with multi-echoes
      
        
        json_list = dir(fullfile(nifti_dir, '*_e1.json'));            % magnitude echo 1
        if ~isempty(json_list) %% philisps data
            filename_prefix_mag = extractBefore(json_list(1).name, '_e1.json');     % prefix without "_" now
        else %% SIEMENS data
            json_list_ph = dir(fullfile(nifti_dir, '*_e*_ph.json'));
            echotimes=zeros(1,length(json_list_ph));
            for i = 1:length(json_list_ph)
                echotimes(i) = extractNumberFromFilename(json_list_ph(i));
                % fprintf('文件 "%s" 中的数字: %f\n', filenames{i}, num);
            end
            [~, idx] = sort(echotimes);
            echotimes_sorted = echotimes(idx);
            json_list_ph_sorted = json_list_ph(idx);

            % 重命名文件
            fprintf('开始重命名文件...\n');
            % 获取文件夹中的所有文件
            files = dir(nifti_dir);
            % 存储结果
            results = [];

            for i = 1:length(files)
                fileName = files(i).name;
                
                % 跳过目录
                if files(i).isdir
                    continue;
                end
                
                % 使用正则表达式提取文件名中的 echotime 数值
                pattern = 'e(\d+\.?\d*)';
                tokens = regexp(fileName, pattern, 'tokens');
                
                if ~isempty(tokens)
                    echotime_str = tokens{1}{1};
                    echotime_value = str2double(echotime_str);
                    
                    % 查找匹配的索引
                    [found, index] = ismember(echotime_value, echotimes_sorted);
                    
                    if found
                        % 构造新的文件名，将 echotime 替换为索引
                        newFileName = regexprep(fileName, pattern, ['e' num2str(index)]);
                        
                        filename_prefix = extractBefore(json_list_ph(1).name, '_(MSV21)_'); 
                        filename_postfix = extractAfter(newFileName, ['e' num2str(index)]);
                        newFileName=[filename_prefix '_e' num2str(index)  filename_postfix];
                        % 重命名文件
                        oldFilePath = fullfile(nifti_dir, fileName);
                        newFilePath = fullfile(nifti_dir, newFileName);
                        movefile(oldFilePath, newFilePath);
                        
                        fprintf('已重命名: %s -> %s (索引: %d)\n', fileName, newFileName, index);
                        
                        % 存储结果
                        results = [results; struct('originalName', fileName, 'newName', newFileName, 'echotime', echotime_value, 'index', index)];
                    else
                        fprintf('文件: %s, Echotime: %.2f, 在 echotimes_sorted 中未找到匹配\n', fileName, echotime_value);
                    end
                end
            end
            
            
        end
    
        json_list_ph = dir(fullfile(nifti_dir, '*_e*_ph.json')); 
        num_echo = length(json_list_ph);
        filename_prefix_phase = extractBefore(json_list_ph(1).name, '_e1_ph.json');     % prefix without "_" now
        filename_prefix_mag=filename_prefix_phase;
    elseif length(json_test) == 1
        % if with single-echo OR with eDICOM data
        num_echo = 1;
        json_list_ph = json_test;
        filename_prefix_phase = extractBefore(json_list_ph(1).name, '_ph.json');  % prefix without "_"

        json_list = dir(fullfile(nifti_dir, '*.json'));
        for json_ii = 1:length(json_list)
            if ~contains(json_list(json_ii).name, json_list_ph(1).name)
                filename_prefix_mag = extractBefore(json_list(json_ii).name, '.json');
            end
        end

    else
        error('There is no phase data.')
    end

    % magnitude & phase
    img_mag     = [];
    img_phase   = [];

    for kecho = 1:num_echo
        if num_echo > 1
            curr_echo_mag_filename      = strcat(filename_prefix_mag, '_e', num2str(kecho),'.nii.gz');
            curr_echo_phase_filename    = strcat(filename_prefix_phase, '_e', num2str(kecho),'_ph.nii.gz');
        else
            curr_echo_mag_filename      = strcat(filename_prefix_mag, '.nii.gz'); % prefix without "_"
            curr_echo_phase_filename    = strcat(filename_prefix_phase, '_ph.nii.gz');
        end

        img_mag     = cat(4,img_mag, load_nii_img_only(fullfile(nifti_dir,curr_echo_mag_filename)));
        img_phase   = cat(4,img_phase, load_nii_img_only(fullfile(nifti_dir,curr_echo_phase_filename)));
        
    end

    % check phase scaling 
    phase_range = max(img_phase(:)) - min(img_phase(:));
    if phase_range > 2*pi + 100*eps
        disp('correct phase scaling ...')
        img_phase = img_phase./phase_range*(2*pi);
    end

    % check num_echo in case of eDICOM data
    if size(img_phase, 4) > num_echo
        flag_eDICOM = 1;
    else
        flag_eDICOM = 0;    % default
    end

    % find longest common prefix & save
    for prefix_i = 1:length(filename_prefix_mag)
        prefix_pat = filename_prefix_mag(1:prefix_i);
            if contains(filename_prefix_phase, prefix_pat)
                filename_prefix = prefix_pat;
            else
                continue;
            end
    end
    % intersect will remove repetitions
    % filename_prefix = intersect(filename_prefix_mag, filename_prefix_phase, 'stable'); 

    output_mag_filename     = strcat(filename_prefix, '_GRE_mag.nii.gz');
    output_phase_filename   = strcat(filename_prefix, '_GRE_phase.nii.gz');

    save_nii_img_only(fullfile(nifti_dir,curr_echo_mag_filename),fullfile(nifti_dir,output_mag_filename),img_mag);
    save_nii_img_only(fullfile(nifti_dir,curr_echo_phase_filename),fullfile(nifti_dir,output_phase_filename),img_phase);

    disp(['GRE mag & phase data saved at ', nifti_dir, '.'])

    % Extract Params from json files
    nii_phase = load_untouch_nii(fullfile(nifti_dir,output_phase_filename));

    % Params only for NIFTI output, thus should be single volume
    Params.nifti_hdr = nii_phase.hdr;   % nifti head for output, with multi-layer
    Params.nifti_hdr.dime.dim(5) = 1;   % should be echo combined
    Params.nifti_hdr.dime.pixdim(5) = 0;
    Params.nifti_hdr.dime.dim(1) = 3;

    if flag_eDICOM > 0
        % read Params.TEs/B0/TR from eDICOM
        [~, dicomheader] = dicomeread(dcm_name);
        Params_temp = [];
        Params_temp = readparamsfromdicom(dicomheader, Params_temp);
        Params.TEs = Params_temp.TEs;
        Params.TR = Params_temp.TR;
        Params.B0 = Params_temp.B0;

    else
        % read TE from json file
        json_list = dir(fullfile(nifti_dir, [filename_prefix_phase, '*_ph.json']));
        json_filenames = cell(length(json_list), 1);
        for json_ii = 1:length(json_list)
            json_filenames{json_ii} = fullfile(json_list(json_ii).folder, json_list(json_ii).name);
        end
        Params.TEs = readTE_dcm2niix_JSON(json_filenames);

        % read other params from json file, B0, TR etc.
        Params = readParams_dcm2niix_JSON(json_filenames{1}, Params);
    end

    Params.nEchoes = length(Params.TEs);
    Params.sizeVol = nii_phase.hdr.dime.dim(2:4);  % 2-4, pixdim, 5:echoes, 6:dynamics?
    Params.voxSize = nii_phase.hdr.dime.pixdim(2:4);
    Params.fov = Params.sizeVol.*Params.voxSize;
    Params.nDynamics = nii_phase.hdr.dime.dim(6);  % need check 

    [Params.b0dir, Params.TAng] = get_B0_dir_from_nifti(nii_phase);
    % in case we got complex number
    Params.TAng = real(Params.TAng);
    Params.b0dir = real(Params.b0dir);

    % par/rec to NIFTI test, convert to LAS, b0 needs correction
    if parrecflag == 1
       disp('usning NIFTI file converted from par/rec, check with caution...')
       negz = diag([1, 1, -1]);   %
       Params.b0dir = negz*Params.b0dir; Params.TAng = Params.TAng*negz;
    end

    if flag_eDICOM > 0
       disp('usning NIFTI file converted from enhanced DICOM, check with caution...')
       negxz = diag([-1, 1, -1]);   %
       Params.b0dir = negxz*Params.b0dir; Params.TAng = Params.TAng*negxz;
    end

    % save header .mat file
    output_header_filename     = strcat(filename_prefix, '_header.mat');
    save(fullfile(nifti_dir, output_header_filename), "Params", '-v7.3');

    disp('GRE header .mat file saved.')

    % clean up if selected
    if cleanup
        if num_echo > 1
            % if multi-echo
            delete(fullfile(nifti_dir, [filename_prefix_mag, '_e*.json']));
            delete(fullfile(nifti_dir, [filename_prefix_mag, '_e*.nii.gz']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_e*.json']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_e*.nii.gz']));
        else
            % if single-echo
            delete(fullfile(nifti_dir, [filename_prefix_mag, '.json']));
            delete(fullfile(nifti_dir, [filename_prefix_mag, '.nii.gz']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_ph.json']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_ph.nii.gz']));
        end

    end

    end
end


function echo_time = extractNumberFromFilename(filename)
    % 提取 e 和 .json 之间的数字（处理可能的 '_ph' 后缀）
    pattern = 'e(\d+(?:\.\d+)?)(?:_ph)?\.json';
    tokens = regexp(filename.name, pattern, 'tokens');
    
    if ~isempty(tokens)
        echo_time = str2double(tokens{1}{1});
    else
        echo_time = NaN;  % 未找到匹配时返回 NaN
    end
end

function newFilename = replaceNumberInFilename(filename, newNumber)
    % 将文件名中的数字替换为新数字
    % 例如: Axial_3D_e6.71_ph.json → Axial_3D_e1_ph.json
    
    % 找到 e 和 .json 的位置
    eIdx = strfind(filename.name, 'e');
    jsonIdx = strfind(filename.name, '.json');
    
    if ~isempty(eIdx) && ~isempty(jsonIdx) && eIdx < jsonIdx
        % 提取 e 之前和 .json 之后的部分
        prefix = filename.name(1:eIdx);
        suffix = filename.name(jsonIdx:end);
        
        % 构建新文件名
        newFilename = [prefix, num2str(newNumber), suffix];
    else
        % 如果无法处理，保持原文件名
        newFilename = filename;
        warning('无法处理文件名: %s', filename);
    end
end