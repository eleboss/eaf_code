%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is an unoptimized version, with square complexity. But it is still quite fast. We will release the optimized version with linear complexity once we have done the preparation.
% Author: Shijie Lin
% Email: lsj2048@connect.hku.hk
% License: MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
format longG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timescale = 1e6;
image_height = 260;
image_width = 346;
start_path = "..\EAD"; % your path to the EAD
spliter = ";"; % for win
% spliter = ":" % for linux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gt_path = fullfile(start_path, "gt_focus_points.txt");
fileID_score = fopen('result_egs_detail.txt','w'); % resutlt txt
allSubFolders = genpath(start_path);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, spliter);
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);
% Process all image files in those folders.
for k = 1 : numberOfFolders
	% Get this folder and print it out.
	thisFolder = listOfFolderNames{k};
	% Get MAT files.
	filePattern = sprintf('%s/*.mat', thisFolder);
	baseFileNames = dir(filePattern);
	numberOfImageFiles = length(baseFileNames);
	% Now we have a list of all files in this folder.
	if numberOfImageFiles >= 1
    thisFolder_sep=regexp(thisFolder,filesep,'split');
    sequence_name = thisFolder_sep{end};
		% Go through all those mat files.
      for f = 1 : numberOfImageFiles
        fullFileName = fullfile(thisFolder, baseFileNames(f).name);
        load(fullFileName);
        % load the focus position
        focus_position_file = fullfile(thisFolder, "focus_position.txt");
        fileID = fopen(focus_position_file);
        focus_position_ts_list = textscan(fileID, '%f %f');
        fclose(fileID);
        focus_position_list = focus_position_ts_list{2};
        focus_timestamp_list = focus_position_ts_list{1}./ timescale;
        % load event data 
        y_o = double(aedat.data.polarity.y);
        x_o = double(aedat.data.polarity.x);
        pol_o = double(aedat.data.polarity.polarity);
        focus_o = double(aedat.data.polarity.focusPosition);
        pol_o(pol_o==0) = -1;
        t_o = double(aedat.data.polarity.timeStamp) ./ timescale;
        t_o = t_o - t_o(1);
        % load Groundtruth
        gt_focus_file = gt_path;
        gt_focus_fileID = fopen(gt_focus_file);
        gt_focus = textscan(gt_focus_fileID,  '%s %f');
        gt_focus_seq_name_list = gt_focus{1};
        find_return = strfind(gt_focus_seq_name_list, sequence_name);
        index_gt_foc = 0;
        for i=1:length(find_return)
          if find_return{i} == 1
            index_gt_foc = i;
          end
        end
        gt_focus_pos_list = gt_focus{2};
        gt_focus_pos = gt_focus_pos_list(index_gt_foc);
        fclose(gt_focus_fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimal_focus_time = realmax;
optimal_focus_pos = realmax;
prev_optimal_focus_time = 0;
uptaded_focus_moving_step = realmax;
stopping_threshold = 0.001;
eventstart_time = t_o(1);
eventend_time = t_o(end);
golden_search_range = t_o(end) - t_o(1);
while golden_search_range > stopping_threshold
  % -------------------------------------------- %
  % load data in a smaller subset
  % -------------------------------------------- %
  x = x_o; y = y_o; pol = pol_o; t = t_o; focus = focus_o;
  idx = (t>=eventstart_time)&(t<=eventend_time);
  y(idx~=1)=[];
  x(idx~=1)=[];
  pol(idx~=1)=[];
  t(idx~=1)=[];
  focus(idx~=1)=[];
  % -------------------------------------------- %
  % computer ER to determine next best interval
  % -------------------------------------------- %
  event_rate = 0;
  event_rate_sec = 0;
  event_rate_list = [];
  event_rate_ts_list = [];
  time_range = t(end) - t(1);
  golden_search_range = time_range;
  delta_time = time_range * 0.618;
  event_size =  length(t);
  fst_start_time = t(1);
  second_start_time = t(1) + time_range * 0.381;
  CAPTURE = true;
  CAPTURE_SEC = true;
  EV_FLAG = true;
  EV_FLAG_SEC = true;
  corres_focus_pos_list = [];
  for i = 1:event_size
    ev_t = t(i);
    ev_x = x(i) + 1;
    ev_y = y(i) + 1;
    ev_p = pol(i);
    ev_focus = focus(i);
    % -------------------------------------------- %
    % the first interval
    % -------------------------------------------- %
    if ev_t > fst_start_time && ev_t <= (fst_start_time + delta_time) && EV_FLAG
      event_rate = event_rate + 1;
    elseif ev_t > fst_start_time + delta_time && EV_FLAG
      event_rate = event_rate / delta_time;
      event_rate_list = [event_rate_list, event_rate];
      event_rate_ts_list = [event_rate_ts_list, ev_t];
      event_rate = 0;
      fst_start_time = ev_t;
      EV_FLAG = false;
    end
    if ev_t >= (fst_start_time + delta_time/2) && CAPTURE
      CAPTURE = false;
      corres_focus_pos_list(end+1) = ev_focus;
    end
    % -------------------------------------------- %
    % the second interval
    % -------------------------------------------- %
    if ev_t > second_start_time && ev_t <= (second_start_time + delta_time) && EV_FLAG_SEC
      event_rate_sec = event_rate_sec + 1;
    elseif ev_t > second_start_time + delta_time  && EV_FLAG_SEC
      event_rate_sec = event_rate_sec / delta_time;
      event_rate_list = [event_rate_list, event_rate_sec];
      event_rate_ts_list = [event_rate_ts_list, ev_t];
      event_rate_sec = 0;
      second_start_time = ev_t;
      EV_FLAG_SEC = false;
    end
    if ev_t >= (second_start_time + delta_time/2) && CAPTURE_SEC
      CAPTURE_SEC = false;
      corres_focus_pos_list(end+1) = ev_focus;
    end
  end
  % find the optimal focus point
  [max_ev_rate, max_ev_rate_index] = max(event_rate_list);
  optimal_focus_time = event_rate_ts_list(max_ev_rate_index);
  optimal_focus_pos = corres_focus_pos_list(max_ev_rate_index);
  uptaded_focus_moving_step = abs(optimal_focus_time - prev_optimal_focus_time);
  prev_optimal_focus_time = optimal_focus_time;
  % give next best event interval
  eventstart_time = optimal_focus_time - delta_time;
  eventend_time = optimal_focus_time;
end
error_cal =  (optimal_focus_pos - gt_focus_pos);
% print the error and step size
fprintf('     Processing mat file %s\n', fullFileName);
disp(["optimal_f_pos:", optimal_focus_pos, "steps:", uptaded_focus_moving_step , "e_cal:", error_cal])
% record results in the text file
fprintf(fileID_score,'%s & %.1f\n',sequence_name, error_cal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clearvars aedat
      end
  end
end
fclose(fileID_score);