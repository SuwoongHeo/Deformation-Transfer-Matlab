close all;

addpath('./1. BasicModules/kd_tree');
addpath('./1. BasicModules');
addpath('./2. Non_rigid_registration');
%  distcomp.feature( 'LocalUseMpiexec', false )
%% Example
% if ~exist('VS','var')
%     [VS FS NS] = read_obj_file('./3. Data/face-poses/face-reference.obj');
%     [VS2 FS2 NS2] = read_obj_file('./3. Data/face-poses/face-09-surprise.obj');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %     [VT FT NT] = read_obj_file('./3. Data/head-poses/head-reference.obj');
%     [VT FT NT] = read_obj_file('./3. Data/face-sw-morph/face-sw-reference.obj');
% end

% marker = init_marker(VS, FS, VT, FT, 'Face_Marker_sw.mat');
%% Using ICIP paper results
if ~exist('VS','var')
    [VS FS NS] = read_obj_file('./3. Data/FW-FACS/FW_au00.obj');
    [VS2 FS2 NS2] = read_obj_file('./3. Data/FW-FACS/FW_au01.obj');
%     [VS FS NS] = read_obj_file('./3. Data/ICIP-recon/Udara_res_adap.obj');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     [VT FT NT] = read_obj_file('./3. Data/head-poses/head-reference.obj');
%     [VT FT NT] = read_obj_file('./3. Data/FW-FACS/FW_au00.obj');
%     [VT FT NT] = read_obj_file('./3. Data/ICIP-recon/Udara_res_adap.obj');
    load('MM.mat');
    VT = MM.vertices; FT = MM.faces2; NT = MM.normals;
end

% load './3. Data/ICIP-recon/MM_lands_FW.mat'
% load './3. Data/FW-FACS/FWland.mat'
load 'MM_FWland.mat';
load 'MM_land.mat';
% 
marker = [FWland landmark(:,1)];
% marker = [landmark FWland];

% [ VS_Reg VT_Reg ] = non_rigid_registration(VS, FS, VT, FT, 1.0, 0.01, [1 500 3000 5000], marker, 'DF_reg_phase2.mat');
[ VT_Reg VS_Reg ] = non_rigid_registration(VT, FT, VS, FS, 10, 0.001, [1 500 3000 5000], [marker(:,2) marker(:,1)], 'DF_reg_phase2.mat');
corres = build_correspondence(VS_Reg, FS, VT_Reg, FT, 10, 0.05, 'Face_ICIP_corres.mat');
[ x nx ] = deformation_transfer(VS, FS, VT, FT, VS2, FS2, corres);
% write_obj_file(x, FS, nx, './4. Result/objfile/head-03-fury(result).obj');
write_obj_file(x, FS, nx, './3. Data/face-sw-morph/face-sw-09-surprise.obj');

fprintf('End of demo..\n');
system('pause');

clear VS VT S_factor T_factor FS FT NS NT maker;