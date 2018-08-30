clear

addpath('/data/p_nmc002/matlab_tools/niftitools')
addpath('/data/p_nmc002/matlab_tools/spinconv')

%fpath = '/data/pt_nmc002/ROIs/frequency_localizer/max_beta_per_subject_zscored/';
fpath = '/data/pt_nmc002/ROIs/frequency_localizer/max_beta_per_subject_zscored_normed/warped/';

maps = {
'sub-01_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii'
'sub-02_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-03_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-04_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-05_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-06_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-09_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-10_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-11_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-12_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-13_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-15_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-16_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-18_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-20_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-22_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-24_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-28_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-29_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-30_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-31_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-32_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-35_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-37_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-39_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-40_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-41_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii' 
'sub-43_beta_max_F1-10_R_MGB_corr_none_thr_0.05_masked_r5_zscored_pred_trans.nii'
};

% Convert pitch (around y-axis) rotation from Euler angles to rotation
% matrices.
t_d = [0, -45, 0];
Ry = SpinConv('EA321toDCM', t_d);
Ry = [Ry; 0 0 0];
Ry = [Ry(1,:) 0; Ry(2,:) 0; Ry(3,:) 0; Ry(4,:) 1];
tform = affine3d(Ry);

slice_step = 1;  % Step used to sample slices.
plot_slices = 0;

for imap = 1:numel(maps)
    % Load map. load_nii uses radiological orientation, which means that we
    % need to left-right flip the images.
    freq_map = load_nii([fpath, maps{imap}]);  
    % Flip the rows.
    freq_map.img = flipdim(double(freq_map.img),1);
    % Rotate -45 deg around
    freq_map.imgRy = imwarp(freq_map.img,tform);
    freq_map.imgRy(freq_map.imgRy==0) = NaN;
    [i1, i2, i3] = ind2sub(size(freq_map.imgRy), find(~isnan(freq_map.imgRy)));
    slices = min(i3):slice_step:max(i3);    
    figc = 1;
    %figure
    angle = [];
    norms = [];
    
    for islice = 1:length(slices)
         
        slice = freq_map.imgRy(:, :, slices(islice));
        [Fx,Fy] = gradient(slice);
        Fx_r = reshape(Fx, size(Fx, 1) * size(Fx, 2), 1);
        Fy_r = reshape(Fy, size(Fy, 1) * size(Fy, 2), 1);
        indexvox = find(~isnan(Fx_r) & ~isnan(Fy_r));
        
        sliceAngles = [];
        slice_norms = [];
        if indexvox
            for i=indexvox'
                % Then you still need to adjust angles on the right side of the 
                % coordinate system, such that angles [0-360] deg:
                ang_norm = sqrt(Fy_r(i)^2+Fx_r(i)^2);
                if Fx_r(i) < 0
                    sliceAngles = [sliceAngles, 360 - acosd(Fy_r(i)/ang_norm)];  
                else
                    sliceAngles = [sliceAngles, acosd(Fy_r(i)/ang_norm)];  
                end
                slice_norms = [slice_norms, ang_norm];
            end
        end
        
        if plot_slices
            subplot(ceil(length(slices)/2), 2, islice)
            % This needs to be transposed otherwise it looks wrong. But
            % only for viewing! imagesc switches x-y axes.  
            imagesc(slice'); hold on
            tit = sprintf('slice %i', slices(islice));
            title(tit)
            % because transposed
            xlim([min(i1), max(i1)])            
            ylim([min(i2), max(i2)])
            set(gca,'YDir','normal')
            colormap summer
            quiver(Fx,Fy); hold off        
            figc = figc + 1;
        end
        
        angle = [angle, sliceAngles];
        norms = [norms, slice_norms];
    end
    %suptitle(maps{imap})
    
    % Put them all in one vector.
    angles{imap,:} = sort(angle);
    edges = 0:5:360;
    [a,b] = histc(angles{imap}, edges);
    for i = 1:length(a)
        weightedHist(imap, i) = sum(norms(b==i));
    end
    angles_binned(imap,:) = a;
%     figure, plot(angles{imap},y{imap}/length(y{imap}))
%     ptit = sprintf('Slices %s', maps{imap}(1:6));
%     title(ptit)
%     xlabel('Direction of increasing frequency')
%     ylabel('Percentage of voxels')
%     xlim([0, 360])
    
end
save('angles_per_sub_masked_r5_normed_right_unc0p05.mat', 'angles', 'angles_binned', ...
    'weightedHist')
disp('done')