clear

addpath('/data/p_nmc002/matlab_tools/niftitools')

fpath = '/data/pt_nmc002/ROIs/frequency_localizer/max_beta_per_subject_zscored_normed/';
%fpath = '/data/pt_nmc002/ROIs/frequency_localizer/max_beta_per_subjectnormed/';
%fpath = '/data/pt_nmc002/ROIs/frequency_localizer/';

maps = {
    'mean_normed_tonotopy_zscored_pred_left_45deg_resl.nii'
    'mean_normed_tonotopy_zscored_pred_right_45deg_resl.nii'
    %'mean_normed_files_added_-45deg_resl.nii'
    %'group_LMGB_z-beta_max_F1_10_masked_45deg_resl.nii'
    };
map_tit = {'left MGB', 'right MGB'};

% Convert pitch (around y-axis) rotation from Euler angles to rotation
% matrices.
t_d = [0, -45, 0];
Ry = SpinConv('EA321toDCM', t_d);
Ry = [Ry; 0 0 0];
Ry = [Ry(1,:) 0; Ry(2,:) 0; Ry(3,:) 0; Ry(4,:) 1];
tform = affine3d(Ry);gi

RefVec = [0, 1, 0];  % Reference Vector for angle calculation.
slice_step = 1;  % Step used to sample slices.

for imap = 1:numel(maps)
    % Load map
    freq_map = load_nii([fpath, maps{imap}]);   
    freq_map.img = double(freq_map.img);
    % Rotate -45 deg around
    %freq_map.imgRy = imwarp(freq_map.img,tform);
    freq_map.imgRy = freq_map.img;
    freq_map.imgRy(freq_map.imgRy==0) = NaN;
    % Find indices that actually have meaningful values and slice there.
    [i1, i2, i3] = ind2sub(size(freq_map.imgRy), find(~isnan(freq_map.imgRy)));
    slices = min(i3):slice_step:max(i3);    
    figc = 1;
    figure
    angle = [];
    
    for islice = 1:length(slices)
        
        slice = freq_map.imgRy(:, :, slices(islice))';
        slice(slice==0) = NaN;
        [Fx,Fy] = gradient(slice);
        Fx_r = reshape(Fx, size(Fx, 1) * size(Fx, 2), 1);
        Fy_r = reshape(Fy, size(Fy, 1) * size(Fy, 2), 1);
        indexvox = find(~isnan(Fx_r) & ~isnan(Fy_r));
        
        sliceAngles = [];
        if indexvox
            for i=indexvox'
                % Then you still need to adjust angles on the right side of the 
                % coordinate system, such that angles [0-360] deg:
                if Fx_r(i) < 0
                    sliceAngles = [sliceAngles, 360 - acosd(Fy_r(i)/sqrt(Fy_r(i)^2+Fx_r(i)^2))];  
                else
                    sliceAngles = [sliceAngles, acosd(Fy_r(i)/sqrt(Fy_r(i)^2+Fx_r(i)^2))];  
                end
            end
        end
        
        subplot(ceil(length(slices)/2), 2, islice)
        imagesc(slice); hold on
        tit = sprintf('slice %i', slices(islice));
        title(tit)
        % because slice is transposed
        xlim([min(i1), max(i1)])
        ylim([min(i2), max(i2)])
        set(gca,'YDir','normal')
        colormap summer
        quiver(Fx,Fy); hold off        
        figc = figc + 1;        
        angle = [angle, sliceAngles];
    end
    suptitle(map_tit{imap})
    
    % Put them all in one vector.
    angles{imap} = sort(angle);
    [a,b] = histc(angles{imap}, 0:10:360);
    y{imap} = a(b);
    figure, plot(angles{imap},y{imap}/length(y{imap}))
    ptit = sprintf('Slices in %s, average across subjects', map_tit{imap});
    title(ptit)
    xlabel('Direction of increasing frequency')
    ylabel('Percentage of voxels')
    xlim([0, 360])
    
end
save('angles.mat', 'angles', 'y')
disp('done')