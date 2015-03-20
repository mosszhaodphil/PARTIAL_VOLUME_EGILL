% This function performs partial volume correction on M0 tissue
% Input:
% m0_file: averaged control image obtained from quasil
% Output 
% m0t_[tissue_type]_[kernel_size].nii.gz: corrected M0t image


function csf_correct(m0t_file, csf_file, gm_file, wm_file, mask_file, kernel_size)


    % Output file names:
    file_name_m0t_csf = strcat('M0t_csf_k', num2str(kernel_size));
    file_name_m0t_gm  = strcat('M0t_gm_k', num2str(kernel_size));
    file_name_m0t_wm  = strcat('M0t_wm_k', num2str(kernel_size));
    
    % Load data files
    [data,dims,scales] = ra(m0t_file);
    csf                = ra(csf_file);
    gm                 = ra(gm_file);
    wm                 = ra(wm_file);
    mask               = ra(mask_file);

    % Concatinate along the time dimension
    pve = cat(4, csf, gm, wm);

    % Create two empty matrices to save results
    CSFdata = zeros(size(data));
    GMdata  = zeros(size(data));
    WMdata  = zeros(size(data));

    % Check kernel size
    if (length(kernel_size) == 1)
        nsel = kernel_size;nzsel=kernel_size;

    elseif (length(kernel_size) == 2)
        nsel = kernel_size(1);
        nzsel = kernel_size(2);

    else
        error('Size must be scalar or have two entries (xy dimension and z dimension)');
    end

    % Get the dimension of mask
    xsize = size(mask,1);
    ysize = size(mask,2);
    zsize = size(mask,3);

    display('Performing partial volume correction...');

    count=1; % Used to monitor number of iterations

    for i = 1 : size(mask,1)
        for j = 1 : size(mask,2)
            for k = 1 : size(mask,3)

                % Only work on voxels with positive intensity
                if mask(i, j, k) > 0
                    % Create a submask of current kernel size
                    submask = mask(max(i - nsel, 1) : min(i + nsel, xsize), max(j - nsel, 1) : min(j + nsel, ysize), max(k - nzsel, 1) : min(k + nzsel, zsize));
                    
                    % calculate the sum of all elements in submask
                    % proceed if sum is greater than 5 (arbitrary threshold)
                    if (sum(sum(sum(submask))) > 5)
                        % Extract sub perfusion and PV map within the current kernel
                        subdata = vols2matrix(data(max(i - nsel, 1) : min(i + nsel, xsize), max(j - nsel, 1) : min(j + nsel, ysize), max(k - nzsel, 1) : min(k + nzsel, zsize), :), submask);
                        subpve  = vols2matrix(pve(max(i - nsel, 1) : min(i + nsel, xsize), max(j - nsel, 1) : min(j + nsel, ysize), max(k - nzsel, 1) : min(k + nzsel, zsize), :), submask);
                        
                        % Get pseudo inversion matrix
                        pveinv = pinv(subpve);
                        % Ger averaged PV value of current kernel
                        pveprop = sum(subpve) / size(subpve, 1);
                        
                        % Loop through all time points
                        for ti = 1 : size(data, 4)
                            % Calcuated corrected perfusion value for GM and WM
                            subd = pveinv * subdata(:, ti);
                            % Assign values to CSF, GM, and WM result matrix
                            CSFdata(i ,j ,k , ti) = subd(1);
                            GMdata(i ,j ,k , ti)  = subd(2);
                            WMdata(i ,j ,k , ti)  = subd(3);
                        end

                        % Deal with cases where there is very little GM or WM
                        % Within the sample volume
                        % Then assign zero to GM and WM perfusion map
                        if (pveprop(1) < 0.01)
                            GMdata(i, j, k, :) = 0;
                        end
                        if (pveprop(2) < 0.01)
                            WMdata(i, j, k, :) = 0;
                        end
                        if (pveprop(3) < 0.01)
                            WMdata(i, j, k, :) = 0;
                        end
                       
                    end
                end
            end

            count = count + 1;
            % Display iteration process every 100 times
            if (rem(count, 100) == 0)
                disp(count);
            end

        end
    end

    % Save results
    save_avw(CSFdata, file_name_m0t_csf, 'f', scales);
    save_avw(GMdata, file_name_m0t_gm, 'f', scales);
    save_avw(GMdata, file_name_m0t_wm, 'f', scales);

    % UAT Moss
    %gm_file_handle.img = GMdata;
    %wm_file_handle.img = WMdata;

    % UAT Moss
    %save_nii(gm_file_handle, strcat('perfusion_gm_k', num2str(kernel_size), '.nii.gz'));
    %save_nii(wm_file_handle, strcat('perfusion_wm_k', num2str(kernel_size), '.nii.gz'));

    display(['Results saved in ' file_name_m0t_csf ', ' file_name_m0t_gm ', and ' file_name_m0t_wm]);
    display('Finish');



end

