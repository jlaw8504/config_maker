seed = 42;
cond_num = 88;
cd('/home1/lawrimor/config_maker');
system(sprintf('/home1/lawrimor/source/brownianMotion/chromoShake/chromoShake -openCL_dir $HOME/source/brownianMotion/chromoShake/ -save $HOME/config_maker/kb30_s%d_c%d_0000.out 5000 10 $HOME/config_maker/kb30_s%d_c%d.cfg',seed, cond_num, seed, cond_num));
system(sprintf('grep spring kb30_s%d_c%d_0000.out > springs_kb30_s%d_c%d.txt', seed, cond_num, seed, cond_num));
n = 1;
while n < 1001
	i = n - 1;
	if n < 10
		str_pad = '000';
	elseif n < 100
		str_pad = '00';
	elseif n < 1000
		str_pad = '0';
	else
		str_pad = [];
	end
	if i < 10
		str_padi = '000';
	elseif i < 100
		str_padi = '00';
	elseif i < 1000
		str_padi = '0';
	else
		str_padi = [];
	end

	condensin_step_ver3(sprintf('kb30_s%d_c%d_%s%d.out', seed, cond_num, str_padi, i),...
        sprintf('kb30_s%d_c%d_%s%d.out', seed, cond_num, str_pad, n));
	system(sprintf('/home1/lawrimor/source/brownianMotion/chromoShake/chromoShake -openCL_dir $HOME/source/brownianMotion/chromoShake/ -save $HOME/config_maker/kb30_s%d_c%d_%s%d.out 5000 10 -continue', seed, cond_num, str_pad, n));
	command = sprintf('grep nan kb30_s%d_c%d_%s%d.out | wc -l', seed, cond_num, str_pad, n);
	[~,result] = system(command);
	result_list = strsplit(result);
	nan_num = str2double(result_list{end-1});
	disp(strcat('Number of NaNs:',num2str(nan_num)));
	if nan_num > 1
		%remove the offending simulation iteration
		system(sprintf('rm kb30_s%d_c%d_%s%d.out', seed, cond_num, str_pad, n));
		[~,rseed_string] = system(sprintf('grep seed kb30_s%d_c%d_%s%d.out', seed, cond_num, str_padi, i));
		rseed_cell = strsplit(rseed_string);
		rseed = str2double(rseed_cell{end-1});
		new_rseed = rseed + 42;
		system(sprintf('sed -i ''s/random_number_seed %d/random_number_seed %d/g'' kb30_s%d_c%d_%s%d.out', rseed, new_rseed, seed, cond_num, str_padi, i))
	else
		n = n + 1;
		system(sprintf('grep spring kb30_s%d_c%d_%s%d.out >> springs_kb30_s%d_c%d.txt', seed, cond_num, str_padi, i, seed, cond_num));
		system(sprintf('rm kb30_s%d_c%d_%s%d.out', seed, cond_num, str_padi, i));
	end
end
