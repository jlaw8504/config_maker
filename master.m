cd('/home1/lawrimor/add_condensin')
add_condensin_ver3(1, 'pT431_equilibrium.out', 'pt431_s1_0000.out', 3)
system('/home1/lawrimor/source/brownianMotion/chromoShake/chromoShake -openCL_dir $HOME/source/brownianMotion/chromoShake/ -save $HOME/add_condensin/pt431_s1_0000.out 5000 10 -continue')
for n = 1:1000
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

	condensin_step_ver3(strcat('pt431_s1_',str_padi,num2str(i),'.out'),...
        strcat('pt431_s1_',str_pad,num2str(n),'.out'))
	system(strcat('/home1/lawrimor/source/brownianMotion/chromoShake/chromoShake -openCL_dir $HOME/source/brownianMotion/chromoShake/ -save $HOME/add_condensin/pt431_s1_',str_pad,num2str(n),'.out 5000 10 -continue'))
end
