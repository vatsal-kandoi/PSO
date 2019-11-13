function cost = CostFunction(signal,z,l)
%% Function to calculate the SNR 
%%
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db13');
    [c,ll]=wavedec(z,l(1),Lo_D,Hi_D);
    A=wrcoef('a',c,ll,Lo_R,Hi_R,l(1));
    mod_sig=A;
    for i=1:l(1)
        D = wrcoef('d',c,ll,Lo_R,Hi_R,i);
        tD = wthresh(D,'s',l(2));
        mod_sig=mod_sig+tD;
    end
	cost = 20*log10(norm(signal(:)) / norm (signal(:)-mod_sig(:)));
end