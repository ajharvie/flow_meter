function [fEst, phEst, intM] = Simulator(y1,y2,timestep,windowsize)
    freqspace = (1/timestep)*(1/1500);
    a = zeros(1,windowsize);
    twiddle = complex(a,0);
    comb = zeros(1, windowsize);
    kEst = zeros(1, length(y1)-windowsize);
    fEst = zeros(1, length(y1)-windowsize);
    intM = zeros(1, length(y1)-windowsize);
    phEst = zeros(1, length(y1)-windowsize);
    filtered_dft1 = zeros(1, windowsize);
    filtered_dft2 = zeros(1, windowsize);
    %phDiff = zeros(1, length(y1));


    %circular buffers for time domain values
    live_samples1 = zeros(windowsize,1);
    dft1 = zeros(windowsize,1);
    
    live_samples2 = zeros(windowsize,1);
    dft2 = zeros(windowsize,1);
    
    sample_index = 1;
    r = 0.99999999999999; %damping factor
    rton = power(r,windowsize);

    %compute twiddle factors
    for k = 1:windowsize
        factor = (2*pi*(k-1))/windowsize;
        twiddle(k) = exp(1i*factor);
    end
    
    for k = 1:length(y1)

        %dft time domain values
        old_1 = live_samples1(sample_index);
        live_samples1(sample_index) = y1(k);
        old_2 = live_samples2(sample_index);
        live_samples2(sample_index) = y2(k);


        %update
        for i = 1:windowsize
            dft1(i) = twiddle(i) * r * (dft1(i) - rton*old_1 + y1(k));
            dft2(i) = twiddle(i) * r * (dft2(i) - rton*old_2 + y2(k));
        end

        %next element in circular buffer
        sample_index = sample_index + 1;
        if sample_index>windowsize
            sample_index=1;
        end    

        dft_out1=transpose(dft1);
        dft_out2=transpose(dft2);

        %filtering

        filtered_dft1(1) = -0.25*(dft_out1(length(dft_out1)) + dft_out1(length(dft_out1))) + 0.5*dft_out1(1);
        filtered_dft1(length(dft_out1)) = -0.25*(dft_out1(1) + dft_out1(length(dft_out1)-1)) + 0.5*dft_out1(length(dft_out1));
        for i = 2:length(dft_out1)-1
            filtered_dft1(i) = -0.25*(dft_out1(i-1) + dft_out1(i+1)) + 0.5*dft_out1(i);
        end

        filtered_dft2(1) = -0.25*(dft_out2(length(dft_out2)) + dft_out2(length(dft_out2))) + 0.5*dft_out2(1);
        filtered_dft2(length(dft_out2)) = -0.25*(dft_out2(1) + dft_out2(length(dft_out2)-1)) + 0.5*dft_out2(length(dft_out2));
        for i = 2:length(dft_out2)-1
            filtered_dft2(i) = -0.25*(dft_out2(i-1) + dft_out2(i+1)) + 0.5*dft_out2(i);
        end
        
        ph1 = angle(filtered_dft1);
        ph2 = angle(filtered_dft2);
        phdiff = ph2-ph1;
        
        for j = 1:length(phdiff)
            if phdiff(j) < 0
                phdiff(j) = phdiff(j) + 2*pi;
            end
        end
        
        
        R1 = abs(filtered_dft1);
        R2 = abs(filtered_dft2);
        
        if k>windowsize
            comb = sqrt(R1.*R2);
        
            intM(k-windowsize) = sum(comb(1:200));

            [~,I] = max(comb(3:round(windowsize/2))); %ignoring DC
            I = I + 2;
            delta = 2* (comb(I+1) - comb(I-1))/(2*comb(I) + comb(I-1) + comb(I+1));
            kEst(k-windowsize) = I - 1 + delta;  
            fEst(k-windowsize) = freqspace*kEst(k-windowsize);
            
            phEst(k-windowsize) = interp1(1:windowsize,phdiff,kEst(k-windowsize),'linear');
            
            if phEst(k-windowsize)<0
                phEst(k-windowsize) = phEst(k-windowsize) + 2*pi;
            end
            
        end



    end
       %%crunch


    
end





