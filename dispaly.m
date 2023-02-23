
    for j=1:9
        pathname=strcat(num2str(j),'.mat');
        ccdgray=load num2str(j).mat;
        I(:,:,jn)=ccdgray(m1:m2,n1:n2);
        jn=jn+1;
    end;