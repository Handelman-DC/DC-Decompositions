function y=lex_index_multi(monom,n_alpha,d_p2,N_sim,N_tot)
lex_indeces=zeros(N_sim,1);
y=0;
for i=1:N_sim
    if i==1
        index=1;
    else
        index=sum(n_alpha(1:i-1))+1;
    end
    index2=sum(n_alpha(1:i));
    lex_indeces(i)=lex_index(monom(index:index2),n_alpha(i),d_p2(i));
    lex_indeces(i)=N_tot(i)-lex_indeces(i)+1;
end
comb=ones(N_sim-1,1);
for i=1:N_sim-1
    for j=i:N_sim-1
        comb(i)=comb(i)*N_tot(j+1);
    end
end

for j=1:N_sim-1
    y=y+(lex_indeces(j)-1)*comb(j);
end
y=y+lex_indeces(end);
end


function I=lex_index(alfa,n,d)
I=0;
for j=1:n-1
    J=0;
    s=0;
    for k=1:j-1
        s=s+alfa(k);
    end
    for i=1:alfa(j)
        J=J+factorial(n-j+d+1-s-i-1)/(factorial(d+1-s-i)*factorial(n-j-1));
    end
    I=I+J;
end
I=I+1;
end
