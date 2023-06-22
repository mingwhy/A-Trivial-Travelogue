library(infotheo)

entropy(c(1,0,0,0,0))
entropy(c(0,1,1,1,1))
-1*(0.2*log(0.2)+0.8*log(0.8))


sex.label=c(rep(0,5),rep(1,5));
sex.label

# H(S)
H_naive=entropy(sex.label)

gene1=c(1,1,1,0,0,1,0,0,0,0)
# H(S|G)=H(S,G)/H(G)
H_g=entropy(gene1)

H_gs=sum(tapply(sex.label,gene1,entropy)*table(gene1)/length(gene1))
H_gs #0.606


# mutual information: I(S;G)=H(S)-H(S|G)
H_naive - H_gs 
# 0.0863 #information gain

