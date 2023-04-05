using Random,Statistics,LinearAlgebra,Distributions,DataFrames,DelimitedFiles;
#data_path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/MTSBayesC/data/"
Random.seed!(123)
geno_path = ARGS[1]
data_path = ARGS[2]
geno=readdlm(geno_path * "gene_QC_n503_p18588.txt")
n,p=size(geno)

nQTL=Int(floor(p*0.01))
QTLpos=sample(MersenneTwister(0),collect(1:p),nQTL,replace=false, ordered=true)
@show nQTL
writedlm(data_path * "QTLpos.txt",QTLpos)

################################################################################
#Step1. simulate marker effects
################################################################################
μ=[0,0]
Σ=[1     0.54
   0.54  1]           #covariance of marker effects
writedlm(data_path * "QTLEffCov.txt",Σ)
betaAll=rand(MersenneTwister(1), MvNormal(μ,Σ), nQTL)
beta1=betaAll[1,:]
beta2=betaAll[2,:]
cov(beta1,beta2)
writedlm(data_path * "QTLeff1.txt",beta1)
writedlm(data_path * "QTLeff2.txt",beta2)
@show cov(beta1,beta2)

h21=0.7
h22=0.3
bv1=geno[:,QTLpos]*beta1
bv1=bv1/std(bv1)*sqrt(h21)
bv2=geno[:,QTLpos]*beta2
bv2=bv2/std(bv2)*sqrt(h22)
println("varg1 = ", var(bv1))
println("varg2 = ", var(bv2))
println("cov(bv1,bv2) = ", cov(bv1,bv2))
println("cor(bv1,bv2) = ", cor(bv1,bv2))
writedlm(data_path * "bv1.txt",bv1)
writedlm(data_path * "bv2.txt",bv2)

################################################################################
#Step2. simulate residuals
################################################################################
μ=[0,0]
Σ=[0.3   0.235
   0.235 0.8] # covariance of residual effects
writedlm(data_path * "residEffCov.txt",Σ)
eAll=rand(MersenneTwister(1), MvNormal(μ,Σ), n)
e1=eAll[1,:]
e2=eAll[2,:]
e1=e1/std(e1)*sqrt(1-h21)
e2=e2/std(e2)*sqrt(1-h22)
println("vare1 = ", var(e1))
println("vare2 = ", var(e2))
println("cov(e1,e2) = ", cov(e1,e2))
println("cor(e1,e2) = ", cor(e1,e2))

y1=bv1+e1
y1=(y1.-mean(y1))/std(y1)
y2=bv2+e2
y2=(y2.-mean(y2))/std(y2)
writedlm(data_path * "y1.txt",y1)
writedlm(data_path * "y2.txt",y2)

#center and scale genotypes
X_scale=deepcopy(geno)
for j in 1:p
	X_scale[:,j] = (geno[:,j] .- mean(geno[:,j]))/std(geno[:,j])
end

######### Summary Statistics ########
D=diagm(repeat([n-1],p))*1.0
b1=inv(D)*X_scale'y1                          # bhat1
b2=inv(D)*X_scale'y2                          # bhat2
writedlm(data_path * "bhat1.txt",b1)
writedlm(data_path * "bhat2.txt",b2)

B=inv(sqrt(D))*X_scale'X_scale*inv(sqrt(D));  # LD matrix 
writedlm(data_path * "LDmatrix.txt",B)


# eigen-decomposition of B matrix (selected eigen-values) 
eigen_values = eigvals(B)
selected = findall(x -> x > 1e-3, eigen_values)
lambda = eigen_values[selected]
eigen_vectors = eigvecs(B)
U = eigen_vectors[:,selected];   # Uq of dimension n_snp x n_eigen
writedlm(data_path * "eigen_values.txt",lambda)
writedlm(data_path * "eigen_vectors.txt",U)



#=
# eigen-decomposition of B matrix (complete)
eigen_values = eigvals(B)
eigen_vectors = eigvecs(B)

writedlm(data_path * "eigen_values_complete.txt",eigen_values)
writedlm(data_path * "eigen_vectors_complete.txt",eigen_vectors)
@show length(eigen_values)
@show size(eigen_vectors)
=#