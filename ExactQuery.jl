include("graph.jl")
include("core.jl")

using LinearAlgebra
using Laplacians
using SparseArrays

####用凸包求所有点的有效电阻偏心率
fname = open("filename.txt", "r")
str   = readline(fname);
nn     = parse(Int, str);


for nnnn=1:nn #nn个网络
    #######读入数据
    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
    on=G.n;om=G.m;
    Gc=findconnect(G)
    G=Gc;
    n=G.n;m=G.m;
    L=lapsp(G);
    A=adjsp(G);
    J=ones(n,n)./n;
    invL=inv(L+J)-J;

    dis1=zeros(n);
    for i=1:n
        for j=1:n
            tmp_dis=invL[i,i]+invL[j,j]-2*invL[i,j]
            if dis1[i] < tmp_dis
                dis1[i]=tmp_dis;
            end
        end
   end
    
   filename = str[1] * "output.txt"
    open(filename, "w") do file
        for i in 1:n
            println(file, dis1[i])
        end
    end

end
close(fname)
