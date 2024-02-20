include("graph.jl")
include("core.jl")
include("bb.jl")


using LinearAlgebra
using Laplacians
using SparseArrays

fname = open("filename.txt", "r")
str   = readline(fname);
nn     = parse(Int, str);


for nnnn=1:nn 
    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
    on=G.n;om=G.m;
    Gc=findconnect(G)
    G=Gc;
    n=G.n;m=G.m;
    L=lapsp(G);
    A=adjsp(G);
    B=getB(G); 
    eps=0.01;
    t=Int(round(log(n)/eps^2/5)+1);
    Z=zeros(t,n);
    f=approxchol_lap(A);
    for i=1:t
         Q=randn(1,m);
         QB=Q*B;
         Z[i,:]=f(QB[1,:]);
    end

    new_lab=Z';
    dis1=zeros(n);

    for i=1:n
        for j=1:n
            tmp_dis=sum((new_lab[i,:] .- new_lab[j,:]).^2)/t;
            if dis1[i] < tmp_dis
                dis1[i]=tmp_dis;
            end
        end
    end


end
close(fname)
