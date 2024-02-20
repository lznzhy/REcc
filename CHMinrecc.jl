include("graph.jl")
include("core.jl")
include("meth.jl")
include("bb.jl")
include("bb2.jl")


using LinearAlgebra
using Laplacians
using SparseArrays



nnn=1;

for ttttt=1:nnn

        fname = open("filename.txt", "r")
        str   = readline(fname);
        nn     = parse(Int, str);

        str = readline(fname);
        str = split(str);
        G   = get_graph(str[1]);
        on=G.n;om=G.m;
        Gc=findconnect(G)
        G=Gc;
        n=G.n;m=G.m;

        eps=0.1;
        t=Int(round(log(n)/eps^2/5)+1);

        L=lapsp(G);
        A=adjsp(G);
        B=getB(G); 
        Q=randn(t,m); 
        QB=Q*B; 
        Z=zeros(t,n);
        f=approxchol_lap(A);
        for i=1:t
            Z[i,:]=f(QB[i,:]);
        end
        new_lab=Z';
        K=500;
        index_this=bb(new_lab,K);
        ll=length(index_this);
    
        s=105;
        kkk=50;
        for i=1:kkk
            disfi2=0;
            chdis=0;
            chxx=0;chyy=0;
            nodedis=new_lab;
            for ch1 in index_this
                for ch2 in index_this
                    tmp_disch=sum((nodedis[ch1,:] .- nodedis[ch2,:]).^2)/t;
                    if chdis < tmp_disch
                        chdis=tmp_disch;
                        chxx=ch1;chyy=ch2;
                    end
                end
            end

            Gv=G.v;
            push!(Gv,chxx);
            Gu=G.u;
            push!(Gu,chyy);
            G=Graph(G.n, G.m+1, Gv, Gu,G.nbr);
            n=G.n;m=G.m;
            if i%5==0
                A=adjsp(G);
                disfi2=accRecc(A,n,s);
                println(disfi2);
            end
            t=150;
            L=lapsp(G);
            A=adjsp(G);
            B=getB(G); 
            Q=randn(t,m); 
            QB=Q*B; 
            Z=zeros(t,n);
            f=approxchol_lap(A);
            for i=1:t
                Z[i,:]=f(QB[i,:]);
            end
            new_lab=Z'; 
            K=500;
            keach=0.05;
            keach2=18-keach*i;
            index_this=bb2(new_lab,K,keach2);
            ll=length(index_this);

        end
        

    end
 

