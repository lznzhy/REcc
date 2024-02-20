include("graph.jl")
include("core.jl")
include("meth.jl")
include("bb2.jl")
include("bb.jl")

using LinearAlgebra
using Laplacians
using SparseArrays


for iii=1:1
    ####对单个点连上最远的k个点
    fname = open("filename.txt", "r")
    str   = readline(fname);
    nn     = parse(Int, str);
    #######读入数据
    str = readline(fname);
    str = split(str);
    G   = get_graph(str[1]);
    on=G.n;om=G.m;
    Gc=findconnect(G)
    G=Gc;
    
    s=100;
    allk=50;

    n=G.n;m=G.m;
    L=lapsp(G);
  

    #####降维
    eps=0.1;
    t=Int(round(log(n)/eps^2/10)+1);
    #t=150;
        
    A=adjsp(G);
    B=getB(G); #m*n
    Q=randn(t,m); #t*m
    QB=Q*B; #t*n
    Z=zeros(t,n);
    f=approxchol_lap(A);
    for i=1:t
        Z[i,:]=f(QB[i,:]);
    end
    new_lab=Z'; #####n*t n个点的新坐标

    K=500;
    kkkk=40;
    index_this=bb(new_lab,K);
    ll=length(index_this);
    println("凸包完成",ll);
    println("降维t  ",t)

    for numedg=1:allk 
        nodedis=new_lab;
        nn,mm=size(nodedis);
        farnode=farsingle(s,nodedis,nn,mm);
        farch1,farch2=farch(nodedis,index_this,mm);
        #println(G.m);#89429
        reccSfar,G=ReccS(s,farnode,G,t,s);
        #println(G.m);#89429
        #println(length(G.v));#89430
        reccSch,G=ReccS(farch1,farch2,G,t,s);
        xx=0;yy=0;
        reccSnew=0;
        if reccSfar < reccSch
            xx=s;yy=farnode;
            reccSnew=reccSfar;
            #println("最远好");
        else
            xx=farch1;yy=farch2;
            reccSnew=reccSch;
        end
        if numedg%5==0
            A[xx,yy]=1;A[yy,xx]=1;
            disfizuiyou=accRecc(A,n,s);
            #println(numedg,' ',disfizuiyou)
            println(disfizuiyou)
        end
        Gv=G.v;
        push!(Gv,xx);
        Gu=G.u;
        push!(Gu,yy);
        G=Graph(G.n, G.m+1, Gv, Gu,G.nbr);
        n=G.n;m=G.m;
        A=adjsp(G);
        B=getB(G); #m*n
        Q=randn(t,m); #t*m
        QB=Q*B; #t*n
        Z=zeros(t,n);
        f=approxchol_lap(A);
        for i=1:t
            Z[i,:]=f(QB[i,:]);
        end
        new_lab=Z'; #####n*t n个点的新坐标

        K2=500;
        kkkk=kkkk-0.08;
        index_this=bb2(new_lab,K2,kkkk);
        ll=length(index_this);
        println("凸包完成",ll);
        println("降维t  ",t)
    end
    


    # filename = "output.txt"
    # open(filename, "w") do file
    #     for i in 1:10
    #         println(file, disall100[i])
    #     end
    # end
end

    
   



 
    



