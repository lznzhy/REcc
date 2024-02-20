include("graph.jl")
include("core.jl")
include("bb.jl")
include("meth.jl")


using LinearAlgebra
using Laplacians
using SparseArrays
using Graphs

function turntype(G)
    g = Graphs.SimpleGraph(G.n)
    for i = 1:G.m
        Graphs.add_edge!(g,G.u[i],G.v[i])
    end
    return g;
end

function pagerank50(G)
    GG = turntype(G);
    pagerank_score = Graphs.pagerank(GG);
    pk=sortperm(pagerank_score);
    pkxx=pk[1];pkyy=pk[2];
    return pkxx,pkyy;
end

function betweenness50(G)
    GG = turntype(G);
    between_score = Graphs.betweenness_centrality(GG);
    be=sortperm(between_score);
    bexx=be[1];beyy=be[2];
    return bexx,beyy;
end

function closeness50(G)
    GG = turntype(G);
    close_score = Graphs.closeness_centrality(GG);
    cl=sortperm(close_score);
    clxx=cl[1];clyy=cl[2];
    return clxx,clyy;
end

function degree50(G)
    GG = turntype(G);
    degree_score = Graphs.degree_centrality(GG); 
    de=sortperm(degree_score);  
    dexx=de[1];deyy=de[2];
    return dexx,deyy;
end

function path50(G)
    GG = turntype(G);
    path_score = Graphs.floyd_warshall_shortest_paths(GG); 
    #println(path_score);
    #println(size(path_score.dists));
    #pa=sortperm(path_score);  
    aaa=argmax(path_score.dists);
    paxx=aaa[1];payy=aaa[2];
    #println(paxx);
    return paxx,payy;
end



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
    
    
    

    #pkall=[];beall=[];clall=[];deall=[];
    s=5791;
    # A=adjsp(G);
    # disfi=accRecc(A,n,s);
    # println(disfi);
    for k in 1:50
        xx,yy = path50(G);
        #println(xx,yy);
        #push!(pkall,(xx,yy));
        Gv=G.v;
        push!(Gv,xx);
        Gu=G.u;
        push!(Gu,yy);
        G=Graph(G.n, G.m+1, Gv, Gu, G.nbr);
        #Gc=findconnect(G)
        #G=Gc;
        n=G.n;m=G.m;
        if k%5==0
            A=adjsp(G);
            disfi2=accRecc(A,n,s);
            println(disfi2);
        end
    end

    # G1=G;
    # for k in 1:50
    #     xx,yy = betweenness50(G1);
    #     push!(beall,(xx,yy));
    #     G1v=G1.v;
    #     push!(G1v,xx);
    #     G1u=G1.u;
    #     push!(G1u,yy);
    #     G1=Graph(G1.n, G1.m+1, G1v, G1u,G1.nbr);
    #     n=G1.n;m=G1.m;
    # end

    # G1=G;
    # for k in 1:50
    #     xx,yy = closeness50(G1);
    #     push!(clall,(xx,yy));
    #     G1v=G1.v;
    #     push!(G1v,xx);
    #     G1u=G1.u;
    #     push!(G1u,yy);
    #     G1=Graph(G1.n, G1.m+1, G1v, G1u,G1.nbr);
    #     n=G1.n;m=G1.m;
    # end

    # G1=G;
    # for k in 1:50
    #     xx,yy = degree50(G1);
    #     push!(deall,(xx,yy));
    #     G1v=G1.v;
    #     push!(G1v,xx);
    #     G1u=G1.u;
    #     push!(G1u,yy);
    #     G1=Graph(G1.n, G1.m+1, G1v, G1u,G1.nbr);
    #     n=G1.n;m=G1.m;
    # end


    #println(pkall);
    # println(beall);
    # println(clall);
    # println(deall);

end
close(fname)
