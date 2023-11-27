function classes = optimized_kmeans_clustering(X, nmax)

klist=2:nmax;%the number of clusters you want to try
myfunc = @(X,K)(kmeans(X, K));
eva = evalclusters(X,myfunc,'CalinskiHarabasz','klist',klist);
classes=kmeans(X,eva.OptimalK);