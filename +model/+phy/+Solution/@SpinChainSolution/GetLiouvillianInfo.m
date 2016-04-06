function  [liouvillian_sorted,state_list_sorted]=GetLiouvillianInfo( obj )
%GETLIOUVILLIANINFO Summary of this function goes here
%   Detailed explanation goes here
    import model.phy.SpinCollection.Iterator.ClusterIterator
    import model.phy.SpinCollection.Iterator.ClusterIteratorGen.RSS_Clustering
    
    clu_para.cutoff=obj.parameters.CutOff;
    clu_para.max_order=obj.parameters.MaxOrder;
    spin_collection=obj.keyVariables('spin_collection');    
    rss=RSS_Clustering(spin_collection,clu_para);
    [~]=rss.generate_clusters();
    state_list=rss.generate_state_list;
    
    liou_mat=sparse(obj.result.liouvillian.matrix);
    IST_ts_mat=obj.result.IST_transform_operator.getMatrix;
    IST_ts_mat=sparse(IST_ts_mat);
    liou_mat=IST_ts_mat*liou_mat*IST_ts_mat';
    [liouvillian_sorted,state_list_sorted]=resort_liouvillian( liou_mat,state_list );
    obj.StoreKeyVariables(liouvillian_sorted,state_list_sorted);
end

function [liouvillian_sorted,state_list_sorted]=resort_liouvillian( liou_mat,state_list )
%  sort the states according to the none zero numbers of the states
[nr,nc]=size(state_list);
nz_list=zeros(nr,1);%none zero list of state list
for n=1:nr
    nz_list(n,1)=nnz(state_list(n,:));
end

A=[state_list,nz_list];
col_index=fliplr(1:(nc+1));
[A,cols]=sortrows(A,col_index);

state_list_sorted=A(:,1:nc);

liouvillian_sorted=liou_mat(cols,cols);
end

