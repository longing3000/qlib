function calculate_tilde_coh_matrix(obj,cohmat,cluster_iter)
    subcluster_list=cluster_iter.cluster_info.subcluster_list;
    cluster_number_list=[cluster_iter.cluster_info.cluster_number_list{:,2}];

    ncluster=length(subcluster_list);
    ntime=length(cohmat(1,:));
    coh_tilde_mat=zeros(ncluster,ntime);
    coh_total=ones(1,ntime);
    coh=struct();

    cceorder=1;
    endpoints=cumsum(cluster_number_list);
    for m=1:ncluster
        subcluster=subcluster_list{m};
        nsubcluster=length(subcluster);
        coh_tilde=cohmat(m,:);
        if nsubcluster==0
            coh_tilde_mat(m,:)= coh_tilde;
        elseif nsubcluster>0
            for n=1:nsubcluster;
                coh_tilde_sub=coh_tilde_mat(subcluster(n),:);
                coh_tilde=coh_tilde./coh_tilde_sub;
            end
            coh_tilde_mat(m,:)=coh_tilde;
        end

        coh_total=coh_total.*coh_tilde;
        if m==endpoints(1,cceorder)
            field_name=strcat('coherence_cce_',num2str(cceorder));
            coh.(field_name)=coh_total;
            cceorder=cceorder+1;
        end
    end
    coh.('coherence')= coh_total;            
   if ncluster<20000          
       obj.keyVariables('coherence_tilde_matrix')=coh_tilde_mat;
   else
       timeTag=datestr(clock,'yyyymmdd_HHMMSS');
       save([OUTPUT_FILE_PATH, 'coherence_tilde_matrix', timeTag, '.mat'],'coh_tilde_mat');
       clear coh_tilde_mat;
   end
   coh.('timelist')=obj.parameters.TimeList;
    obj.keyVariables('coherence')=coh;
end