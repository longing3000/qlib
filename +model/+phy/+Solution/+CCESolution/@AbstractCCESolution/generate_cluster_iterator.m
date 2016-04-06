function cluster_iterator=generate_cluster_iterator(obj)
    import model.phy.SpinCollection.Strategy.FromFile
    import model.phy.SpinCollection.Strategy.FromSpinList
    import model.phy.SpinCollection.SpinCollection            
    import model.phy.SpinCollection.Iterator.ClusterIterator
    import model.phy.SpinCollection.Iterator.ClusterIteratorGen.CCE_Clustering

    if obj.parameters.load_cluster_iter             
         cluster_iterator=obj.keyVariables('cluster_iterator');
         if obj.parameters.SetBathSpins.SetSpin;
             paraCell=obj.parameters.SetBathSpins.BathSpinsSettingCell;
             cluster_iterator.spin_collection.set_spin(paraCell);
         else
             cluster_iterator.spin_collection.set_spin();
         end
    else             
          switch obj.parameters.SpinCollectionStrategy
               case 'File'
                   spin_collection=SpinCollection( FromFile(...
                       [INPUT_FILE_PATH, obj.parameters.InputFile]));
                   if obj.parameters.SetBathSpins.SetSpin;
                       paraCell=obj.parameters.SetBathSpins.BathSpinsSettingCell;
                       spin_collection.set_spin(paraCell);
                   else
                       spin_collection.set_spin();
                   end
               case 'SpinList'
                   error('not surported so far.');
           end

           add_gradient_field=obj.parameters.add_gradient_field;
           if add_gradient_field
               field_gradient=obj.parameters.field_gradient;
               magnetic_field=obj.parameters.MagneticField;
               field_direction=magnetic_field/norm(magnetic_field);%get the direction of the external homogenous field
               spin_collection=set_gradient_field(spin_collection,field_gradient,field_direction);
           end
           obj.keyVariables('spin_collection')=spin_collection;

           clu_para.cutoff=obj.parameters.CutOff;
           clu_para.max_order=obj.parameters.MaxOrder;
           disp('clustering begins...')
           tic
           % the cce strategy can be change
           cce=CCE_Clustering(spin_collection, clu_para);

           cluster_iterator=ClusterIterator(spin_collection,cce);
           [~]=cluster_iterator.cross_relation_gen();
           obj.keyVariables('cluster_iterator')=cluster_iterator;
           save([OUTPUT_FILE_PATH, 'cluster_iterator', obj.timeTag, '.mat'],'cluster_iterator');
           toc
     end
end

function spin_collection=set_gradient_field(spin_collection,field_gradient,field_direction)  
    nspin=spin_collection.getLength;
    for kk=1:nspin
        spin=spin_collection.spin_list{kk};
        local_field=spin.coordinate*field_gradient';
        spin.local_field=local_field*field_direction;
    end
end