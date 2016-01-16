function [hami1, hami2, liou1,liou2]=GetHamiltonianLiouvillian(obj, spin_collection,matrix_strategy)
    import model.phy.QuantumOperator.SpinOperator.Hamiltonian
    import model.phy.SpinInteraction.SpinChainInteraction.OnSiteEnergy
    import model.phy.SpinInteraction.SpinChainInteraction.DQTInteraction
    import model.phy.SpinInteraction.SpinChainInteraction.XYInteraction
    import model.phy.SpinInteraction.SpinChainInteraction.dipInteraction
    
    hami1=Hamiltonian(spin_collection, matrix_strategy);
    if obj.parameters.AddOnSite1>0
        on_site_para.interaction=obj.parameters.onSite1;
        hami1.addInteraction( OnSiteEnergy(spin_collection, on_site_para));
    end
    if obj.parameters.AddDqtInt1>0
        dqt_para.interaction=obj.parameters.dqtInt1;
        hami1.addInteraction( DQTInteraction(spin_collection, dqt_para));
    end
    
    if obj.parameters.AddXYInt1>0
        xy_para.interaction=obj.parameters.xyInt1;
        hami1.addInteraction( XYInteraction(spin_collection, xy_para));
    end
    
    if obj.parameters.AddDipInt1>0
        dip_para.interaction=obj.parameters.dipInt1;
        hami1.addInteraction( dipInteraction(spin_collection, dip_para));
    end
    hami1.generate_matrix();
    
    hami2=Hamiltonian(spin_collection, matrix_strategy);
    if obj.parameters.AddOnSite2>0
        on_site_para.interaction=obj.parameters.onSite2;
        hami2.addInteraction( OnSiteEnergy(spin_collection, on_site_para));
    end
    if obj.parameters.AddDqtInt2>0
        dqt_para.interaction=obj.parameters.dqtInt2;
        hami2.addInteraction( DQTInteraction(spin_collection, dqt_para));
    end
    
    if obj.parameters.AddXYInt2>0
        xy_para.interaction=obj.parameters.xyInt2;
        hami2.addInteraction( XYInteraction(spin_collection, xy_para));
    end
    
    if obj.parameters.AddDipInt2>0
        dip_para.interaction=obj.parameters.dipInt2;
        hami2.addInteraction( dipInteraction(spin_collection, dip_para));
    end
    hami2.generate_matrix();

    if strcmp(obj.parameters.InitialStateType, 'MixedState')
        liou1=hami1.circleC();
        liou2=hami2.circleC();
    else
        liou1=[];
        liou2=[];
    end
end