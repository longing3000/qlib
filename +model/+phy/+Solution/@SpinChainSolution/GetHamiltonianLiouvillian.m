function [hami, liou]=GetHamiltonianLiouvillian(obj, spin_collection,matrix_strategy)
    import model.phy.QuantumOperator.SpinOperator.Hamiltonian
    import model.phy.SpinInteraction.SpinChainInteraction.OnSiteEnergy
    import model.phy.SpinInteraction.SpinChainInteraction.DQTInteraction
    import model.phy.SpinInteraction.SpinChainInteraction.XYInteraction
    import model.phy.SpinInteraction.SpinChainInteraction.dipInteraction
    
    hami=Hamiltonian(spin_collection, matrix_strategy);
    if obj.parameters.AddOnSite>0
        on_site_para.interaction=obj.parameters.onSite;
        hami.addInteraction( OnSiteEnergy(spin_collection, on_site_para));
    end
    if obj.parameters.AddDqtInt>0
        dqt_para.interaction=obj.parameters.dqtInt;
        hami.addInteraction( DQTInteraction(spin_collection, dqt_para));
    end
    
    if obj.parameters.AddXYInt>0
        xy_para.interaction=obj.parameters.xyInt;
        hami.addInteraction( XYInteraction(spin_collection, xy_para));
    end
    
    if obj.parameters.AddDipInt>0
        dip_para.interaction=obj.parameters.dipInt;
        hami.addInteraction( dipInteraction(spin_collection, dip_para));
    end
    hami.generate_matrix();

    if strcmp(obj.parameters.InitialStateType, 'MixedState')
        liou=hami.circleC();
    else
        liou=[];
    end
end