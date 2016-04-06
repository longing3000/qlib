function dynamics = StateEvolve( obj,hami1,hami2,liou1,liou2,state )
%STATEEVOLVE Summary of this function goes here
%   every evolution interval dt is split into N section: dt=N*tau
%   tau=tau1+tau2; tau2=p*tau;tau1=(1-p)*tau;
%   The time evolution for a single interval dt is determined by:
%   [exp(-1i*L2*tau2)exp(-1i*L1*tau1)]^N
%   This kind of evolution generates the quantum Zeno effect.

    import model.phy.Dynamics.QuantumDynamics
    import model.phy.Dynamics.EvolutionKernel.MatrixVectorEvolution
    import model.phy.Dynamics.EvolutionKernel.DensityMatrixEvolution
    para=obj.parameters;
    
    if strcmp(para.InitialStateType, 'MixedState')
        ops={liou1,liou2};
        nsection=para.nsection;
        p=para.proportion;        
        prefactors=ones(1,2*nsection);
        prefactors(1,1:2:2*nsection-1)=-1i*(1-p)/nsection;
        prefactors(1,2:2:2*nsection)=-1i*p/nsection;
        kern=MatrixVectorEvolution(ops, para.InitialStateType,prefactors);
    else
        ops={hami1,hami2};
        prefactors=[-1i,1i];
        kern=DensityMatrixEvolution(ops, para.InitialStateType,prefactors);

    end
    
    dynamics=QuantumDynamics( kern );
    dynamics.set_initial_state(state,'Liouville');
    dynamics.set_time_sequence(para.TimeList);
    
    dynamics.evolve();

end
