

% sol=model.phy.Solution.DipolarCoupledSpinEvolution('DipolarSpinDynamicsButane.xml');
sol=model.phy.Solution.EnsembleCCESolution('EnsembleCCE_hBNLayer.xml');
sol.perform();
sol.save_solution();