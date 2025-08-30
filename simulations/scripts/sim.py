import crpropa

RESULTS_DIR = "../results"

# ----------------------------------------------------------------------------------------------------
def simulation():

	scatter_velocity = 0.1 * crpropa.c_light
	step_length = 0.5 * crpropa.parsec

	flow_direction = crpropa.Vector3d(1., 0, 0) * scatter_velocity
	yzsize = 100. * crpropa.parsec

	upstream_velocity = flow_direction
	upstreamSize = 10000. * crpropa.parsec
	upstreamGeometry = crpropa.ParaxialBox(crpropa.Vector3d(-1 * upstreamSize, -.5 * yzsize, -.5 * yzsize),
		crpropa.Vector3d(upstreamSize, yzsize, yzsize))
	upstream_scatter_module = crpropa.DirectedFlowScattering(upstream_velocity, step_length)
	upstream = crpropa.RestrictToRegion(upstream_scatter_module, upstreamGeometry)

	downstreamSize = 100. * crpropa.parsec
	downstream_velocity = flow_direction * 1./4
	downstreamGeometry = crpropa.ParaxialBox(crpropa.Vector3d(0, -.5 * yzsize, -.5 * yzsize),
		crpropa.Vector3d(downstreamSize, yzsize, yzsize))
	downstream_scatter_module = crpropa.DirectedFlowScattering(downstream_velocity, step_length)
	downstream = crpropa.RestrictToRegion(downstream_scatter_module, downstreamGeometry)

	simulation = crpropa.ModuleList()
	simulation.add(upstream)
	simulation.add(downstream)
	simulation.add(crpropa.ReflectiveBox(crpropa.Vector3d(-upstreamSize * 2, -yzsize /2, -yzsize /2),
		crpropa.Vector3d(upstreamSize * 2 + downstreamSize * 2, yzsize, yzsize)))

	simulation.add(crpropa.SimplePropagation(1E-4 *crpropa.parsec, .5 *crpropa.parsec))
	obs1 = crpropa.Observer()
	obs1.add(crpropa.ObserverSurface(crpropa.Plane(crpropa.Vector3d(-upstreamSize, 0, 0), crpropa.Vector3d(1., 0, 0))))
	obs1.setDeactivateOnDetection(True)
	output1 = crpropa.TextOutput(f'{RESULTS_DIR}/shock_upstream.txt', crpropa.Output.Event3D)
	obs1.onDetection(output1)
	simulation.add(obs1)

	obs2 = crpropa.Observer()
	obs2.add(crpropa.ObserverSurface(crpropa.Plane(crpropa.Vector3d(downstreamSize, 0, 0), crpropa.Vector3d(1., 0, 0))))
	obs2.setDeactivateOnDetection(True)
	output2 = crpropa.TextOutput(f'{RESULTS_DIR}/shock_downstream.txt', crpropa.Output.Event3D)
	obs2.onDetection(output2)
	simulation.add(obs2)

	source = crpropa.Source()
	source.add(crpropa.SourcePosition(crpropa.Vector3d(-10. * crpropa.parsec, 0, 0)))
	source.add(crpropa.SourceParticleType(crpropa.nucleusId(1, 1)))
	source.add(crpropa.SourceEnergy(1E16 * crpropa.eV))
	source.add(crpropa.SourceIsotropicEmission())

	# Execute simulation
	simulation.setShowProgress(True)
	simulation.run(source, 10000)
	output1.close()
	output2.close()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

	simulation()

# ----------------------------------------------------------------------------------------------------
