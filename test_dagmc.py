from IPython.display import Image
import openmc
import math
from matplotlib import pyplot as plt

from warnings import warn
import openmc.lib

if not openmc.lib.feature_enabled("dagmc"):
    raise RuntimeError("DAGMC is not enabled.")

print("DAGMC is enabled!")

import urllib.request

fuel_pin_url = 'https://tinyurl.com/y3ugwz6w' 
teapot_url = 'https://tinyurl.com/y4mcmc3u'
paramak_w_graveyard = 'https://tinyurl.com/y7w238th' 

def download(url, filename="dagmc.h5m"):
    """
    Helper function for retrieving dagmc models
    """
    u = urllib.request.urlopen(url)
    
    if u.status != 200:
        raise RuntimeError("Failed to download file.")
    
    # save file as dagmc.h5m
    with open(filename, 'wb') as f:
        f.write(u.read())

 # materials
u235 = openmc.Material(name="fuel")
u235.add_nuclide('U235', 1.0, 'ao')
u235.set_density('g/cc', 11)
u235.id = 40

water = openmc.Material(name="water")
water.add_nuclide('H1', 2.0, 'ao')
water.add_nuclide('O16', 1.0, 'ao')
water.set_density('g/cc', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')
water.id = 41

materials = openmc.Materials([u235, water])

download(fuel_pin_url)
dagmc_univ = openmc.DAGMCUniverse(filename="dagmc.h5m")
geometry = openmc.Geometry(root=dagmc_univ)

settings = openmc.Settings()
settings.batches = 10
settings.inactive = 2
settings.particles = 5000

model = openmc.Model(geometry=geometry, materials=materials, settings=settings)
model.export_to_model_xml()

p = openmc.Plot()
p.width = (25.0, 25.0)
p.pixels = (400, 400)
p.color_by = 'material'
p.colors = {u235: 'yellow', water: 'blue'}
openmc.plot_inline(
    p,
    openmc_exec="/Users/chrispattlesmith/dagmc-deps/openmc/bin/openmc"
)

model.settings.source = openmc.IndependentSource(space=openmc.stats.Box([-4., -4., -4.],
                                                       [ 4.,  4.,  4.]))
model.export_to_model_xml()

tally = openmc.Tally()
tally.scores = ['total']
tally.filters = [openmc.CellFilter(1)]
tallies = openmc.Tallies([tally])

model.tallies = tallies
model.export_to_model_xml()

download(teapot_url, "teapot.h5m")
dagmc_univ = openmc.DAGMCUniverse(filename="teapot.h5m")
geometry = openmc.Geometry(root=dagmc_univ)

print("Number of DAGMC cells:", dagmc_univ.n_cells)
print("Number of DAGMC surfaces:", dagmc_univ.n_surfaces)
print("Materials",dagmc_univ.material_names)

iron = openmc.Material(name="iron")
iron.add_nuclide("Fe54", 0.0564555822608)
iron.add_nuclide("Fe56", 0.919015287728)
iron.add_nuclide("Fe57", 0.0216036861685)
iron.add_nuclide("Fe58", 0.00292544384231)
iron.set_density("g/cm3", 7.874)

water = openmc.Material(name="water")
water.add_nuclide('H1', 2.0, 'ao')
water.add_nuclide('O16', 1.0, 'ao')
water.set_density('g/cc', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')
water.id = 41

materials = openmc.Materials([iron, water])

settings = openmc.Settings()
settings.batches = 10
settings.particles = 5000
settings.run_mode = "fixed source"
settings.source_rejection_fraction = 0.01

space = openmc.stats.Box(*geometry.bounding_box)
my_source = openmc.IndependentSource(space=space, domains=[water])
my_source.energy = openmc.stats.Discrete(x=[12.0,], p=[1.0,])

settings.source = my_source

mesh = openmc.RegularMesh()  # could use from_domain here to avid setting lower_left and upper_right
mesh.dimension = (120, 1, 40)
mesh.lower_left = (-20.0, 0.0, -10.0)
mesh.upper_right = (20.0, 1.0, 4.0)

mesh_filter = openmc.MeshFilter(mesh)

pot_filter = openmc.CellFilter([1])
pot_tally = openmc.Tally()
pot_tally.filters = [mesh_filter, pot_filter]
pot_tally.scores = ['flux']

water_filter = openmc.CellFilter([5])
water_tally = openmc.Tally()
water_tally.filters = [mesh_filter, water_filter]
water_tally.scores = ['flux']


tallies = openmc.Tallies([pot_tally, water_tally])

	
model = openmc.Model(materials=materials, geometry=geometry, settings=settings, tallies=tallies)
model.export_to_model_xml()

p = openmc.Plot()
p.basis = 'xz'
p.origin = (0.0, -0.1, 0.0)
p.width = (30.0, 20.0)
p.pixels = (450, 300)
p.color_by = 'material'
p.colors = {iron: 'gray', water: 'blue'}
openmc.plot_inline(
    p,
    openmc_exec="/Users/chrispattlesmith/dagmc-deps/openmc/bin/openmc"
)

openmc.run(
    openmc_exec="/Users/chrispattlesmith/dagmc-deps/openmc/bin/openmc"
)

sp = openmc.StatePoint("statepoint.10.h5")

water_tally = sp.get_tally(scores=['flux'], id=water_tally.id)
water_flux = water_tally.mean
water_flux.shape = (40, 120)
water_flux = water_flux[::-1, :]

pot_tally = sp.get_tally(scores=['flux'], id=pot_tally.id)
pot_flux = pot_tally.mean
pot_flux.shape = (40, 120)
pot_flux = pot_flux[::-1, :]

sp.close()

from matplotlib import pyplot as plt

fig = plt.figure(figsize=(18, 16))

sub_plot1 = plt.subplot(121, title="Kettle Flux")
sub_plot1.imshow(pot_flux)

sub_plot2 = plt.subplot(122, title="Water Flux")
sub_plot2.imshow(water_flux)

plt.tight_layout()
plt.savefig("flux_results.png", dpi=300)
plt.show()

