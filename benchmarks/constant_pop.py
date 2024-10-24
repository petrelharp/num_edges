import stdpopsim
import tskit

engine = stdpopsim.get_engine("msprime")

species = stdpopsim.get_species("CanFam")
contig = species.get_contig('1', genetic_map="Campbell2016_CanFam3_1")

model = stdpopsim.PiecewiseConstantSize(1e4)

samples = {"pop_0": 10000}

ts = engine.simulate(model, contig, samples)

ts.dump("constant_pop.trees")

