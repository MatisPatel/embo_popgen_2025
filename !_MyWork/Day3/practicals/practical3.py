import msprime
import demesdraw
from matplotlib import pyplot as plt

plotdir = "/home/patel/embo_popgen_2025/!_MyWork/Day3/practicals/"

demography = msprime.Demography()
demography.add_population(name="Pop1", initial_size=125_000)
demography.add_population(name="Pop2", initial_size=40_000)
demography.add_population(name="Anc", initial_size=7_000_000)
demography.add_population_split(time=5000, derived=["Pop1", "Pop2"], ancestral="Anc")
print(demography)

ts = msprime.sim_ancestry({"Pop1": 50, "Pop2": 50}, demography=demography, random_seed=42)
print(ts)

plt.clf()
demesdraw.tubes(demography.to_demes(), ax=plt.gca(), seed=1, log_time=True)
plt.savefig(plotdir + "demes_plot_example.png")

ts = msprime.sim_ancestry(
        {"Pop1": 50, "Pop2": 50}, 
        demography=demography, 
        recombination_rate=8.4E-9, # as in humans
        sequence_length=1_000,
        random_seed=124234)
print(ts)

mts = msprime.sim_mutations(ts, rate=3.5E-9, random_seed=1234)
print(mts.tables.sites)

# visualise the site frequency spectrum
plt.clf()
afs = mts.allele_frequency_spectrum()
plt.bar(range(mts.num_samples + 1), afs)
plt.title("Allele frequency spectrum")
plt.savefig(plotdir+"All_freq_spec.png")

# Define the samples between which Fst will be calculated
pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
sample_sets=[mts.samples(pop_id["Pop1"]), mts.samples(pop_id["Pop2"])]

print(mts.Fst(sample_sets))