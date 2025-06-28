import msprime
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import random
import scipy.stats as scst
import numpy as np 
import pandas as pd

plotdir = "/home/patel/embo_popgen_2025/!_MyWork/Day3/practicals/"

def make_pop_hist(pop1_size=125_000, pop2_size=50_000, anc_size=7_000_000, split_time=5000, seed=42, sample_size=50, recombination=8.4E-9, mutation=3.5E-9):
        demography = msprime.Demography()
        demography.add_population(name="Pop1", initial_size=pop1_size)
        demography.add_population(name="Pop2", initial_size=pop2_size)
        demography.add_population(name="Anc", initial_size=anc_size)
        demography.add_population_split(time=split_time, derived=["Pop1", "Pop2"], ancestral="Anc")
        # print(demography)
        
        ts = msprime.sim_ancestry(
        {"Pop1": sample_size, "Pop2": sample_size}, 
        demography=demography, 
        recombination_rate=recombination, # as in humans
        sequence_length=1_000,
        random_seed=seed)
        # print(ts)
        
        mts = msprime.sim_mutations(ts, rate=mutation, random_seed=seed)
        
        return mts

anc_size = 7_000_000
N1 = 150_000
N2 = 5000 
T_split = scst.randint(1000, 8000)

# sample estiamtes 
sample_fst = 0.2134

runs = 1000
tsplit_draws = np.zeros(runs)
fst_draws = np.zeros(runs)
seg1_draws = np.zeros(runs)
seg2_draws = np.zeros(runs)
tajima1_draws = np.zeros(runs)
tajima2_draws = np.zeros(runs)
pi1_draws = np.zeros(runs)
pi2_draws = np.zeros(runs)
        
for i in tqdm(range(runs)):
        t_split_sample = T_split.rvs()
        mts = make_pop_hist(
                pop1_size=N1,
                pop2_size=N2,
                split_time=t_split_sample
        )

        pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
        sample_sets=[mts.samples(pop_id["Pop1"]), mts.samples(pop_id["Pop2"])]

        Fst = mts.Fst(sample_sets)
        dxy = mts.divergence(sample_sets)
        segs = mts.segregating_sites(sample_sets)
        seg1 = segs[0] 
        seg2 = segs[1]
        tajimas = mts.Tajimas_D(sample_sets)
        tajima1 = tajimas[0] 
        tajima2 = tajimas[1]
        pis = mts.diversity(sample_sets)
        pi1 = pis[0] 
        pi2 = pis[1]
        # check for similarity
        # if (np.sqrt((Fst - sample_fst)**2) < 0.01):
        tsplit_draws[i] = t_split_sample
        fst_draws[i] = Fst
        seg1_draws[i] = seg1 
        seg2_draws[i] = seg2 
        tajima1_draws[i] = seg1 
        tajima2_draws[i] = seg2 
        pi1_draws[i] = pi1 
        pi2_draws[i] = pi2
        
data = {
        "t_split" : tsplit_draws,
        "fst" : fst_draws,
        "seg1" : seg1_draws,
        "seg2" : seg2_draws,
        "tajima1" : tajima1_draws,
        "tajima2" : tajima2_draws,
        "pi1" : pi1_draws,
        "pi2" : pi2_draws
} 

dat = pd.DataFrame(data)
dat.to_csv("tplit_samples.csv")
