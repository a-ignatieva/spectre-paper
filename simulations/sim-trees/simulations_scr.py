#!/usr/bin/env python
# coding: utf-8

import stdpopsim
import tskit
import tszip
import sys
import os

import simulations

chromosome = sys.argv[1]  # e.g. chr21

mask = "../mask/20140520." + chromosome + ".pilot_mask.fasta.gz"
rec_map_file = "../HapMapII_GRCh37/genetic_map_GRCh37_" + chromosome + ".txt"
relate_loc = "~/relate"
bcftools_loc = "~/bcftools-1.21/bin/bcftools"
memory = 40
threads = 8

species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("AmericanAdmixture_4B18")
contig = species.get_contig(chromosome=chromosome, genetic_map="HapMapII_GRCh37", mutation_rate=model.mutation_rate)
mutation_rate = contig.mutation_rate
recombination_map = contig.recombination_map
samples = {"EUR": int((2*20000 + 1006)/2)}
Ne = species.population_size  # diploids

print("Simulating " + chromosome, flush=True)
print("mutation rate " + str(mutation_rate), flush=True)
print("mean recombination rate " + str(recombination_map.mean_rate), flush=True)
print("total samples " + str(samples), flush=True)
ts = simulations.simulate_data(
    contig,
    model,
    samples,
)

tszip.compress(
    ts,
    "sim_" + chromosome + ".trees.tsz",
)

# Split into three sets:
# ts_relate will be the ARG of 1006 samples used for spectre testing
# ts_data will be the data used in interaction testing
# ts_replicate will be used for replicating the interaction testing results
print("Subsetting simulated ts")
ts_relate = ts.simplify(samples=[i for i in range(1006)], record_provenance=False)
ts_data = ts.simplify(samples=[i for i in range(1006, 1006 + int((ts.num_samples - 1006)/2))], record_provenance=False)
ts_replicate = ts.simplify(samples=[i for i in range(1006 + int((ts.num_samples - 1006)/2), ts.num_samples)], record_provenance=False)
print("number of samples for Relate:", ts_relate.num_samples)
print("number of samples for interactions test:", ts_data.num_samples)
print("number of samples for replication test:", ts_replicate.num_samples)

# Store everything
tszip.compress(
    ts_relate,
    "sim_ARG_" + chromosome + ".trees.tsz",
)
tszip.compress(
    ts_data,
    "sim_data_" + chromosome + ".trees.tsz",
)
tszip.compress(
    ts_replicate,
    "sim_replicate_" + chromosome + ".trees.tsz",
)
with open("sim_data_" + chromosome + ".vcf", "w") as vcf_file:
    ts_data.write_vcf(vcf_file)
os.system(bcftools_loc + " view -O b sim_data_" + chromosome + ".vcf > sim_data_" + chromosome + ".bcf")
os.system("rm sim_data_" + chromosome + ".vcf")
with open("sim_replicate_" + chromosome + ".vcf", "w") as vcf_file:
    ts_replicate.write_vcf(vcf_file)
os.system(bcftools_loc + " view -O b sim_replicate_" + chromosome + ".vcf > sim_replicate_" + chromosome + ".bcf")
os.system("rm sim_replicate_" + chromosome + ".vcf")

# Run Relate on the ARG of 1006 samples
print("Running Relate", flush=True)
simulations.run_relate(
    ts_relate,
    rec_map_file,
    mutation_rate,
    2 * Ne,
    relate_loc,
    mask=mask,
    memory=int(memory / threads),
    threads=threads,
)
os.system(
    "mv relate.anc "
    + "sim_relate_"
    + chromosome
    + ".anc"
)
os.system(
    "mv relate.mut "
    + "sim_relate_"
    + chromosome
    + ".mut"
)
os.system(
    "gzip "
    + "sim_relate_"
    + chromosome
    + ".anc"
)
os.system(
    "gzip "
    + "sim_relate_"
    + chromosome
    + ".mut"
)
os.system(
    relate_loc
    + "/bin/RelateFileFormats --mode ConvertToTreeSequence -i sim_relate_"
    + chromosome
    + " -o sim_relate_"
    + chromosome
    + ";"
)
print("Running final cleanup", flush=True)
simulations.clean_relate()
print("Done", flush=True)
