import stdpopsim
import numpy as np
import sys
import os
import tskit


def simulate_data(contig, model, samples):
    engine = stdpopsim.get_engine("msprime")

    print("simulating...", end="", flush=True)
    ts = engine.simulate(
        model, contig, samples, random_seed=1,
    )
    print("done", flush=True)

    print("trees: ", ts.num_trees)
    print("length: ", ts.sequence_length)
    print("mutations: ", ts.num_mutations)
    print("sample sequences: ", ts.num_samples)

    return ts


def clean_relate():
    """
    Clears old Relate files
    :return:
    """
    os.system("rm -r relate*")
    os.system("rm relate_sim*")
    os.system("rm relate.*")
    os.system("rm ancestral_genome.fa")
    os.system("rm samples.vcf")
    os.system("rm sample_ages.txt")
    os.system("rm out.txt")


def run_relate(
    ts,
    rec_map_file,
    mutation_rate,
    Ne,
    relate_loc,
    mask=None,
    memory=5,
    threads=1,
):
    """
    Runs Relate on the input tree sequence
    """
    clean_relate()

    print("Writing files...", end="", flush=True)
    with open("samples.vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file)

    if mask is not None:
        mask = " --mask " + mask
    else:
        mask = ""

    # Hacky: ancestral genome is A except where we have mutations so know the ancestral state from the ts.
    ancestral_genome = ["A"] * int(ts.sequence_length)
    with open("ancestral_genome.fa", "w") as file:
        file.write(">ancestral_sequence\n")
        for s in ts.sites():
            ancestral_genome[int(s.position) - 1] = s.ancestral_state
        file.write("".join(ancestral_genome))
        file.write("\n")

    sample_ages = np.zeros(ts.num_samples)
    np.savetxt("sample_ages.txt", sample_ages, delimiter="\n")
    print("done", flush=True)

    S = "relate"
    Si = S + "_input"

    os.system(
        relate_loc
        + "/bin/RelateFileFormats --mode ConvertFromVcf --haps "
        + S
        + ".haps --sample "
        + S
        + ".sample -i samples;"
    )
    os.system(
        relate_loc
        + "/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps "
        + S
        + ".haps --sample "
        + S
        + ".sample --ancestor ancestral_genome.fa -o "
        + Si
        + mask
        + ";"
    )
    if threads == 1:
        os.system(
            relate_loc
            + "/bin/Relate --mode All --sample_ages sample_ages.txt "
            + "--memory "
            + str(memory)
            + " -m "
            + str(mutation_rate)
            + " -N "
            + str(Ne)
            + " --haps "
            + Si
            + ".haps.gz --sample "
            + Si
            + ".sample.gz --map "
            + str(rec_map_file)
            + " --seed 1 -o "
            + S
            + ";"
        )
    else:
        os.system(
            relate_loc
            + "/scripts/RelateParallel/RelateParallel.sh --sample_ages sample_ages.txt "
            + "--memory "
            + str(memory)
            + " -m "
            + str(mutation_rate)
            + " -N "
            + str(Ne)
            + " --haps "
            + Si
            + ".haps.gz --sample "
            + Si
            + ".sample.gz --map "
            + str(rec_map_file)
            + " --seed 1 -o "
            + S
            + " --threads "
            + str(threads)
            + " ;"
        )


