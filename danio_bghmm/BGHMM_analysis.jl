#JOB FILEPATHS
sample_output = "/bench/PhD/NGS_binaries/BGHMM/BGHMM_samples"
hmm_output = "/bench/PhD/NGS_binaries/BGHMM/hmmchains"
selected_hmm_output = "/bench/PhD/NGS_binaries/BGHMM/selected_BGHMMs"

#GENERAL SETUP
@info "Loading libraries..."
using BGHMM, BioSequences, DataFrames, Distributions, CLHMM, ProgressMeter, Plots, Serialization, Statistics

const replicates = 4 #repeat optimisation from this many seperately initialised samples from the prior
const Ks = [1,2,4,6] #mosaic class #s to test
const order_nos = [0,1,2] #DNA kmer order #s to test
const partitions = ["exon", "intergenic", "periexonic"]
const run_samples = 10000

#LOAD SAMPLES
@info "Loading samples from $sample_output..."
sample_dfs = deserialize(sample_output)

#LOAD HMM CHAINS
@info "Loading hmm chains from $hmm_output"
hmm_results_dict = deserialize(hmm_output)

#BUILD TRAINING AND TEST SETS FROM SAMPLES
training_sets, test_sets = BGHMM.split_obs_sets(sample_dfs)

naive_likelihood_dict = Dict()
naive_hmm = HMM(ones(1,1),[Categorical(4)])

#For all partitions, calculate the test set's probability given the naive model
@showprogress 1 "Calculating naive model test set likelihoods..." for partition in partitions
    naive_likelihood_dict[partition] = BGHMM.test_hmm(naive_hmm, test_sets[partition], 0)
end

#COMPOSE DICT OF TEST SET LIKELIHOODS GIVEN LAST HMM IN EACH CHAIN, BY JOBID
hmm_likelihoods_dict = Dict()
@showprogress 1 "Testing HMMs..." for (jobid, hmm_chain) in hmm_results_dict
    last_hmm = hmm_chain[end][2]
    partition_seqs = test_sets[jobid[1]]
    order = jobid[3]
    hmm_likelihoods_dict[jobid] = BGHMM.test_hmm(last_hmm, partition_seqs, order)
end

#COMPOSE DICT OF MAX STATE RUN LENGTH GIVEN LAST HMM DIAGONALS IN EACH CHAIN, BY JOBID
hmm_max_run_length_dict = Dict()
@showprogress 1 "Simulating run lengths..." for (jobid, hmm_chain) in hmm_results_dict
    last_hmm = hmm_chain[end][2]
    diagonal = BGHMM.get_diagonal_array(last_hmm)
    diagonal_mrls = BGHMM.sim_run_lengths(diagonal, run_samples)
    hmm_max_run_length_dict[jobid] = maximum(diagonal_mrls)
end

#INITIALIZE DATA MATRICES AND COMPOSE VALUES FOR PLOTTING
data_matrix_dict = Dict()
for order in order_nos
    data_matrix_dict[order] = Array{Union{Float64, Int64, Symbol}}(undef, length(Ks)*replicates, length(partitions), 4)
end
#Iterating over hmm_likelihoods_dict, compose the data matrices for plots
for (jobid, likelihood) in hmm_likelihoods_dict
    partition, K, order, replicate = jobid
    converged = hmm_results_dict[jobid][end][5]

    entry_coord = replicate + ((findfirst(isequal(K),Ks) - 1) * replicates)
    partition_series = findfirst(isequal(partition),partitions)

    data_matrix_dict[order][entry_coord,partition_series,1] = likelihood - naive_likelihood_dict[partition]
    data_matrix_dict[order][entry_coord,partition_series,2] = K
    data_matrix_dict[order][entry_coord,partition_series,3] = hmm_max_run_length_dict[jobid]

    if converged
        data_matrix_dict[order][entry_coord,partition_series,4] = :black
    else
        data_matrix_dict[order][entry_coord,partition_series,4] = :white
    end
end

#SETUP PLOTS
zeroOlh = plot(data_matrix_dict[0][:,:,2], title="0th order", label=["exon", "intergenic", "periexonic"], data_matrix_dict[0][:,:,1],markercolors=data_matrix_dict[0][:,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], legend=:bottomright, xlims=[0.9,6.1], xticks=[1,2,4,6], xlabel="classes", ylabel="likelihood vs naive")

zeroOrun = plot(data_matrix_dict[0][(replicates+1):end,:,2], title="0th order", label=["exon", "intergenic", "periexonic"], data_matrix_dict[0][(replicates+1):end,:,3],markercolors=data_matrix_dict[0][(replicates+1):end,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], xlims=[1.9,6.1], xticks=[2,4,6], xlabel = "classes", ylabel="model maximum mean(SRL)")

firstOlh = plot(data_matrix_dict[1][:,:,2],  title="1st order", label=["exon", "intergenic", "periexonic"],data_matrix_dict[1][:,:,1],markercolors=data_matrix_dict[1][:,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], legend=:bottomright,  xlims=[0.9,6.1], xticks=[1,2,4,6], xlabel="classes", ylabel="likelihood vs naive")

firstOrun = plot(data_matrix_dict[1][(replicates+1):end,:,2], title="1st order", label=["exon", "intergenic", "periexonic"], data_matrix_dict[1][(replicates+1):end,:,3],markercolors=data_matrix_dict[1][(replicates+1):end,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle],  xlims=[1.9,6.1], xticks=[2,4,6], xlabel = "classes", ylabel="model maximum mean(SRL)")

secondOlh = plot(data_matrix_dict[2][:,:,2], title="2nd order", label=["exon", "intergenic", "periexonic"], data_matrix_dict[2][:,:,1],markercolors=data_matrix_dict[2][:,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle],  xlims=[0.9,6.1], xticks=[1,2,4,6], legend=:bottomright, xlabel="classes", ylabel="likelihood vs naive")

secondOrun = plot(data_matrix_dict[2][(replicates+1):end,:,2], title="2nd order", label=["exon", "intergenic", "periexonic"], data_matrix_dict[2][(replicates+1):end,:,3],markercolors=data_matrix_dict[2][(replicates+1):end,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], xlims=[1.9,6.1], xticks=[2,4,6], xlabel = "classes", ylabel="model maximum mean(SRL)")

#SAVE PLOTS
plots_to_print = Dict("zeroOlh"=>zeroOlh, "zeroOrun"=>zeroOrun, "firstOlh"=>firstOlh, "firstOrun"=>firstOrun, "secondOlh"=>secondOlh, "secondOrun"=>secondOrun)

for (filename, plot) in plots_to_print
    png(plot, filename)
end
exp
#SELECT BEST MODELS FOR PARTITIONS
BGHMM_dict = Dict{String,Tuple{HMM, Int64, Float64}}()
for (jobid, likelihood) in hmm_likelihoods_dict
    partition = jobid[1]
    order = jobid[3]
    if !haskey(BGHMM_dict, partition)
        BGHMM_dict[partition] = (hmm_results_dict[jobid][end][2], order, likelihood)
    else
        if likelihood > BGHMM_dict[partition][3]
            BGHMM_dict[partition] = (hmm_results_dict[jobid][end][2], order, likelihood)
        end
    end
end

@info "Best models by partition:"
for partition in partitions
    hmm = BGHMM_dict[partition][1]
    println("$partition :")
    println(hmm)
end

#SERIALIZE MODELS FOR USE IN NESTED SAMPLING
serialize(selected_hmm_output, BGHMM_dict)
