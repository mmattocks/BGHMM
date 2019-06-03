#JOB FILEPATHS
Sys.islinux() ? hmm_output = "/media/main/Bench/PhD/NGS_binaries/BGHMM/hmmchains" : hmm_output = "F:\\PhD\\NGS_binaries\\BGHMM\\hmmchains"
Sys.islinux() ? sample_output = "/media/main/Bench/PhD/NGS_binaries/BGHMM/BGHMM_samples" : sample_output = "F:\\PhD\\NGS_binaries\\BGHMM\\BGHMM_samples" #path to sequence samples for learning
Sys.islinux() ? selected_hmm_output = "/media/main/Bench/PhD/NGS_binaries/BGHMM/selected_BGHMMs" : selected_hmm_output = "F:\\PhD\\NGS_binaries\\BGHMM\\selected_BGHMMs"

#GENERAL SETUP
@info "Loading libraries..."
using BGHMM, BioSequences, DataFrames, Distributions, MS_HMMBase, ProgressMeter, Plots, Serialization, Statistics

const replicates = 3 #repeat optimisation from this many seperately initialised samples from the prior
const Ks = [1,2,4] #mosaic class #s to test
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
training_sets, test_sets = Main.BGHMM.split_obs_sets(sample_dfs)

naive_likelihood_dict = Dict()
naive_hmm = MS_HMMBase.HMM(ones(1,1),[Categorical(4)])

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

#COMPOSE DICT OF MEAN STATE RUN LENGTH GIVEN LAST HMM DIAGONALS IN EACH CHAIN, BY JOBID
hmm_mean_run_length_dict = Dict()
@showprogress 1 "Simulating run lengths..." for (jobid, hmm_chain) in hmm_results_dict
    last_hmm = hmm_chain[end][2]
    diagonal = BGHMM.get_diagonal_array(last_hmm)
    diagonal_mrls = BGHMM.sim_run_lengths(diagonal, run_samples)
    hmm_mean_run_length_dict[jobid] = mean(diagonal_mrls)
end

#INITIALIZE DATA MATRICES AND COMPOSE VALUES FOR PLOTTING
data_matrix_dict = Dict()
for partition in partitions
    data_matrix_dict[partition] = Array{Union{Float64, Int64, Symbol}}(undef, length(Ks)*replicates, length(order_nos), 4)
end
#Iterating over hmm_likelihoods_dict, compose the data matrices for plots
for (jobid, likelihood) in hmm_likelihoods_dict
    partition, K, order, replicate = jobid
    converged = hmm_results_dict[jobid][end][5]

    entry_coord = replicate + ((findfirst(isequal(K),Ks) - 1) * replicates)
    order_series = findfirst(isequal(order),order_nos)

    data_matrix_dict[partition][entry_coord,order_series,1] = likelihood - naive_likelihood_dict[partition]
    data_matrix_dict[partition][entry_coord,order_series,2] = K
    data_matrix_dict[partition][entry_coord,order_series,3] = hmm_mean_run_length_dict[jobid]

    if converged
        data_matrix_dict[partition][entry_coord,order_series,4] = :black
    else
        data_matrix_dict[partition][entry_coord,order_series,4] = :white
    end
end

pyplot()
#PLOT LAYOUT
l = @layout [a{.5w} b{.2w} c{.2w}
             d e f
             g h i]

#SETUP PLOTS
exl = plot(data_matrix_dict["exon"][:,:,2], label=["0th order", "1st order", "2nd order"], data_matrix_dict["exon"][:,:,1],markercolors=data_matrix_dict["exon"][:,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], legend=:bottomright, xlabel="classes", ylabel="likelihood vs naive")
exlinset = scatter(data_matrix_dict["exon"][(replicates+1):end,:,2], data_matrix_dict["exon"][(replicates+1):end,:,1],markercolors=data_matrix_dict["exon"][(replicates+1):end,:,4], markershape=[:circle :rect :utriangle], xlims=[1.9,4.1], xticks=[2,4], xlabel="classes")
exr = plot(data_matrix_dict["exon"][(replicates+1):end,:,2],data_matrix_dict["exon"][(replicates+1):end,:,3],markercolors=data_matrix_dict["exon"][(replicates+1):end,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], xlims=[1.9,4.1], xticks=[2,4], xlabel = "classes", ylabel="mean SRL")
inl = plot(data_matrix_dict["intergenic"][:,:,2], label=["0th order", "1st order", "2nd order"], data_matrix_dict["intergenic"][:,:,1],markercolors=data_matrix_dict["intergenic"][:,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], legend=false, xlabel="classes", ylabel="likelihood vs naive")
inlinset = scatter(data_matrix_dict["intergenic"][(replicates+1):end,:,2], data_matrix_dict["intergenic"][(replicates+1):end,:,1],markercolors=data_matrix_dict["intergenic"][(replicates+1):end,:,4], markershape=[:circle :rect :utriangle], xlims=[1.9,4.1], xticks=[2,4], xlabel="classes")
inr = plot(data_matrix_dict["intergenic"][(replicates+1):end,:,2],data_matrix_dict["intergenic"][(replicates+1):end,:,3],markercolors=data_matrix_dict["intergenic"][(replicates+1):end,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], xlabel = "classes", ylabel="mean SRL")
pel = plot(data_matrix_dict["periexonic"][:,:,2], label=["0th order", "1st order", "2nd order"], data_matrix_dict["periexonic"][:,:,1],markercolors=data_matrix_dict["periexonic"][:,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], legend=false, xlabel="classes", ylabel="likelihood vs naive")
pelinset = scatter(data_matrix_dict["periexonic"][(replicates+1):end,:,2], data_matrix_dict["periexonic"][(replicates+1):end,:,1],markercolors=data_matrix_dict["periexonic"][(replicates+1):end,:,4], markershape=[:circle :rect :utriangle], xlims=[1.9,4.1], xticks=[2,4], xlabel="classes")
per = plot(data_matrix_dict["periexonic"][(replicates+1):end,:,2],data_matrix_dict["periexonic"][(replicates+1):end,:,3],markercolors=data_matrix_dict["periexonic"][(replicates+1):end,:,4], line=(3, [:solid :dash :dot]), markershape=[:circle :rect :utriangle], xlims=[1.9,4.1], xticks=[2,4], xlabel = "classes", ylabel="mean SRL")

#ASSEMBLE & SAVE DIAGNOSTIC PLOT
plot(exl, exlinset, exr,
    inl, inlinset, inr,
    pel, pelinset, per,
    layout=l, guidefontsize=7)

#SELECT BEST MODELS FOR PARTITIONS
BGHMM_dict = Dict()
for (jobid, likelihood) in hmm_likelihoods_dict
    partition = jobid[1]
    order = jobid[3]
    if !haskey(BGHMM_dict, partition)
        BGHMM_dict[partition] = (hmm_results_dict[jobid][end][2], order, likelihood)
    else
        if likelihood > BGHMM_dict[partition][2]
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
