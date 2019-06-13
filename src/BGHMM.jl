"""
Utility functions for learning and using background genomic hidden markov models
"""
module BGHMM
    using BioSequences, DataFrames, GenomicFeatures, ProgressMeter
    import Distributed: RemoteChannel
    import Distributions: Dirichlet, Categorical
    import MS_HMMBase: HMM, obs_set_likelihood
    import ProgressMeter: AbstractProgress
    import Printf: @sprintf
    import Random: rand
    import Statistics: mean

    #function to split random sample dataframe into training and test sets (divide total sequence length by half)
    function split_obs_sets(sample_dfs::Dict{String,DataFrame})
        training_sets = Dict{String,Vector{BioSequence{DNAAlphabet{4}}}}()
        test_sets = Dict{String,Vector{BioSequence{DNAAlphabet{4}}}}()

        for (partition_id, partition) in sample_dfs
            partition.sampleLength = (partition.SampleEnd - partition.SampleStart) .+ 1
            midway = sum(partition.sampleLength)÷2
            split_index = 0
            counter = 0
            while split_index == 0
                counter += 1
                length_sum = sum(partition.sampleLength[1:counter])
                if length_sum > midway
                    split_index = counter
                end
            end

            training_sets[partition_id]  = partition.SampleSequence[1:split_index-1]
            test_sets[partition_id] = partition.SampleSequence[split_index:end]
        end
        return training_sets, test_sets
    end

    #function to construct HMM transition matrix with strong priors on auto-transition
    function generate_transition_matrix(states::Int64, prior_dope::Float64=(states*250.0), prior_background::Float64=.1)
        transition_matrix=zeros(states,states)
        for k in 1:states
            dirichlet_params = fill(prior_background, states)
            dirichlet_params[k] = prior_dope
            transition_matrix[k,:] = rand(Dirichlet(dirichlet_params))
        end
        return transition_matrix
    end

    #function to construct HMM state emission distribution from uninformative dirichlet over the alphabet size
    function generate_emission_dist(no_emission_symbols, prior=Dirichlet(ones(no_emission_symbols)/no_emission_symbols))
        return Categorical(rand(prior))
    end

    #function to setup an HMM results dictionary and RemoteChannel for learning jobs, given a vector of state #s, order_nos, replicates to train, the dictionary to fill, the RemoteChannel and the training sequences
    #resumes any existing non-converged chains, otherwise initialises hmms for new chains given provided constants
    function HMM_setup!(order_nos::Array{Int64}, Ks::Array{Int64}, replicates::Int64, hmm_results_dict::Dict, input_hmms::RemoteChannel, training_sets::Dict{String,Vector{BioSequence{DNAAlphabet{4}}}}, base_alphabet_size::Int64)
        no_input_hmms = length(Ks)*length(order_nos)*replicates*length(training_sets)
        code_dict = Dict{Tuple{String,Int64}, Array{Int64}}()

        @showprogress 1 "Encoding observations..." for order_no in order_nos, (partition_id, partition) in training_sets #build the appropriate sample sets once
            order_seqs = get_order_n_seqs(partition,order_no) #get the kmer sequences at the appropriate order
            coded_seqs = code_seqs(order_seqs) #numerically code the sequences in trainable format
            code_dict[(partition_id, order_no)] = coded_seqs
        end

        #extremely primitive cluster memory management: building the queue with this iterator mostly keeps one machine from getting all the 6th order exonic jobs
        @showprogress 1 "Setting up HMMs..." for i in 1:replicates, order_no in order_nos, K in Ks, (partition_id, partition) in training_sets #for each combination of order and mosaic state number to test for each partition, init HMMs for workers
            jobid = (partition_id, K, order_no, i)
            if haskey(hmm_results_dict, jobid) && length(hmm_results_dict[jobid]) > 0 #true if resuming from incomplete chain
                lastindex = length(hmm_results_dict[jobid])
                job_convergence = hmm_results_dict[jobid][lastindex][5]
                if !job_convergence #push the last hmm iterate for nonconverged chains to the input channel
                    iterate = hmm_results_dict[jobid][lastindex][1]
                    hmm =  hmm_results_dict[jobid][lastindex][2]
                    put!(input_hmms, (jobid, iterate, hmm, code_dict[(partition_id, order_no)]))
                else #skip any jobs that have converged from previous runs
                    no_input_hmms -= 1
                end

            else #initialise first HMM in chain
                iterate = 1
                π0 = rand(Dirichlet(ones(K)/K)) #uninformative prior on initial state probabilities
                π = generate_transition_matrix(K)
                no_emission_symbols = Int(base_alphabet_size^(order_no+1)) #alphabet size for the order
                emission_dists = [generate_emission_dist(no_emission_symbols) for i in 1:K]
                #generate the HMM with the appropriate transition matrix and emissions distributions
                hmm = HMM(π0, π, emission_dists)
                hmm_results_dict[jobid] = [] #initialise the relevant results array
                put!(input_hmms, (jobid, iterate, hmm, code_dict[(partition_id, order_no)]))
            end
        end

        return no_input_hmms
    end

    include("hmm_tests.jl")
    include("hmm_progressmeter.jl")
    include("order_coding.jl")
    include("partition_masker.jl")
    include("sequence_sampler.jl")

end #module

