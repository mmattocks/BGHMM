#FUNCTIONS TO ANALYSE LEARNT HMMS
    #returns likelihood of appropriately coded observation set given hmm
    function test_hmm(hmm::HMM, test_set, order)
        order_seqs = get_order_n_seqs(test_set,order) #get the kmer sequences at the appropriate order
        coded_seqs = code_seqs(order_seqs) #numerically code the sequences in trainable format
        return MS_HMMBase.obs_set_likelihood(hmm, coded_seqs)
    end

    function get_diagonal_array(hmm::HMM)
        k = length(hmm.π0)
        if k > 1
            diagonal = zeros(k)
            for i in 1:k
                diagonal[i] = hmm.π[i,i]
            end
            return diagonal
        else
            return [NaN]
        end
    end

    #function to simulate run lengths for vector of diagonal values
    function sim_run_lengths(diagonal_value::Array{Float64}, samples::Int64)
        if isnan(diagonal_value[1])
            return [NaN]
        else
            mean_run_lengths = zeros(length(diagonal_value))
            for (i, value) in enumerate(diagonal_value)
                runlengths = zeros(Int64, samples)
                for s in 1:samples
                    run = true
                    runlength = 0
                    while run
                        runlength += 1
                        if rand(1)[1] > value
                            run = false
                        end
                    end
                    runlengths[s] = runlength
                end
                mean_run_lengths[i] = mean(runlengths)
            end
            return mean_run_lengths
        end
    end
