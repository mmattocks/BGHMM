function print_synthetic_fasta(path::String,line_length::Int64=80)
	header=">CM002885.2.1 BGHMM Synthetic chromosome for tests\n"
	write(path, header, format_lines(generate_synthetic_seq(),line_length))
end

function print_synthetic_index(path::String)
	write(path, "CM002885.2	1000	5	80	81\n")
end

function print_synthetic_position(path::String)
	write(path, ">1 start 501 end 850 smt_pos 111 smt_value 205.11821 fuzziness_score 43.53698\n",generate_synthetic_seq()[501:800],"\n")
end

function print_synthetic_gff(path::String)
	gff_txt = "##gff-version \n
	##sequence-region   1 1 1000\n
	#!genome-build BGHMM Synthetic\n
	#!genome-version Synth1.0\n
	#!genome-date 2019-08\n
	#!genome-build-accession NCBI:FAKE\n
	#!genebuild-last-updated 2019-08\n
	1	Ensembl	chromosome	1	1000	.	.	.	ID=chromosome:1;Alias=CM002885.2,NC_000001.1\n
	###\n
	1	ensembl_havana	gene	501	1000	.	-	.	ID=gene:ENSDARG00000000001;Name=fake01;biotype=protein_coding;description=fake gene [Source:FAKESRC];gene_id=ENSDARG00000000001;logic_name=ensembl_havana_gene;version=1\n
	1	ensembl_havana	mRNA	501	1000	.	-	.	ID=transcript:ENSDART00000000001;Parent=gene:ENSDARG00000000001;Name=fake01-203;biotype=protein_coding;transcript_id=ENSDART00000000001;version=1\n
	1	ensembl_havana	three_prime_UTR	501	510	.	-	.	Parent=transcript:ENSDART00000000001\n
	1	ensembl_havana	exon	501	570	.	-	.	Parent=transcript:ENSDART00000000001;Name=ENSDARE00000000001;constitutive=0;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSDARE00000000001;rank=5;version=1\n
	1	ensembl_havana	CDS	511	570	.	-	0	ID=CDS:ENSDARP00000000001;Parent=transcript:ENSDART00000000001;protein_id=ENSDARP00000000001\n
	1	ensembl_havana	exon	601	660	.	-	.	Parent=transcript:ENSDART00000000001;Name=ENSDARE00000000002;constitutive=0;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSDARE00000000002;rank=4;version=1\n
	1	ensembl_havana	CDS	601	660	.	-	1	ID=CDS:ENSDARP00000000001;Parent=transcript:ENSDART00000000001;protein_id=ENSDARP00000000001\n
	1	ensembl_havana	exon	701	760	.	-	.	Parent=transcript:ENSDART00000000001;Name=ENSDARE00000000003;constitutive=0;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSDARE00000000003;rank=3;version=1\n
	1	ensembl_havana	CDS	701	760	.	-	0	ID=CDS:ENSDARP00000000001;Parent=transcript:ENSDART00000000001;protein_id=ENSDARP00000000001\n
	1	ensembl_havana	exon	801	860	.	-	.	Parent=transcript:ENSDART00000000001;Name=ENSDARE00000000004;constitutive=0;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSDARE00000000004;rank=2;version=1\n
	1	ensembl_havana	CDS	801	860	.	-	0	ID=CDS:ENSDARP00000000001;Parent=transcript:ENSDART00000000001;protein_id=ENSDARP00000000001\n
	1	ensembl_havana	exon	901	975	.	-	.	Parent=transcript:ENSDART00000000001;Name=ENSDARE00000000005;constitutive=0;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSDARE00000000005;rank=1;version=1\n
	1	ensembl_havana	CDS	901	960	.	-	1	ID=CDS:ENSDARP00000000001;Parent=transcript:ENSDART00000000001;protein_id=ENSDARP00000000001\n
	1	ensembl_havana	five_prime_UTR	961	975	.	-	.	Parent=transcript:ENSDART00000000001\n
	###\n"
	write(path,gff_txt)
end

function generate_synthetic_seq(gene_start::Int64=501, UTR3L::Int64=10, exon_length::Int64=60, intron_length::Int64=40, no_exons::Int64=5, UTR5L::Int64=25, OAL::Int64=1000, (iseq,pseq,eseq)::Tuple{String,String,String}=("AT","CG","CAT"); verbose::Bool=false)
	@assert length(iseq) == 2
	@assert length(pseq) == 2
	@assert length(eseq) == 3
	@assert mod(gene_start-1,2) == 0
	@assert mod(exon_length,3) == 0
	@assert mod(intron_length,2) == 0
	gene_length=gene_start+UTR3L+(no_exons*exon_length)+((no_exons-1)*intron_length)+UTR5L
	@assert gene_length <=OAL
	seq=""
	for i in 1:((gene_start-1) / 2)
		seq *= iseq
	end
	verbose && @info "Intergenic length $(length(seq))"
	for i in 1:UTR3L / 2
		seq *= pseq
	end
	verbose && @info "Added 3'UTR... $(length(seq))"
	for ex in 1:no_exons-1
		for i in 1:(exon_length/3)
			seq*=eseq
		end
		ex == 1 ? ilength = 30 : ilength = 40
		for i in 1:(ilength/2)
			seq*=pseq
		end
	end
	for i in 1:(exon_length/3)
		seq*=eseq
	end
	verbose && @info "Added exons and introns... $(length(seq))"
	for i in 1:UTR5L / 2
		seq *= pseq
	end
	verbose && @info "Added 5'UTR... $(length(seq))"

	for i in 1:((OAL-length(seq))/2)
		seq*=iseq
	end

	verbose && @info "Generated synthetic genome sequence of length $(length(seq))"

	return seq
end

  function format_lines(seq::String, line_length::Int64)
	a=join((SubString(seq,i,min(i+line_length-1,length(seq))) for i=1:line_length:length(seq)),'\n')
	return(a)
  end