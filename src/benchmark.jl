
function time_all_methods(db_fasta_path, query_tsv_path, query_fasta_path;
                            receptor = "IG", uchime_mode = "balanced", CHMMAIRRa_bin = "CHMMAIRRa", disable_internal_dedup = true)
    println("Running CHMMAIRRa PC")
    CHMMAIRRa_PC = PC_time_chmmairra(db_fasta_path, query_tsv_path, receptor = receptor, detailed = false, disable_internal_dedup = disable_internal_dedup, CHMMAIRRa_bin = CHMMAIRRa_bin)
    println("Running CHMMAIRRa PC detailed")
    CHMMAIRRa_PC_detailed = PC_time_chmmairra(db_fasta_path, query_tsv_path, receptor = receptor, detailed = true, disable_internal_dedup = disable_internal_dedup, CHMMAIRRa_bin = CHMMAIRRa_bin)
    println("Running CHMMAIRRa")
    CHMMAIRRa = julia_time_chmmairra(db_fasta_path, query_tsv_path, receptor = receptor, detailed = false, disable_internal_dedup = disable_internal_dedup)
    println("Running CHMMAIRRa detailed")
    CHMMAIRRa_detailed = julia_time_chmmairra(db_fasta_path, query_tsv_path, receptor = receptor, detailed = true, disable_internal_dedup = disable_internal_dedup)
    println("Running USEARCH uchime2 ref")
    USEARCH_uchime2_ref = linux_time_usearch_uchime2_ref(query_fasta_path, db_fasta_path, mode = uchime_mode)
    println("Running VSEARCH uchime ref")
    VSEARCH_uchime_ref = linux_time_vsearch_uchime_ref(query_fasta_path, db_fasta_path)
    return NamedTuple{(:CHMMAIRRa, :CHMMAIRRa_detailed, :CHMMAIRRa_PC, :CHMMAIRRa_PC_detailed, :USEARCH_uchime2_ref, :VSEARCH_uchime_ref,)}((CHMMAIRRa, CHMMAIRRa_detailed, CHMMAIRRa_PC, CHMMAIRRa_PC_detailed, USEARCH_uchime2_ref, VSEARCH_uchime_ref,))
end

function julia_time_chmmairra(db_fasta_path::String, query_tsv_path::String;
                            receptor = "IG", detailed = false, disable_internal_dedup = true)
    return mktempdir() do dir
        t = @timed CHMMAIRRa.detect_chimeras_from_files(
            db_fasta_path, 
            query_tsv_path, 
            joinpath(dir, "out.tsv"),
            receptor = receptor,
            align_database = true,
            detailed = detailed,
            disable_internal_dedup = disable_internal_dedup
        )
        return t.time
    end
end

function PC_time_chmmairra(db_fasta_path::String, query_tsv_path::String; receptor = "IG", detailed = false, threads = Base.Threads.nthreads(), CHMMAIRRa_bin = "CHMMAIRRa", disable_internal_dedup = true)
    return mktempdir() do dir
        # weird command building block to prevent Julia's Cmd parsing from causing issues
        cmd = ["/usr/bin/time", "--format", "%e", CHMMAIRRa_bin, "--receptor", receptor, "--align-database", "--quiet"]
        if detailed
            push!(cmd, "--detailed")
        end
        if disable_internal_dedup
            push!(cmd, "--disable-internal-dedup")
        end
        # add positional args at the end
        cmd = vcat(cmd, ["--V-fasta", db_fasta_path, "--assignments", query_tsv_path, "--out", joinpath(dir, "out.tsv")])
        println(cmd)
        err = IOBuffer()
        ENV["JULIA_NUM_THREADS"] = threads
        run(pipeline(Cmd(cmd), stdout=devnull, stderr=err))
        return parse_linux_time_from_stderr(err)
    end
end

function linux_time_vsearch_uchime_ref(query_fasta_path::String, db_fasta_path::String, threads = Base.Threads.nthreads())
    return mktempdir() do dir
        cmd = `/usr/bin/time --format '%e' /home/mchernys/software/vsearch-2.22.1/bin/vsearch -uchime_ref $(query_fasta_path) -db $(db_fasta_path) -uchimeout $(joinpath(dir, "out.txt")) --threads $(threads)`
        println(cmd)
        err = IOBuffer()    
        run(pipeline(cmd, stdout=devnull, stderr=err))
        return parse_linux_time_from_stderr(err)
    end
end

function linux_time_usearch_uchime2_ref(query_fasta_path::String, db_fasta_path::String; mode = "balanced", threads = Base.Threads.nthreads())
    return mktempdir() do dir
        cmd = `/usr/bin/time --format '%e' usearch -uchime2_ref $(query_fasta_path) -db $(db_fasta_path) -uchimeout $(joinpath(dir, "out.txt")) -strand plus -mode $(mode) -threads $(threads)`
        println(cmd)
        err = IOBuffer()
        run(pipeline(cmd, stdout=devnull, stderr=err))
        return parse_linux_time_from_stderr(err)
    end
end

function parse_linux_time_from_stderr(err)
    return parse(Float64, split(strip(String(take!(err)), '\n'))[end])
end

function precompile_chmmairra(db_fasta_path::String, assignments_path::String)
    println("Precompiling CHMMAIRRa")    
    julia_time_chmmairra(db_fasta_path, assignments_path, receptor = "IG", detailed = false)
    julia_time_chmmairra(db_fasta_path, assignments_path, receptor = "IG", detailed = true)
    julia_time_chmmairra(db_fasta_path, assignments_path, receptor = "TCR", detailed = false)
    julia_time_chmmairra(db_fasta_path, assignments_path, receptor = "TCR", detailed = true)
end