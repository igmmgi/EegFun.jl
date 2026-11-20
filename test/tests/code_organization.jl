# Automated structural, include-tree, docstring, and public API checks across `src/`, `ext/`, `docs/`, and `test/`.

using Test
using EegFun

@testset "Code organization & public API invariants" begin
    pkg_root = pkgdir(EegFun)
    JS = Base.JuliaSyntax

    banner_types = """
    ######################################
    ######### Types & Constants ##########
    ######################################"""

    banner_user = """
    ######################################
    ########### User Functions ###########
    ######################################"""

    banner_internal = """
    ######################################
    ######### Internal Functions #########
    ######################################"""

    special_files = Set(["src/EegFun.jl", "src/public_api.jl", "src/precompile.jl"])
    distinct_return_detector_families = Set(["is_extreme_value", "is_flatline", "is_peak_to_peak", "is_step_value"])

    function _is_func_def_expr(x)
        x isa Expr || return false
        x.head === :function && return true
        return x.head === :(=) && !isempty(x.args) && x.args[1] isa Expr && x.args[1].head in (:call, :where, :(::))
    end

    function _unwrap_expr_with_doc(ex)
        doc_str = nothing
        cur = ex
        while cur isa Expr && cur.head === :macrocall && !isempty(cur.args)
            mname = cur.args[1]
            if mname === GlobalRef(Core, Symbol("@doc")) || mname === Symbol("@doc")
                doc_str = string(cur.args[end-1])
                cur = cur.args[end]
            elseif _is_func_def_expr(cur.args[end]) || (cur.args[end] isa Expr && cur.args[end].head === :macrocall)
                cur = cur.args[end]
            else
                break
            end
        end
        return cur, doc_str
    end

    _unwrap_expr(ex) = first(_unwrap_expr_with_doc(ex))

    function _func_sig_info(ex)
        cur, doc_str = _unwrap_expr_with_doc(ex)
        if cur isa Expr && cur.head === :macrocall && !isempty(cur.args) && occursin("kernel", string(cur.args[1]))
            for a in cur.args[2:end]
                if _is_func_def_expr(a)
                    cur = a
                    break
                end
            end
        end
        _is_func_def_expr(cur) || return nothing
        sig = cur.args[1]
        where_params = Dict{Symbol,Any}()
        while sig isa Expr && sig.head === :where
            for wp in sig.args[2:end]
                wp isa Expr && wp.head === :(<:) && wp.args[1] isa Symbol && (where_params[wp.args[1]] = wp.args[2])
            end
            sig = sig.args[1]
        end
        sig isa Expr && sig.head === :(::) && (sig = sig.args[1])
        (sig isa Expr && sig.head === :call) || return nothing

        callee = sig.args[1]
        qual_mod = nothing
        fname = nothing
        if callee isa Symbol
            fname = callee
        elseif callee isa Expr && callee.head === :. && length(callee.args) == 2 && callee.args[2] isa QuoteNode
            qual_mod = callee.args[1] isa Symbol ? callee.args[1] : Symbol(sprint(print, callee.args[1]))
            fname = callee.args[2].value
        elseif callee isa Expr && callee.head === :curly
            base = callee.args[1]
            if base isa Symbol
                fname = base
            elseif base isa Expr && base.head === :. && length(base.args) == 2 && base.args[2] isa QuoteNode
                qual_mod = base.args[1] isa Symbol ? base.args[1] : Symbol(sprint(print, base.args[1]))
                fname = base.args[2].value
            end
        elseif callee isa Expr && callee.head === :(::)
            t_target = callee.args[end]
            fname = t_target isa Symbol ? t_target : Symbol(sprint(print, t_target))
        end
        fname === nothing && return nothing

        pos_args = [a for a in sig.args[2:end] if !(a isa Expr && a.head === :parameters)]
        arg1_type = ""
        if !isempty(pos_args)
            a1 = pos_args[1]
            a1 isa Expr && a1.head === :kw && (a1 = a1.args[1])
            if a1 isa Expr && a1.head === :(::)
                tex = a1.args[end]
                tex isa Symbol && haskey(where_params, tex) && (tex = where_params[tex])
                arg1_type = replace(sprint(print, tex), " " => "")
            end
        end
        return (qual_mod = qual_mod, fname = fname, arg1_type = arg1_type, doc_str = doc_str)
    end

    function _collect_syms!(s::Set{Symbol}, x)
        if x isa Symbol
            push!(s, x)
        elseif x isa Expr
            for a in x.args
                _collect_syms!(s, a)
            end
        end
        return s
    end

    function _collect_includes!(acc::Vector{String}, abs_file::String)
        parsed = Meta.parseall(read(abs_file, String))
        function _walk(x, cur_dir::String)
            x isa Expr || return
            if x.head === :call && length(x.args) == 2 && x.args[1] === :include && x.args[2] isa String
                child = normpath(joinpath(cur_dir, x.args[2]))
                push!(acc, child)
                isfile(child) && _walk(Meta.parseall(read(child, String)), dirname(child))
            else
                for a in x.args
                    _walk(a, cur_dir)
                end
            end
        end
        _walk(parsed, dirname(abs_file))
        return acc
    end

    is_str_arg(t) = t in (
        "String",
        "AbstractString",
        "Vector{String}",
        "AbstractVector{String}",
        "AbstractVector{<:AbstractString}",
        "Union{String,Vector{String}}",
    )
    is_meta_arg(t) = t == "Function" || startswith(t, "Type{") || t == "DataType" || t == "Type"

    pub_path = joinpath(pkg_root, "src", "public_api.jl")
    pub_syms = Set{Symbol}()
    current_cat = String[]
    cat_order_violations = String[]

    for line in eachline(pub_path)
        s = strip(line)
        if startswith(s, "# ===")
            empty!(current_cat)
        elseif startswith(s, "public ")
            sym_str = strip(s[8:end])
            push!(pub_syms, Symbol(sym_str))
            if !isempty(current_cat) && lowercase(sym_str) < lowercase(current_cat[end])
                push!(cat_order_violations, "$(current_cat[end]) > $sym_str")
            end
            push!(current_cat, sym_str)
        end
    end

    # 1. Public API <-> leading-underscore function invariant
    @testset "public_api.jl <-> underscore function invariant" begin
        @test isempty(cat_order_violations)
        undefined_pubs = [sym for sym in sort!(collect(pub_syms)) if !isdefined(EegFun, sym)]
        @test isempty(undefined_pubs)

        unlisted_non_underscore_funcs = String[]
        for sym in names(EegFun; all = true)
            str = string(sym)
            (startswith(str, "#") || startswith(str, "@") || startswith(str, "_")) && continue
            (sym in (:EegFun, :eval, :include) || startswith(str, "cpu__") || startswith(str, "gpu__")) && continue
            val = getfield(EegFun, sym)
            if val isa Function && parentmodule(val) === EegFun && !(sym in pub_syms)
                push!(unlisted_non_underscore_funcs, str)
            end
        end
        @test isempty(unlisted_non_underscore_funcs)
    end

    # 2. Subsystem include-tree completeness in src/EegFun.jl and test/runtests.jl
    @testset "src/ and test/ include-tree completeness" begin
        src_entry = normpath(joinpath(pkg_root, "src", "EegFun.jl"))
        included_src = _collect_includes!(String[src_entry], src_entry)
        all_src = String[]
        for (r, _, fs) in walkdir(joinpath(pkg_root, "src"))
            for f in fs
                endswith(f, ".jl") && push!(all_src, normpath(joinpath(r, f)))
            end
        end
        @test length(included_src) == length(unique(included_src))
        @test Set(included_src) == Set(all_src)

        runtests_path = normpath(joinpath(pkg_root, "test", "runtests.jl"))
        included_tests = _collect_includes!(String[], runtests_path)
        all_test_files = String[]
        for (r, _, fs) in walkdir(joinpath(pkg_root, "test", "tests"))
            for f in fs
                endswith(f, ".jl") && push!(all_test_files, normpath(joinpath(r, f)))
            end
        end
        @test Set(included_tests) == Set(all_test_files)
    end

    # 3. Intra-file headers, section banners, section ordering, dispatch ladder, and primary docstring alignment
    @testset "File headers, section banners, dispatch ladder, and docstring order" begin
        all_files = String[]
        for root_dir in ("src", "ext")
            abs_root = joinpath(pkg_root, root_dir)
            isdir(abs_root) || continue
            for (r, _, fs) in walkdir(abs_root)
                for f in fs
                    endswith(f, ".jl") && push!(all_files, relpath(joinpath(r, f), pkg_root))
                end
            end
        end
        sort!(all_files)

        header_violations = String[]
        rank_violations = String[]
        banner_violations = String[]
        ladder_violations = String[]
        primary_doc_violations = String[]

        for relf in all_files
            rel_norm = replace(relf, '\\' => '/')
            content = read(joinpath(pkg_root, relf), String)
            first_line = strip(first(split(content, '\n'; limit = 2)))
            if !startswith(first_line, "#") || startswith(first_line, "####") || startswith(first_line, "# fmt:")
                push!(header_violations, rel_norm)
            end

            rel_norm in special_files && continue

            tree = JS.parseall(JS.SyntaxNode, content; ignore_errors = false)
            children = JS.children(tree)
            if children !== nothing && length(children) == 1 && JS.kind(children[1]) == Base.JuliaSyntax.K"module"
                children = JS.children(JS.children(children[1])[2])
            end
            (children === nothing || isempty(children)) && continue

            exprs = [Expr(ch) for ch in children]
            # Pure include(...) hub files have no code sections or banners
            if all(ex -> (c = _unwrap_expr(ex); c isa Expr && c.head === :call && c.args[1] === :include), exprs)
                if occursin("######################################", content)
                    push!(banner_violations, "$rel_norm (unexpected banner in include hub)")
                end
                continue
            end

            local_internals = Set{Symbol}()
            imported_ext = Set{Symbol}()
            for ex in exprs
                cur = _unwrap_expr(ex)
                if cur isa Expr && cur.head === :import
                    for path in cur.args
                        if path isa Expr && path.head === :(:)
                            for sn in path.args[2:end]
                                if sn isa Expr && sn.head === :. && !isempty(sn.args)
                                    sym = sn.args[end]
                                    sym isa Symbol && sym !== :fit_mass_lmm && push!(imported_ext, sym)
                                end
                            end
                        end
                    end
                end
                si = _func_sig_info(ex)
                si !== nothing && si.qual_mod === nothing && startswith(string(si.fname), "_") && push!(local_internals, si.fname)
            end

            ranks = Int[]
            r2_sigs = NamedTuple[]
            for ex in exprs
                cur = _unwrap_expr(ex)
                rk = if cur isa Expr && (cur.head === :using || cur.head === :import)
                    0
                elseif cur isa Expr && cur.head === :call && cur.args[1] === :include
                    4
                elseif cur isa Expr && cur.head in (:struct, :abstract, :primitive, :macro)
                    1
                elseif cur isa Expr && cur.head === :const
                    syms = _collect_syms!(Set{Symbol}(), cur)
                    lhs = cur.args[1]
                    ls = lhs isa Symbol ? lhs : (lhs isa Expr && lhs.head === :(=) && lhs.args[1] isa Symbol ? lhs.args[1] : nothing)
                    ls !== nothing && delete!(syms, ls)
                    isempty(intersect(syms, local_internals)) ? 1 : 3
                elseif cur isa Expr &&
                       cur.head === :macrocall &&
                       !isempty(cur.args) &&
                       (
                           occursin("kwdef", string(cur.args[1])) ||
                           occursin("enum", string(cur.args[1])) ||
                           any(a -> a isa Expr && a.head === :struct, cur.args)
                       )
                    1
                else
                    si = _func_sig_info(ex)
                    if si === nothing
                        -1
                    elseif si.fname === :__init__ || (si.qual_mod === nothing && startswith(string(si.fname), "_"))
                        3
                    else
                        push!(r2_sigs, si)
                        2
                    end
                end
                push!(ranks, rk)
            end

            if any(==(-1), ranks) || !issorted(ranks)
                push!(rank_violations, rel_norm)
            end

            has_r1 = 1 in ranks
            has_r2 = 2 in ranks
            has_r3 = 3 in ranks
            expected_border_lines = 2 * (Int(has_r1) + Int(has_r2) + Int(has_r3))
            actual_border_lines = count(l -> strip(l) == "######################################", split(content, '\n'))
            if (has_r1 != occursin(banner_types, content)) ||
               (has_r2 != occursin(banner_user, content)) ||
               (has_r3 != occursin(banner_internal, content)) ||
               (count("######### Types & Constants ##########", content) != Int(has_r1)) ||
               (count("########### User Functions ###########", content) != Int(has_r2)) ||
               (count("######### Internal Functions #########", content) != Int(has_r3)) ||
               (actual_border_lines != expected_border_lines)
                push!(banner_violations, rel_norm)
            end

            ext_flags = [(si.qual_mod !== nothing && si.qual_mod !== :EegFun) || (si.fname in imported_ext) for si in r2_sigs]
            if !issorted(ext_flags)
                push!(ladder_violations, "$rel_norm (interface methods before primary functions)")
            else
                prim_sigs = [si for (si, is_ext) in zip(r2_sigs, ext_flags) if !is_ext]
                by_fam = Dict{String,Vector{NamedTuple}}()
                for si in prim_sigs
                    fam = rstrip(string(si.fname), '!')
                    push!(get!(() -> NamedTuple[], by_fam, fam), si)
                end
                for (fam, fsigs) in by_fam
                    if length(fsigs) > 1
                        has_inmem = any(si -> !is_str_arg(si.arg1_type) && !is_meta_arg(si.arg1_type), fsigs)
                        tiers = [(has_inmem && is_str_arg(si.arg1_type)) ? 3 : (endswith(string(si.fname), '!') ? 2 : 1) for si in fsigs]
                        if !issorted(tiers)
                            push!(ladder_violations, "$rel_norm (family `$fam`: tiers=$tiers)")
                        end
                    end

                    isuppercase(first(fam)) && continue
                    if any(si -> si.doc_str !== nothing, fsigs) && fsigs[1].doc_str === nothing
                        push!(primary_doc_violations, "$rel_norm (`$fam`: first method lacks docstring)")
                    end
                    if !(fam in distinct_return_detector_families)
                        nb_docs = [length(si.doc_str) for si in fsigs if string(si.fname) == fam && si.doc_str !== nothing]
                        b_docs = [length(si.doc_str) for si in fsigs if string(si.fname) == fam * "!" && si.doc_str !== nothing]
                        if !isempty(nb_docs) && !isempty(b_docs) && maximum(b_docs) > maximum(nb_docs) + 40
                            push!(primary_doc_violations, "$rel_norm (`$fam!` docstring longer than `$fam`)")
                        end
                    end
                end
            end
        end

        @test isempty(header_violations)
        @test isempty(rank_violations)
        @test isempty(banner_violations)
        @test isempty(ladder_violations)
        @test isempty(primary_doc_violations)

        test_header_violations = String[]
        standalone_test_violations = String[]
        for (r, _, fs) in walkdir(joinpath(pkg_root, "test"))
            for f in fs
                endswith(f, ".jl") || continue
                abs_p = joinpath(r, f)
                rel_norm = replace(relpath(abs_p, pkg_root), '\\' => '/')
                txt = read(abs_p, String)
                first_line = strip(first(split(txt, '\n'; limit = 2)))
                if !startswith(first_line, "#") || startswith(first_line, "####") || startswith(first_line, "# fmt:")
                    push!(test_header_violations, rel_norm)
                end
                if startswith(rel_norm, "test/tests/")
                    has_test = occursin(r"^using\s+.*\bTest\b"m, txt)
                    has_eegfun = occursin(r"^using\s+.*\bEegFun\b"m, txt)
                    if !has_test || !has_eegfun
                        push!(standalone_test_violations, "$rel_norm (Test=$has_test, EegFun=$has_eegfun)")
                    end
                end
            end
        end
        @test isempty(test_header_violations)
        @test isempty(standalone_test_violations)
    end

    # 4. Docstring binding integrity, canonical structure, & 100% public API documentation / reference-page coverage
    @testset "Docstring integrity, structure, & reference-page coverage" begin
        dmeta = Docs.meta(EegFun)

        canonical_doc_headings = Dict(
            "Fields" => 1,
            "Arguments" => 1,
            "Keyword Arguments" => 2,
            "Returns" => 3,
            "Details" => 4,
            "Notes" => 4,
            "Examples" => 5,
            "References" => 6,
            "See Also" => 7,
        )
        demoted_subheading_names =
            Set(["see also", "arguments", "keyword arguments", "returns", "examples", "example", "fields", "references"])

        detached_public_docs = String[]
        unindented_docstrings = String[]
        doc_heading_violations = String[]

        for (binding, mdoc) in dmeta
            sym_str = string(binding.var)
            for (_, docstr) in mdoc.docs
                raw_full = join(docstr.text, "")
                lines = split(rstrip(raw_full), '\n')
                first_idx = findfirst(l -> !isempty(strip(l)), lines)
                if first_idx !== nothing && !startswith(lines[first_idx], "    ")
                    push!(unindented_docstrings, "$sym_str ($(docstr.data[:path]):$(docstr.data[:linenumber]))")
                end

                in_fence = false
                headings = String[]
                bad_subheadings = String[]
                for l in lines
                    s = strip(l)
                    if startswith(s, "```")
                        in_fence = !in_fence
                        continue
                    end
                    in_fence && continue
                    if startswith(l, "# ")
                        push!(headings, strip(l[3:end]))
                    elseif startswith(l, "## ")
                        sh = strip(l[4:end])
                        lowercase(sh) in demoted_subheading_names && push!(bad_subheadings, sh)
                    end
                end
                unknown_headings = [h for h in headings if !haskey(canonical_doc_headings, h)]
                heading_ranks = [get(canonical_doc_headings, h, 99) for h in headings]
                if !isempty(unknown_headings) || !issorted(heading_ranks) || !isempty(bad_subheadings)
                    push!(
                        doc_heading_violations,
                        "$sym_str ($(docstr.data[:path]):$(docstr.data[:linenumber])): headings=$headings sub=$bad_subheadings",
                    )
                end

                if startswith(sym_str, "_")
                    first_sig = strip(first(split(strip(raw_full), '\n'; limit = 2)))
                    m = match(r"^([A-Za-z_][A-Za-z0-9_!]*)\s*(?:\(|\{|$)", first_sig)
                    if m !== nothing && !startswith(m.captures[1], "_")
                        push!(detached_public_docs, "$sym_str <- $(m.captures[1])")
                    end
                end
            end
        end
        @test isempty(detached_public_docs)
        @test isempty(unindented_docstrings)
        @test isempty(doc_heading_violations)

        undocumented_public_syms = [sym for sym in sort!(collect(pub_syms)) if !haskey(dmeta, Docs.Binding(EegFun, sym))]
        @test isempty(undocumented_public_syms)

        struct_fields_violations = String[]
        for sym in sort!(collect(pub_syms))
            haskey(dmeta, Docs.Binding(EegFun, sym)) || continue
            dt = Base.unwrap_unionall(getfield(EegFun, sym))
            if dt isa DataType && isstructtype(dt) && !isabstracttype(dt) && !isempty(fieldnames(dt))
                raw_doc = join([join(d.text, "") for (_, d) in dmeta[Docs.Binding(EegFun, sym)].docs], "\n")
                has_fields_sec = occursin("\n# Fields\n", raw_doc) || occursin("\n# Fields\r\n", raw_doc)
                missing_f = [f for f in fieldnames(dt) if !occursin("- `$(f)", raw_doc)]
                if !has_fields_sec || !isempty(missing_f)
                    push!(struct_fields_violations, "$sym (has_section=$has_fields_sec, missing=$missing_f)")
                end
            end
        end
        @test isempty(struct_fields_violations)

        # All public non-Function symbols (types and constants) must appear in docs/src/reference/types.md @docs
        types_md = joinpath(pkg_root, "docs", "src", "reference", "types.md")
        types_ref_syms = Set{Symbol}()
        in_docs = false
        for line in eachline(types_md)
            s = strip(line)
            if startswith(s, "```@docs")
                in_docs = true
            elseif in_docs && startswith(s, "```")
                in_docs = false
            elseif in_docs && startswith(s, "EegFun.")
                rest = split(s[8:end], "(")[1]
                push!(types_ref_syms, Symbol(strip(split(rest, "{")[1])))
            end
        end
        unlisted_public_types_consts = Symbol[]
        for sym in sort!(collect(pub_syms))
            val = getfield(EegFun, sym)
            if !(val isa Function) && !(sym in types_ref_syms)
                push!(unlisted_public_types_consts, sym)
            end
        end
        @test isempty(unlisted_public_types_consts)

        # Every src/ file defining a public Function docstring must appear once in reference/{analysis,plotting}.md Pages
        func_ref_pages = String[]
        for ref_file in ("analysis.md", "plotting.md")
            ref_txt = read(joinpath(pkg_root, "docs", "src", "reference", ref_file), String)
            for m in eachmatch(r"\"(src/[^\"]+\.jl)\"", ref_txt)
                push!(func_ref_pages, m.captures[1])
            end
        end
        @test length(func_ref_pages) == length(unique(func_ref_pages))
        missing_disk_pages = [p for p in func_ref_pages if !isfile(joinpath(pkg_root, p))]
        @test isempty(missing_disk_pages)

        public_func_files = Set{String}()
        for sym in pub_syms
            getfield(EegFun, sym) isa Function || continue
            haskey(dmeta, Docs.Binding(EegFun, sym)) || continue
            for (_, docstr) in dmeta[Docs.Binding(EegFun, sym)].docs
                push!(public_func_files, replace(relpath(docstr.data[:path], pkg_root), '\\' => '/'))
            end
        end
        @test issubset(public_func_files, Set(func_ref_pages))
    end

    # 5. Repository-wide code formatting compliance with .JuliaFormatter.toml
    @testset "JuliaFormatter compliance" begin
        jf_pkg = Base.identify_package("JuliaFormatter")
        if jf_pkg !== nothing
            JF = Base.require(jf_pkg)
            for dir in ("src", "ext", "test", "tutorials", "benchmark")
                abs_dir = joinpath(pkg_root, dir)
                isdir(abs_dir) && @test Base.invokelatest(JF.format, abs_dir; overwrite = false)
            end
        end
    end
end
