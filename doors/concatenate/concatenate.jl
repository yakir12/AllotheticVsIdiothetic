using CSV, DataFrames, ProgressMeter

const target = "AllotheticVsIdiothetic"
const home = "/home/yakir/mnt/Data/Elin/"

todo = ["Project_AllotheticVsIdiothetic_indoors", "Project_AllotheticVsIdiothetic_outdoors_50cm", "Project_AllotheticVsIdiothetic_outdoors"]

dfs = DataFrame[]
for path in todo, file in ["runs", "calibs"]
    path_field = "$(file)_path"
    _df = CSV.read(joinpath(home, path, "$file.csv"), DataFrame, select = [path_field, "file"])
    if path_field ∈ names(_df)
        select!(_df, [path_field, "file"] => ByRow((p, f) -> joinpath(home, path, p, f)) => :fullfile, :file)
    else
        select!(_df, :file, :file => ByRow(f -> joinpath(home, path, f)) => :fullfile)
    end
    push!(dfs, _df)
end
df = vcat(dfs...)

unique!(df, :fullfile)
@assert allunique(df.file)


if isdir(joinpath(home, target))
    rm.(readdir(joinpath(home, target), join = true), force = true, recursive = true)
else
    mkdir(target)
end

@showprogress for (fullfile, file) in zip(df.fullfile, df.file)
    cp(fullfile, joinpath(home, target, file))
end

for file in ["runs", "calibs"]
    csv_sources = joinpath.(home, todo, "$file.csv")
    _dfs = DataFrame[]
    for csv_source in csv_sources
        _df = CSV.read(csv_source, DataFrame)
        _df.csv_source .= csv_source
        push!(_dfs, _df)
    end
    _df = vcat(_dfs..., cols = :union)
    if file == "runs"
        transform!(groupby(_df, [:run_id, :csv_source]), groupindices => :temp_run_id)
        rename!(select!(_df, Not(:run_id)), :temp_run_id => :run_id)
    end
    path_field = "$(file)_path"
    _df[!, path_field] .= missing
    CSV.write(joinpath(home, target, "$file.csv"), _df)
end
