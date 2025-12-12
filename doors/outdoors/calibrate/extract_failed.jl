function get_file(fldr)
    for line in readlines(joinpath(fldr, "error.log"))
        m = match(r"FILE: (.+)", line)
        if isnothing(m)
            continue
        else
            fullfile = m.captures[1]
            _, file = splitdir(fullfile)
            return file
        end
    end
    return missing
end

files = get_file.(filter(startswith("jl_"), readdir(".")))
# resulted in 31 folders, blur was 3.
