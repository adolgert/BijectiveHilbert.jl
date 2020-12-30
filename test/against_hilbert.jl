using Vimes

dir = "/home/adolgert/dev/BijectiveHilbert.jl"

function initialise_noclean(dir)
    (isfile(joinpath(dir, "Project.toml")) && isdir(joinpath(dir, "src"))) ||
      error("No Julia project found at $dir")
    tmp = joinpath(tempdir(), "vimes-$(rand(UInt64))")
    mkdir(tmp)
    for path in readdir(dir)
        if !startswith(path, ".")
            cp(joinpath(dir, path), joinpath(tmp, path))
        end
    end
    return tmp
end
tmp = initialise_noclean(dir)
idx = Vimes.indices(joinpath(tmp, "src"), Vimes.defaults)
Vimes.mutate_and_reset(dir, tmp, idx)
Vimes.diff(dir, tmp)

idx1 = idx[1]
idx1[3]
Vimes.checktests(dir, tmp)
Vimes.diff(dir, tmp)
Vimes.mutate(tmp, idx)
aa="/tmp/vimes-2625159901050647083"
bb = "/tmp/vimes-2625159901050647083/src/hilbert.jl"

# Use the BijectiveHilbert.jl/test project.
using Test
using TestReports
mutation_idx = 2
Vimes.mutate_and_reset(dir, tmp, idx) do
    cd(tmp) do
        ts = @testset ReportingTestSet "" begin
            include("test/runtests.jl")
        end
        open("testlog$(mutation_idx).xml","w") do fh
            print(fh, report(ts))
        end
    end
end

using EzXML
function test_outcomes_of_xml(xmlpath)
    outcomes = Dict{String, Int}()
    for testsuite in findall("//testsuite", xmlpath)
        name = nothing
        failures = nothing
        for attrib in attributes(testsuite)
            if attrib.name == "name"
                name = attrib.content
            elseif attrib.name == "failures"
                failures = parse(Int, attrib.content)
            end
        end
        if !isnothing(name) && !isnothing(failures)
            outcomes[name] = failures
        end
    end
    outcomes
end

base="/tmp/vimes-10290787238000520708/"
log = EzXML.readxml(joinpath(base, "testlog1.xml"))
test_outcomes_of_xml(log)

ts1 = testsuites[1]
attribs = attributes(ts1)
attribs[1].name
attribs[1].content
ts1.failures
fieldnames(ts1)