using Documenter, LARgenerators

makedocs(sitename="LARGenerators.jl",
    pages=[
        "index.md",
        "Page title" => "foo.md",
    ],
    modules=[LARgenerators]
)
deploydocs(
    repo="github.com:paolomazzitti/CPD22-9.git",
)