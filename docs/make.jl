using Documenter, LARgenerators, DocumenterMarkdown
using DocumenterTools: Themes

makedocs(sitename="LARGenerators.jl",
    pages=[
        "Introduzione"=>"index.md",
        "Grafo delle dipendenze" => "grafo.md",
        "Sviluppo" => "sviluppo.md",
        "Documentazione" => "docs.md",
        "Analisi delle prestazioni" => "analisi.md",
    ],
    modules=[LARgenerators]
)

deploydocs(
    repo="github.com/paolomazzitti/CPD22-9.git",
)

Themes.compile("docs/src/documenter-light.scss", "docs/build/assets/themes/documenter-light.css")