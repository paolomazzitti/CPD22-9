using Documenter, LARgenerators, DocumenterMarkdown
using DocumenterTools: Themes

makedocs(
    format = Documenter.HTML(
		prettyurls = get(ENV, "CI", nothing) == "true"
	),
    sitename="LARGenerators.jl",
    pages=[
        "Informazioni generali" => "index.md",
        "Introduzione" => "intro.md",
        "Grafo delle dipendenze" => "grafo.md",
        "Sviluppo" => "sviluppo.md",
        "Documentazione" => [
            "Documentazione (2D)" => "docs2d.md",
            "Documentazione (3D)" => "docs3d.md"
        ],
        "Analisi delle prestazioni" => "analisi.md",
    ],
    modules=[LARgenerators]
)

deploydocs(
    repo="github.com/paolomazzitti/CPD22-9.git"
)

Themes.compile("docs/src/documenter-dark.scss", "docs/build/assets/themes/documenter-dark.css")
Themes.compile("docs/src/documenter-light.scss", "docs/build/assets/themes/documenter-light.css")