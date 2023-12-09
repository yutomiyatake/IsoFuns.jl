using IsoFuns
using Documenter

DocMeta.setdocmeta!(IsoFuns, :DocTestSetup, :(using IsoFuns); recursive=true)

makedocs(;
    modules=[IsoFuns],
    authors="yutomiyatake <miyatake@cas.cmc.osaka-u.ac.jp> and contributors",
    # repo="https://github.com/yutomiyatake/IsoFuns.jl.git",
    sitename="IsoFuns.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yutomiyatake.github.io/IsoFuns.jl",
        edit_link="main",
        assets=String[],
    ),
    # pages=[
    #     "Home" => "index.md",
    # ],
    pages = ["Home" => "index.md", "Functions" => Any["func/iso.md", "func/neariso.md"],  "examples.md"]
)

deploydocs(;
    repo="github.com/yutomiyatake/IsoFuns.jl",
    devbranch="main",
)
