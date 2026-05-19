using Pkg
cd(@__DIR__)
pkg" activate ."
pkg" dev Documenter BifurcationKit"

using Documenter, BifurcationKit, MultiParamContinuation

# using DocThemeIndigo
ENV["GKSwstype"] = "100"

# to display progress
ENV["JULIA_DEBUG"] = Documenter

makedocs(
	modules = [MultiParamContinuation],
	doctest = false,
	pagesonly = false, # this is on Documenter#master, do not compile what is not in pages =
	draft = false,
	warnonly = true,
	sitename = "Multi parameter continuation in Julia",
	authors = "Romain Veltz",
	format = Documenter.HTML(
		collapselevel = 1,
		assets=[
			asset("https://bifurcationkit.github.io/assets/js/documentation.js"),
			asset("https://bifurcationkit.github.io/assets/css/documentation.css"),
				],
	),
	pages = Any[
		"Home" => "index.md",
		"Tutorials" => "tutorials/tutorials.md",
		"Covering methods" => [
			"Henderson" => "henderson.md",
				],
		"Problems" => [
			"Manifold Problems" => "BifProblem.md",
			"BifurcationKit" => "BifProblemBK.md",
		],
		"Library" => "library.md"
	]
	)

deploydocs(;
	repo = "github.com/bifurcationkit/MultiParamContinuation.jl.git",
	push_preview=true, target="build", devbranch="main")
