# using StaticArrays, LaTeXStrings

# ---------------------------------------------------------------------------
# ParamSpec
# ---------------------------------------------------------------------------
# `path` addresses a *nested* model: the tuple of field names to walk from the
# top-level model down to the struct that actually owns `field`. Empty tuple
# means the field lives on the model itself.
#
#   ParamSpec(:H, L"H", 0:0.1:10; group=:vane, path=(:vane,))
#     → slider drives  m.vane.H
#
# `get`/`set` remain as an escape hatch for derived params with no backing field.

struct ParamSpec
    field::Symbol
    label
    range
    index::Union{Nothing,Int}
    group::Symbol
    path::Tuple{Vararg{Symbol}}
    get::Union{Nothing,Function}   # m -> value
    set::Union{Nothing,Function}   # (m, v) -> (possibly new) m
end

# positional form 
ParamSpec(field::Symbol, label, range; group, index=nothing, path=(), get=nothing, set=nothing) =
    ParamSpec(field, label, range, index, group, Tuple(path), get, set)

# lens form: no backing field, addresses nested/derived params
ParamSpec(label, range; group, get, set) =
    ParamSpec(:_, label, range, nothing, group, (), get, set)

"""
    resolve(model, s::ParamSpec)

Walk `s.path` from `model` to the struct that owns `s.field`. Returns `model`
itself when the path is empty.
"""
resolve(model, s::ParamSpec) = foldl(getproperty, s.path; init=model)

"""
    nest(specs, prefix::Symbol...; regroup=nothing)

Re-root a set of specs beneath `prefix`, so specs written against an inner model
address the same fields from an outer one.

`regroup` disambiguates group names when two nested models expose the same one
(e.g. two flagella both with a `:body` group). Pass a `Symbol` to force every
group, or a function `Symbol -> Symbol` to rewrite them.
"""
function nest(specs, prefix::Symbol...; regroup=nothing)
    _g = regroup === nothing ? identity :
         regroup isa Symbol  ? (_ -> regroup) : regroup
    [ParamSpec(s.field, s.label, s.range;
               group = _g(s.group),
               index = s.index,
               path  = (prefix..., s.path...),
               get   = s.get,
               set   = s.set)
     for s in specs]
end

specget(m, s::ParamSpec) = s.get !== nothing ? s.get(m) : begin
    parent = resolve(m, s)
    s.index === nothing ? getproperty(parent, s.field) :
                          getproperty(parent, s.field)[s.index]
end

function specset!(m, s::ParamSpec, v)
    s.set !== nothing && return s.set(m, v)
    parent = resolve(m, s)          # `Vane`/`PlanarVanedFlagellum` are mutable,
    if s.index === nothing          # so plain assignment through the path works
        setproperty!(parent, s.field, v)
    else
        setproperty!(parent, s.field, _setcomp(getproperty(parent, s.field), s.index, v))
    end
    m
end

# rebuild an SVector with component i replaced (SVector is immutable)
_setcomp(v::SVector{N,T}, i, x) where {N,T} =
    SVector{N,T}(ntuple(k -> k == i ? T(x) : v[k], N))
_setcomp(v::AbstractVector, i, x) = (w = copy(v); w[i] = x; w)   # fallback for Vector

# expand any vector-valued field into one spec per component
function expand_specs(m, specs)
    out = ParamSpec[]
    for s in specs
        if s.get !== nothing        # lens specs are opaque — nothing to expand
            push!(out, s); continue
        end
        v = getproperty(resolve(m, s), s.field)   # resolve *before* testing
        if v isa AbstractVector
            for k in eachindex(v)
                push!(out, ParamSpec(s.field, L"%$(s.field)^%$(k)", s.range;
                                     group=s.group, index=k, path=s.path))
            end
        else
            push!(out, s)
        end
    end
    out
end

# ---------------------------------------------------------------------------
# Specs per model
# ---------------------------------------------------------------------------

param_specs(::Type{<:PlanarFlagellum}) = [
    ParamSpec(:L, L"L", 0:0.1:20; group=:flagellum),
    ParamSpec(:C, L"C", -3:0.05:3; group=:flagellum),
    ParamSpec(:R₀, L"R_0", 0:0.05:1.5; group=:flagellum),
    ParamSpec(:R₁, L"R_1", 0:0.05:1.5; group=:flagellum),
    ParamSpec(:k, L"k", 0:0.05:8; group=:flagellum),
    ParamSpec(:ϕ, L"\phi", 0:0.05:8; group=:flagellum),
    ParamSpec(:ω, L"\omega", 0:0.5:100; group=:flagellum),
]

param_specs(::Type{<:ThreeDimensionalFlagellum}) = [
    ParamSpec(:L,   L"L",        0:0.1:20;    group=:body),

    ParamSpec(:fᵩ,  L"f_\phi",   0:0.5:100;   group=:phi),
    ParamSpec(:Aᵩ,  L"A_\phi",   0:0.1:1.5;   group=:phi),
    ParamSpec(:δᵩ,  L"δ_\phi",   0:0.01:4;    group=:phi),
    ParamSpec(:λᵩ,  L"λ_\phi",   0:0.5:30;    group=:phi),
    ParamSpec(:Cᵩ,  L"C_\phi",   -1:0.01:1;   group=:phi),
    ParamSpec(:Δγ,  L"Δϕ",       -1.0π:0.1:π; group=:phi),   # label ≠ field name

    ParamSpec(:f_θ, L"f_θ",      0:0.5:100;   group=:theta),
    ParamSpec(:A_θ, L"A_θ",      0:0.1:1.5;   group=:theta),
    ParamSpec(:δ_θ, L"δ_θ",      0:0.01:4;    group=:theta),
    ParamSpec(:λ_θ, L"λ_θ",      0:0.5:30;    group=:theta),
    ParamSpec(:C_θ, L"C_θ",      -1:0.01:1;   group=:theta),
]

param_specs(::Type{<:PlanarStandingWaveFlagellum}) = [
    ParamSpec(:L, L"L", 0:0.1:20;  group=:body),
    ParamSpec(:ω, L"ω", 0:0.5:100; group=:body),
    ParamSpec(:C, L"C", -3:0.1:3;  group=:body),
    ParamSpec(:R, L"R", -1:0.05:1; group=:R),      # → R_{1..4} in one column
    ParamSpec(:I, L"I", -1:0.05:1; group=:I),      # → I_{1..4} in another
]

param_specs(::Type{<:ThreeDimensionalStandingWaveFlagellum}) = [
    ParamSpec(:L, L"L", 0:0.1:20;  group=:body),
    ParamSpec(:ω, L"ω", 0:0.5:100; group=:body),

    ParamSpec(:C_θ, L"C_θ", -3:0.1:3;  group=:theta),
    ParamSpec(:R_θ, L"R_θ", -1:0.05:1; group=:theta),
    ParamSpec(:I_θ, L"I_θ", -1:0.05:1; group=:theta),
    ParamSpec(:Cᵩ,  L"Cᵩ",  -3:0.1:3;  group=:phi),
    ParamSpec(:Rᵩ,  L"Rᵩ",  -1:0.05:1; group=:phi),
    ParamSpec(:Iᵩ,  L"Iᵩ",  -1:0.05:1; group=:phi),
]

# `direction` is deliberately omitted: "Planar" pins it to -ez, and as an
# SVector{3} it would expand to three sliders that break planarity.
param_specs(::Type{<:Vane}) = [
    ParamSpec(:s_start, L"s_0", 0:0.01:1;  group=:vane),
    ParamSpec(:s_end,   L"s_1", 0:0.01:1;  group=:vane),
    ParamSpec(:H,       L"H",   0:0.1:10;  group=:vane),
]

# Recurse on the type parameter — no per-inner-model methods needed.
param_specs(::Type{<:PlanarVanedFlagellum{FM}}) where {FM} = vcat(
    nest(param_specs(FM), :flagellum),
    nest(param_specs(Vane), :vane),
)

# ---------------------------------------------------------------------------
# design
# ---------------------------------------------------------------------------

"""
    _sliders!(fig, row, m, specs, redraw)

Build one `SliderGrid` per group, wiring each slider to `specset!` + `redraw`.
"""
function _sliders!(fig, row, m, specs, redraw)
    for (col, g) in enumerate(unique(s.group for s in specs))
        gspecs = filter(s -> s.group == g, specs)
        rows = [(label = s.label, range = s.range, startvalue = specget(m, s))
                for s in gspecs]
        sg = SliderGrid(fig[row, col], rows...)
        for (slider, s) in zip(sg.sliders, gspecs)
            on(slider.value) do val
                specset!(m, s, val)
                redraw()
            end
        end
    end
end

function _controls!(fig, animstep, dt)
    run = Button(fig[2,1]; label="Start/Pause", tellwidth=false)
    isrunning = Observable(false)
    on(run.clicks) do _
        isrunning[] = !isrunning[]
        isrunning[] || return
        @async while isrunning[]
            isopen(fig.scene) || break
            animstep()
            sleep(dt)
        end
    end
    fig
end

"""
    design(m::Model; fps=30, specs=param_specs(typeof(m)), limits=..., kwargs...)

Interactive parameter tuner. Builds the visualisation straight from the model —
no `Part`, no `NearestDiscretisation`. Extra `kwargs` are forwarded to `viz!`
(and thence `get_buffer`), so e.g. `design(vf; N=80, N_v=25, N_h=12)`.
"""
function design(m::Model; fps=30, specs=param_specs(typeof(m)),
                limits=(nothing, nothing, nothing, nothing, nothing, nothing),
                kwargs...)
    specs = expand_specs(m, specs)          # ← vectors become component sliders
    fig = Figure()
    ax  = Axis3(fig[1:2,1:3], aspect=:data, limits=limits)
    B   = viz!(ax, m; t=0.0, kwargs...)

    dt = 1/fps; t = Ref(0.0)
    redraw()   = update_buffer_observable!(B, m, t[])
    animstep() = (t[] += dt; redraw())

    _sliders!(fig, 3, m, specs, redraw)
    display(fig)
    _controls!(fig, animstep, dt)
    return fig
end