using Accessors, StaticArrays

struct ParamSpec
    field::Symbol
    label
    range
    index::Union{Nothing,Int}
    group::Symbol
    get::Union{Nothing,Function}   # m -> value
    set::Union{Nothing,Function}   # (m, v) -> (possibly new) m
end

# your existing positional form — unchanged behaviour
ParamSpec(field::Symbol, label, range; group, index=nothing, get=nothing, set=nothing) =
    ParamSpec(field, label, range, index, group, get, set)

# lens form: no backing field, addresses nested/derived params
ParamSpec(label, range; group, get, set) =
    ParamSpec(:_, label, range, nothing, group, get, set)

specget(m, s::ParamSpec) = s.get !== nothing ? s.get(m) :
    (s.index === nothing ? getproperty(m, s.field) : getproperty(m, s.field)[s.index])

function specset!(m, s::ParamSpec, v)
    s.set !== nothing && return s.set(m, v)
    if s.index === nothing
        setproperty!(m, s.field, v)
    else
        setproperty!(m, s.field, _setcomp(getproperty(m, s.field), s.index, v))
    end
    m
end


# struct ParamSpec
#     field::Symbol
#     label::AbstractString
#     range::AbstractRange
#     group::Symbol
#     index::Union{Int,Nothing}    # nothing = scalar field; k = kth component
# end

# ParamSpec(field, label, range; group=:params, index=nothing) =
#     ParamSpec(field, label, range, group, index)

# rebuild an SVector with component i replaced (SVector is immutable)
_setcomp(v::SVector{N,T}, i, x) where {N,T} =
    SVector{N,T}(ntuple(k -> k == i ? T(x) : v[k], N))
_setcomp(v::AbstractVector, i, x) = (w = copy(v); w[i] = x; w)   # fallback for Vector

# expand any vector-valued field into one spec per component
function expand_specs(tdf, specs)
    out = ParamSpec[]
    for s in specs
        v = getproperty(tdf, s.field)
        if v isa AbstractVector
            for k in eachindex(v)
                # lbl = latexstring(string(s.field), "_{", k, "}")   # R_{1}, R_{2}, …
                push!(out, ParamSpec(s.field, L"%$(s.field)^%$(k)", s.range; group=s.group, index=k))
            end
        else
            push!(out, s)
        end
    end
    out
end

param_specs(::Type{<:ThreeDimensionalFlagellum}) = [
    ParamSpec(:L,   L"L",        0:0.1:20;    group=:body),

    ParamSpec(:fᵩ,  L"f_\phi",   0:0.5:100;   group=:phi),
    ParamSpec(:Aᵩ,  L"A_\phi",   0:0.1:1.5;   group=:phi),
    ParamSpec(:δᵩ,  L"δ_\phi",   0:0.01:4;    group=:phi),
    ParamSpec(:λᵩ,  L"λ_\phi",   0:0.5:30;    group=:phi),
    ParamSpec(:Cᵩ,  L"C_\phi",   -3:0.05:3;    group=:phi),
    ParamSpec(:Δγ,  L"Δϕ",       -1.0π:0.1:π; group=:phi),   # label ≠ field name

    ParamSpec(:f_θ, L"f_θ",      0:0.5:100;   group=:theta),
    ParamSpec(:A_θ, L"A_θ",      0:0.1:1.5;   group=:theta),
    ParamSpec(:δ_θ, L"δ_θ",      0:0.01:4;    group=:theta),
    ParamSpec(:λ_θ, L"λ_θ",      0:0.5:30;    group=:theta),
    ParamSpec(:C_θ, L"C_θ",      -3:0.05:3;    group=:theta),
]



param_specs(::Type{<:PlanarStandingWaveFlagellum}) = [
    ParamSpec(:L, L"L", 0:0.1:20;  group=:body),
    ParamSpec(:ω, L"ω", 0:0.5:100; group=:body),
    ParamSpec(:C, L"C", -3:0.1:3;  group=:body),
    ParamSpec(:R, L"R", -1:0.05:1;  group=:R),      # → R_{1..4} in one column
    ParamSpec(:I, L"I", -1:0.05:1;  group=:I),      # → I_{1..4} in another
]


param_specs(::Type{<:ThreeDimensionalStandingWaveFlagellum}) = [
    ParamSpec(:L, L"L", 0:0.1:20; group=:body),
    ParamSpec(:ω, L"ω", 0:0.5:100; group=:body),

    ParamSpec(:C_θ, L"C_θ", -3:0.1:3; group=:theta),
    ParamSpec(:R_θ, L"R_θ", -1:0.05:1; group=:theta),
    ParamSpec(:I_θ, L"I_θ", -1:0.05:1; group=:theta),
    ParamSpec(:Cᵩ, L"Cᵩ", -3:0.1:3; group=:phi),
    ParamSpec(:Rᵩ, L"Rᵩ", -1:0.05:1; group=:phi),
    ParamSpec(:Iᵩ, L"Iᵩ", -1:0.05:1; group=:phi)
]

param_specs(::Type{<:Vane}) = [
    ParamSpec(:s_start, L"s_0", 0:0.01:1;  group=:vane),
    ParamSpec(:s_end,   L"s_1", 0:0.01:1;  group=:vane),
    ParamSpec(:H,       L"H",   0:0.1:10;  group=:vane),
]

param_specs(::Type{<:PlanarVanedFlagellum{FM}}) where {FM} = vcat(
    nest(param_specs(FM), :flagellum),
    nest(param_specs(Vane), :vane),
)


function design(tdf; fps=30, specs=param_specs(typeof(tdf)))
    specs = expand_specs(tdf, specs)          # ← vectors become component sliders
    fig = Figure()
    ax = Axis3(fig[1:2,1:3], aspect=:data, limits=(-5.,5.,-5.,5.,-5.,5.))
    p = Part(tdf, 33, 117)
    B = viz!(ax, p)
    dt = 1/fps; t = Ref(0.0)
    animstep() = (t[] += dt; update_boundary!(p, t[]); update_buffer_observable!(B, p))

    for (col, g) in enumerate(unique(s.group for s in specs))
        gspecs = filter(s -> s.group == g, specs)
        rows = [(label = s.label, range = s.range,
                 startvalue = s.index === nothing ? getproperty(tdf, s.field)
                                                  : getproperty(tdf, s.field)[s.index])
                for s in gspecs]
        sg = SliderGrid(fig[3, col], rows...)

        for (slider, s) in zip(sg.sliders, gspecs)
            field, idx = s.field, s.index
            on(slider.value) do val
                if idx === nothing
                    setproperty!(tdf, field, val)
                else
                    setproperty!(tdf, field, _setcomp(getproperty(tdf, field), idx, val))
                end
                update_boundary!(p, t[]); update_buffer_observable!(B, p)
            end
        end
    end

    display(fig)
    run = Button(fig[2,1]; label="Start/Pause", tellwidth=false)
    isrunning = Observable(false)
    on(run.clicks) do _; isrunning[] = !isrunning[]; end
    on(run.clicks) do _
        @async while isrunning[]
            isopen(fig.scene) || break; animstep(); sleep(dt)
        end
    end
    return fig
end

function design(p::Part; fps=30, specs=param_specs(typeof(p.model)))
    specs = expand_specs(p.model, specs)          # ← vectors become component sliders
    fig = Figure()
    ax = Axis3(fig[1:2,1:3], aspect=:data, limits=(-5.,5.,-5.,5.,-5.,5.))
    B = viz!(ax, p)
    dt = 1/fps; t = Ref(0.0)
    animstep() = (t[] += dt; update_boundary!(p, t[]); update_buffer_observable!(B, p))

    for (col, g) in enumerate(unique(s.group for s in specs))
        gspecs = filter(s -> s.group == g, specs)
        rows = [(label = s.label, range = s.range,
                 startvalue = s.index === nothing ? getproperty(p.model, s.field)
                                                  : getproperty(p.model, s.field)[s.index])
                for s in gspecs]
        sg = SliderGrid(fig[3, col], rows...)

        for (slider, s) in zip(sg.sliders, gspecs)
            field, idx = s.field, s.index
            on(slider.value) do val
                if idx === nothing
                    setproperty!(p.model, field, val)
                else
                    setproperty!(p.model, field, _setcomp(getproperty(p.model, field), idx, val))
                end
                update_boundary!(p, t[]); update_buffer_observable!(B, p)
            end
        end
    end

    display(fig)
    run = Button(fig[2,1]; label="Start/Pause", tellwidth=false)
    isrunning = Observable(false)
    on(run.clicks) do _; isrunning[] = !isrunning[]; end
    on(run.clicks) do _
        @async while isrunning[]
            isopen(fig.scene) || break; animstep(); sleep(dt)
        end
    end
    return fig
end