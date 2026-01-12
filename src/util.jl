
macro import_everything(module_expr)
    @assert module_expr isa Symbol || Meta.isexpr(module_expr, :., 2)
    # The “bare” syntax “X.Y.Z” results in a different expression structure than
    # in an import expression, i.e. “import X.Y.Z”. Since both are needed, we
    # need to convert between the two. The input to the macro is the first form,
    # i.e. without “import”.
    import_args = []
    expr = module_expr
    while Meta.isexpr(expr, :.)
        @assert Meta.isexpr(expr, :., 2)
        last_qn = expr.args[end]
        @assert last_qn isa QuoteNode
        # eval to remove QuoteNode
        last = Core.eval(__module__, last_qn)
        @assert last isa Symbol
        pushfirst!(import_args, last)
        expr = expr.args[begin]
    end
    @assert expr isa Symbol
    pushfirst!(import_args, expr)
    return quote
        using $(import_args...): $(import_args...)
        for name in filter(
            n -> Base.isidentifier(n) && n ∉ (:eval, :include),
            names($(esc(module_expr)); all=true),
        )
            @eval import $(import_args...): $(Expr(:$, :name))
        end
    end
end

